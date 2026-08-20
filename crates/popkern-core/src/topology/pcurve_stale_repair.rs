//! Stale-pcurve detection and improve-or-revert repair.
//!
//! Nested inside `pcurve_recompute` (shares its module context);
//! extracted to hold the parent under its frozen file-size ratchet.

use crate::topology::BRep;

impl BRep {
    /// Detect and rebuild pcurves that no longer agree with the 3D model.
    ///
    /// For every half-edge that stores a pcurve, maps the pcurve's
    /// parameter-range endpoints through the owning face's surface and
    /// compares them against the half-edge's actual start/end vertex
    /// positions, and also measures the interior trace deviation between
    /// the pcurve image and the edge's 3D curve. If either error exceeds
    /// the stale tolerance (`max_imprecision` widened to at least the
    /// linear model tolerance), the pcurve is dropped as stale. All
    /// missing pcurves are then recomputed via `recompute_all_pcurves`,
    /// and `fix_pcurve_endpoints` is run per face so the rebuilt pcurves
    /// land exactly on their loop vertices.
    pub fn fix_stale_pcurves(&mut self) {
        let faces: Vec<_> = self.faces.keys().collect();
        self.fix_stale_pcurves_among(&faces);
    }

    /// As [`BRep::fix_stale_pcurves`], restricted to the half-edges
    /// of `faces` (PK-274): whole-model staleness repair at operator
    /// tails is priced by model size, while an operation only
    /// invalidates pcurves in the neighbourhood it rewired. Callers
    /// pass the dirty faces dilated by one ring (see
    /// `BRep::validation_region`); the whole-model form remains the
    /// serving-boundary and debug cross-check path
    /// ([`BRep::stale_pcurve_handles`]). Same detection gate,
    /// recompute, endpoint fix, and improve-or-revert contract.
    pub fn fix_stale_pcurves_among(&mut self, faces: &[crate::topology::FaceHandle]) {
        // Improve-or-revert contract (PK-125): remember every cleared
        // pcurve with its measured staleness. Recomputation projects
        // walks through the surface chart, and near collapsed or
        // folded chart regions the walk can land on a WRONG path
        // whose trace is orders worse than the "stale" original
        // (measured: pyramid4_mixed sector-1 cap iso, 3.3e-6 over the
        // staleness gate, recomputed into a trace the downstream
        // boundary-row relock then chased into a genuine surface
        // fold). A repair pass must never replace evidence with
        // something it can measure to be worse.
        let region: Vec<crate::topology::FaceHandle> = {
            let mut seen = std::collections::BTreeSet::new();
            faces
                .iter()
                .copied()
                .filter(|fh| self.faces.contains_key(*fh) && seen.insert(*fh))
                .collect()
        };
        let he_keys = self.region_half_edges(&region);
        let cleared = self.stale_pcurve_scan(&he_keys);
        for (heh, _, _) in &cleared {
            self.half_edges[*heh].pcurve = None;
        }

        let _ = self.recompute_pcurves_for(&region);

        // Endpoint fixes: compute per face in parallel (read-only,
        // face-independent — see `compute_pcurve_endpoint_fixes`),
        // apply serially in the unchanged region/loop order (PK-302).
        #[cfg(feature = "parallel")]
        let endpoint_fixes: Vec<Vec<_>> = {
            use rayon::prelude::*;
            let brep: &Self = self;
            region
                .par_iter()
                .map(|&fh| brep.compute_pcurve_endpoint_fixes(fh))
                .collect()
        };
        #[cfg(not(feature = "parallel"))]
        let endpoint_fixes: Vec<Vec<_>> = region
            .iter()
            .map(|&fh| self.compute_pcurve_endpoint_fixes(fh))
            .collect();
        for (heh, fixed) in endpoint_fixes.into_iter().flatten() {
            self.half_edges[heh].pcurve = Some(fixed);
        }

        // Improve-or-revert: the recomputed-deviation measurement is
        // read-only and per-half-edge (each reads only its own stored
        // pcurve, its face surface, and its edge curve), so it runs
        // in parallel with an ordered collect; the keep/snap/revert
        // decisions stay serial in the unchanged order (PK-302).
        let measure_recomputed = |heh: crate::topology::HalfEdgeHandle| -> Option<f64> {
            self.half_edges
                .get(heh)
                .filter(|half_edge| half_edge.pcurve.is_some())?;
            self.measure_half_edge_pcurve_deviation(heh)
                .map(|report| report.max_deviation)
        };
        #[cfg(feature = "parallel")]
        let recomputed_errors: Vec<Option<f64>> = {
            use rayon::prelude::*;
            cleared
                .par_iter()
                .map(|(heh, _, _)| measure_recomputed(*heh))
                .collect()
        };
        #[cfg(not(feature = "parallel"))]
        let recomputed_errors: Vec<Option<f64>> = cleared
            .iter()
            .map(|(heh, _, _)| measure_recomputed(*heh))
            .collect();
        for ((heh, original, original_error), recomputed_error) in
            cleared.into_iter().zip(recomputed_errors)
        {
            if !self.half_edges.contains_key(heh) {
                continue;
            }
            let keep_recomputed =
                recomputed_error.is_some_and(|error| error <= original_error.max(1.0e-12));
            if keep_recomputed {
                // A recomputed pcurve on a periodic surface can land a
                // whole period away from the face's chart window (the
                // projection's principal branch): the 3-D image is
                // identical, so no deviation measure can catch it, but
                // every consumer that reasons in the chart (boundary-
                // row relocks, trim fits, seam measures) is thrown a
                // full period off (PK-125: pyramid4_mixed support
                // pcurve recomputed at u - 2π; the export-repair
                // relock then chased the off-window trace into a
                // genuine cap-surface fold). Translate the recomputed
                // pcurve by whole periods back to the ORIGINAL's
                // branch — an exact, image-preserving operation.
                self.snap_pcurve_to_reference_branch(heh, &original);
            } else if let Some(half_edge) = self.half_edges.get_mut(heh) {
                half_edge.pcurve = Some(original);
            }
        }
    }

    /// Read-only staleness detection over the WHOLE model: every
    /// half-edge whose stored pcurve fails the same endpoint/trace
    /// gate [`BRep::fix_stale_pcurves`] repairs. This is the
    /// cross-check face of the scoped repair (PK-274): after a
    /// scoped tail pass, debug builds assert that no stale pcurve
    /// survives OUTSIDE the scoped region — a miss here means the
    /// dirty-region derivation lost a face, which would otherwise
    /// surface as a silent mesher/export defect.
    pub fn stale_pcurve_handles(&self) -> Vec<crate::topology::HalfEdgeHandle> {
        let he_keys: Vec<_> = self.half_edges.keys().collect();
        self.stale_pcurve_scan(&he_keys)
            .into_iter()
            .map(|(heh, _, _)| heh)
            .collect()
    }

    /// Every half-edge of every loop (outer + inner) of `faces`.
    fn region_half_edges(
        &self,
        faces: &[crate::topology::FaceHandle],
    ) -> Vec<crate::topology::HalfEdgeHandle> {
        let mut out = Vec::new();
        for &fh in faces {
            let Some(face) = self.faces.get(fh) else {
                continue;
            };
            let loops = std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied());
            for loop_handle in loops {
                let Some(loop_data) = self.loops.get(loop_handle) else {
                    continue;
                };
                let start = loop_data.half_edge;
                let mut cursor = start;
                for _ in 0..=self.half_edges.len() {
                    if !self.half_edges.contains_key(cursor) {
                        break;
                    }
                    out.push(cursor);
                    cursor = self.half_edges[cursor].next;
                    if cursor == start {
                        break;
                    }
                }
            }
        }
        out
    }

    /// The detection half of the stale repair: measure each
    /// half-edge's stored pcurve against its vertices (endpoint
    /// error) and its edge's 3-D curve (trace error); return the
    /// ones beyond the stale tolerance with their pcurves and
    /// measured errors. Read-only.
    fn stale_pcurve_scan(
        &self,
        half_edges: &[crate::topology::HalfEdgeHandle],
    ) -> Vec<(
        crate::topology::HalfEdgeHandle,
        crate::geometry::pcurve::Pcurve,
        f64,
    )> {
        // Deduplicate in first-appearance order, then measure each
        // half-edge independently. The measurement is read-only and
        // per-half-edge, so it runs in parallel under the `parallel`
        // feature; the indexed ordered collect keeps the result
        // sequence identical to the serial scan (PK-302).
        let mut seen = std::collections::BTreeSet::new();
        let deduped: Vec<crate::topology::HalfEdgeHandle> = half_edges
            .iter()
            .copied()
            .filter(|heh| seen.insert(*heh))
            .collect();
        let measure = |heh: crate::topology::HalfEdgeHandle| -> Option<(
            crate::topology::HalfEdgeHandle,
            crate::geometry::pcurve::Pcurve,
            f64,
        )> {
            use crate::geometry::surface::SurfaceEval;
            let half_edge = self.half_edges.get(heh)?;
            let pcurve = half_edge.pcurve.clone()?;

            let loop_handle = half_edge.loop_;
            let face_handle = self.loops.get(loop_handle)?.face;
            let surface = self
                .faces
                .get(face_handle)
                .and_then(|f| f.surface.as_ref())?
                .geometry
                .clone();

            let (t0, t1) = pcurve.parameter_range();
            let uv_start = pcurve.point_at(t0);
            let uv_end = pcurve.point_at(t1);
            let start_from_pcurve = surface.point_at(uv_start.x, uv_start.y);
            let end_from_pcurve = surface.point_at(uv_end.x, uv_end.y);
            let actual_start = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let actual_end = self.vertices[self.half_edges[next_heh].origin].position;

            let endpoint_error = start_from_pcurve
                .distance_to(&actual_start)
                .max(end_from_pcurve.distance_to(&actual_end));
            let trace_error = self
                .measure_half_edge_pcurve_deviation(heh)
                .map_or(0.0, |report| report.max_deviation);
            let stale_tolerance = self.max_imprecision.max(self.tolerance.linear());
            (endpoint_error > stale_tolerance || trace_error > stale_tolerance)
                .then(|| (heh, pcurve, endpoint_error.max(trace_error)))
        };
        #[cfg(feature = "parallel")]
        let measured: Vec<_> = {
            use rayon::prelude::*;
            deduped.par_iter().map(|&heh| measure(heh)).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let measured: Vec<_> = deduped.iter().map(|&heh| measure(heh)).collect();
        measured.into_iter().flatten().collect()
    }

    /// Translate `heh`'s pcurve by whole surface periods so its
    /// parameter-space midpoint lands on the same branch as
    /// `reference`. Image-preserving on periodic surfaces; a no-op
    /// everywhere else.
    fn snap_pcurve_to_reference_branch(
        &mut self,
        heh: crate::topology::HalfEdgeHandle,
        reference: &crate::geometry::pcurve::Pcurve,
    ) {
        let Some(half_edge) = self.half_edges.get(heh) else {
            return;
        };
        let Some(face) = self
            .loops
            .get(half_edge.loop_)
            .map(|loop_data| loop_data.face)
        else {
            return;
        };
        let Some(surface) = self.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            return;
        };
        let (period_u, period_v) = surface.geometry.periodicity();
        if period_u.is_none() && period_v.is_none() {
            return;
        }
        let Some(pcurve) = half_edge.pcurve.as_ref() else {
            return;
        };
        let midpoint_of = |pc: &crate::geometry::pcurve::Pcurve| {
            let (t0, t1) = pc.parameter_range();
            pc.point_at(0.5 * (t0 + t1))
        };
        let current = midpoint_of(pcurve);
        let target = midpoint_of(reference);
        let branch_shift = |delta: f64, period: Option<f64>| -> f64 {
            match period {
                Some(period) if period > 0.0 => (delta / period).round() * period,
                _ => 0.0,
            }
        };
        let shift_u = branch_shift(target.x - current.x, period_u);
        let shift_v = branch_shift(target.y - current.y, period_v);
        if shift_u == 0.0 && shift_v == 0.0 {
            return;
        }
        if let Some(half_edge) = self.half_edges.get_mut(heh)
            && let Some(pcurve) = half_edge.pcurve.as_mut()
        {
            *pcurve = pcurve.translated_uv(shift_u, shift_v);
        }
    }
}
