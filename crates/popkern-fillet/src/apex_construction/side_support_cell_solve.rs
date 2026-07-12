//! Support-cell surface solve for cylinder-mediated side-support bridges.
//!
//! The outside-outside side-support continuation closes the re-routed shell
//! with topology-only bridge faces. This pass consumes them: every maximal
//! run of bridge boundary that lies exactly on one local carrier surface
//! (support plane or fillet cylinder) is split off as its own cell via `mef`
//! with an exact carrier-contained closing curve, and the cell receives the
//! carrier surface plus recomputed pcurves. Boundary that fits no carrier
//! stays on the residual wrap face, which remains surface-less — the same
//! tolerated pocket-wrap class the passing direct-apex fixtures ship (see
//! docs/decisions/2026-07-11-side-support-pocket-floor-reconstruction.md).
//!
//! Every face is processed transactionally: any Euler-op failure, curve
//! authoring failure, or cell that fails its own on-carrier postcondition
//! rolls the face's edits back and records a decline reason. The scaffold
//! path then continues exactly as before, so fixtures that never enter this
//! route are untouched.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

#[path = "side_support_cell_solve_geometry.rs"]
mod geometry;
#[path = "side_support_cell_solve_partition.rs"]
mod partition;
use geometry::{author_exact_pcurve, closing_curve, pcurve_agrees_with_curve, point_on};
use partition::{
    assign_labels, collect_carriers, face_has_live_twin, face_surface_misfits_boundary,
    first_carrier_run, worst_boundary_misfit,
};

/// Outcome summary for diagnostics.
#[derive(Debug, Default)]
pub(crate) struct SupportCellSolveReport {
    /// Cells split off and surfaced.
    pub(crate) cells_created: usize,
    /// Bridge faces surfaced whole (single-carrier boundary).
    pub(crate) faces_surfaced: usize,
    /// Per-face decline reasons (face kept as scaffold).
    pub(crate) declined: Vec<String>,
    /// Refit candidates whose committed surface mis-fit their own
    /// boundary and were re-solved into pocket cells. These were never
    /// honest cap faces — the caller must drop them from the measured
    /// cap set (their remaining seams are pocket-internal structure).
    pub(crate) demoted_refit_candidates: Vec<FaceHandle>,
}

#[derive(Clone)]
struct Carrier {
    surface: SurfaceGeometry,
    source: FaceHandle,
}

#[derive(Clone)]
struct HeSample {
    he: HalfEdgeHandle,
    origin: Point3,
    dest: Point3,
    mid: Option<Point3>,
    carrier_mask: u64,
}

/// Late (acceptance-preparation) pocket completion for the direct-apex
/// dispatch: cell-solve every surviving surface-less non-orphan face,
/// plus any `refit_candidates` (the dispatch's created cap faces) whose
/// committed surface measurably fails to contain its own boundary —
/// those are pocket-annulus residue whose rim lies exactly on the local
/// carriers (measured 0..1e-17 on cube_extreme while the annulus Nurbs
/// was 0.02..0.18 off). Stripping + re-solving them is per-face
/// transactional like everything else in this module.
pub(crate) fn resolve_pocket_faces(
    brep: &mut BRep,
    refit_candidates: &[FaceHandle],
) -> SupportCellSolveReport {
    let gate = fit_gate(brep);
    let mut targets: Vec<FaceHandle> = Vec::new();
    for (handle, face) in brep.faces.iter() {
        if face.surface.is_none() && face_has_live_twin(brep, handle) {
            targets.push(handle);
        }
    }
    for &candidate in refit_candidates {
        if targets.contains(&candidate) {
            continue;
        }
        let misfits = face_surface_misfits_boundary(brep, candidate, gate);
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace] pocket refit candidate {candidate:?}: misfits={misfits} \
                 worst={:?}",
                worst_boundary_misfit(brep, candidate)
            );
        }
        if misfits {
            targets.push(candidate);
        }
    }
    let mut report = SupportCellSolveReport::default();
    let carriers = collect_carriers(brep, &targets);
    if carriers.is_empty() {
        return report;
    }
    for &face in &targets {
        let snapshot = brep.clone();
        let was_refit = brep
            .faces
            .get(face)
            .is_some_and(|data| data.surface.is_some());
        if let Some(data) = brep.faces.get_mut(face) {
            data.surface = None;
        } else {
            continue;
        }
        match solve_one_face(brep, face, &carriers, gate) {
            Ok((cells, surfaced)) => {
                report.cells_created += cells;
                report.faces_surfaced += usize::from(surfaced);
                if was_refit {
                    report.demoted_refit_candidates.push(face);
                }
            }
            Err(reason) => {
                *brep = snapshot;
                report.declined.push(format!("{face:?}: {reason}"));
            }
        }
    }
    // Stale-analytic pcurve repair: Steiner splits and doubled slit
    // records confuse downstream repair heuristics into re-authoring
    // full-extent analytic pcurves on the halves (measured: ShellGap
    // 1e-2 on a split slit edge whose Linear pcurves spanned the whole
    // pre-split segment). For half-edges on analytic (Plane/Cylinder)
    // faces, disagreement between the pcurve and the edge curve has an
    // exact repair — re-project.
    let handles: Vec<HalfEdgeHandle> = brep.half_edges.keys().collect();
    let mut repaired = 0usize;
    for he in handles {
        let analytic = brep
            .half_edges
            .get(he)
            .and_then(|d| brep.loops.get(d.loop_))
            .and_then(|l| brep.faces.get(l.face))
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|s| {
                matches!(
                    s.geometry,
                    SurfaceGeometry::Plane { .. } | SurfaceGeometry::Cylinder { .. }
                )
            });
        if !analytic {
            continue;
        }
        if matches!(pcurve_agrees_with_curve(brep, he, gate), Some(false))
            && author_exact_pcurve(brep, he, gate)
        {
            repaired += 1;
        }
    }
    if repaired > 0 && crate::direct_apex_trace::enabled() {
        eprintln!("[direct-apex-trace] stale-analytic pcurve repair: {repaired} half-edge(s)");
    }
    report
}

fn fit_gate(brep: &BRep) -> f64 {
    brep.tolerance.linear().max(1.0e-9) * 10.0
}

fn he_fits(sample: &HeSample, carrier: &Carrier, gate: f64) -> bool {
    point_on(&carrier.surface, &sample.origin, gate)
        && point_on(&carrier.surface, &sample.dest, gate)
        && sample
            .mid
            .as_ref()
            .is_none_or(|m| point_on(&carrier.surface, m, gate))
}

fn loop_samples(brep: &BRep, face: FaceHandle) -> Result<Vec<HeSample>, String> {
    let loop_ = brep
        .faces
        .get(face)
        .and_then(|f| brep.loops.get(f.outer_loop))
        .ok_or_else(|| "face has no outer loop".to_string())?;
    let start = loop_.half_edge;
    let mut out = Vec::new();
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        let Some(data) = brep.half_edges.get(he) else {
            return Err("loop walk hit a stale half-edge".to_string());
        };
        let origin = brep.vertices[data.origin].position;
        let dest = brep
            .half_edges
            .get(data.next)
            .map(|n| brep.vertices[n.origin].position)
            .ok_or_else(|| "loop walk hit a stale next".to_string())?;
        let mid = brep
            .edges
            .get(data.edge)
            .and_then(|e| e.curve.as_ref())
            .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
        out.push(HeSample {
            he,
            origin,
            dest,
            mid,
            carrier_mask: 0,
        });
        he = data.next;
        if he == start {
            return Ok(out);
        }
    }
    Err("loop did not close".to_string())
}

/// Solve one bridge face. Returns (cells split off, whole-face surfaced).
fn solve_one_face(
    brep: &mut BRep,
    face: FaceHandle,
    carriers: &[Carrier],
    gate: f64,
) -> Result<(usize, bool), String> {
    let mut cells = 0usize;
    // Closing half-edges this solve authored into the parent loop. A
    // run made ONLY of these must not be split again — its cell already
    // exists on the other side, and re-splitting would mint zero-area
    // bigons forever.
    let mut authored: std::collections::BTreeSet<HalfEdgeHandle> =
        std::collections::BTreeSet::new();
    let budget = loop_samples(brep, face)?.len() + 2;
    let mut iterations = 0usize;
    while iterations < budget {
        iterations += 1;
        let mut samples = loop_samples(brep, face)?;
        for sample in &mut samples {
            for (ci, carrier) in carriers.iter().enumerate() {
                if he_fits(sample, carrier, gate) {
                    sample.carrier_mask |= 1 << ci;
                }
            }
        }
        let labels = assign_labels(&samples, carriers.len());
        let bare = carriers.len();
        let n = samples.len();
        if labels.iter().all(|&l| l == labels[0]) {
            if labels[0] == bare {
                if cells == 0 {
                    return Err("no boundary run fits any carrier".to_string());
                }
                return Ok((cells, false));
            }
            // Whole (remaining) face on one carrier.
            assign_cell_surface(brep, face, &carriers[labels[0]], gate)?;
            refresh_pcurves(brep, face)?;
            return Ok((cells, true));
        }
        // Find a carrier run to split off: the first maximal run whose
        // label is a carrier. Runs are contiguous circular index ranges.
        let Some((run_start, run_len, run_label)) =
            first_carrier_run(&labels, bare, &samples, &authored)
        else {
            if cells == 0 {
                return Err("no boundary run fits any carrier".to_string());
            }
            return Ok((cells, false));
        };
        if run_len == n {
            unreachable!("uniform label handled above");
        }
        let closing = split_off_run(
            brep,
            &samples,
            run_start,
            run_len,
            &carriers[run_label],
            gate,
        )?;
        authored.insert(closing);
        cells += 1;
        // Loop again on the shrunken parent face.
    }
    // Budget exhausted. When the terminal remainder is pure wrap (no
    // carrier fits any half-edge), the process CONVERGED — the split
    // count was simply larger than the budget guessed (measured on
    // cube_extreme: 13 legitimate alternating-run splits) and the bare
    // wrap is the tolerated outcome. Only a remainder that still has
    // carrier-fitting runs is a genuine anomaly.
    let samples = loop_samples(brep, face)?;
    let all_bare = samples.iter().all(|sample| {
        carriers
            .iter()
            .all(|carrier| !he_fits(sample, carrier, gate))
    });
    if all_bare && cells > 0 {
        return Ok((cells, false));
    }
    let profile = samples
        .iter()
        .map(|sample| {
            let mut mask = 0u64;
            for (index, carrier) in carriers.iter().enumerate() {
                if he_fits(sample, carrier, gate) {
                    mask |= 1 << index;
                }
            }
            format!(
                "{}{mask}",
                if authored.contains(&sample.he) {
                    "A"
                } else {
                    "r"
                },
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    Err(format!(
        "support-cell split budget exhausted with carrier runs remaining; terminal loop \
         [authored?carrier_mask]: {profile}"
    ))
}

fn split_off_run(
    brep: &mut BRep,
    samples: &[HeSample],
    run_start: usize,
    run_len: usize,
    carrier: &Carrier,
    gate: f64,
) -> Result<HalfEdgeHandle, String> {
    let n = samples.len();
    let first_he = samples[run_start].he;
    let after_idx = (run_start + run_len) % n;
    let mut after_he = samples[after_idx].he;
    let run_origin = samples[run_start].origin;
    let mut run_dest = samples[(run_start + run_len - 1) % n].dest;
    let run_origin_vertex = brep
        .half_edges
        .get(first_he)
        .map(|d| d.origin)
        .ok_or_else(|| "run start half-edge is stale".to_string())?;
    let run_dest_vertex = brep
        .half_edges
        .get(after_he)
        .map(|d| d.origin)
        .ok_or_else(|| "run end half-edge is stale".to_string())?;
    if run_origin_vertex == run_dest_vertex {
        // Pinch: the closing needs two distinct vertices. Shrink the run
        // by its WHOLE last half-edge (no Steiner split — split_edge on
        // the pocket's doubled slit records left half-done record
        // surgery: a halved curve whose edge still referenced the
        // pre-split far vertex, ShellGap 1e-2). The closing then joins
        // origin(first) to origin(last), which are distinct unless the
        // run itself is a single self-loop edge.
        if run_len < 2 {
            return Err("single-edge pinch run cannot be shrunk".to_string());
        }
        let last = &samples[(run_start + run_len - 1) % n];
        let last_origin_vertex = brep
            .half_edges
            .get(last.he)
            .map(|d| d.origin)
            .ok_or_else(|| "pinch run's last half-edge is stale".to_string())?;
        if last_origin_vertex == run_origin_vertex {
            return Err("pinch shrink is still degenerate".to_string());
        }
        after_he = last.he;
        run_dest = last.origin;
    }
    // mef's new edge's primary half-edge runs origin(he1) -> origin(he2)
    // = run_origin -> run_dest; author the curve in that direction.
    let closing = closing_curve(
        &carrier.surface,
        run_origin,
        run_dest,
        gate,
        brep.tolerance.linear(),
    )?;
    let mef = brep
        .mef(first_he, after_he)
        .map_err(|e| format!("mef split failed: {e}"))?;
    brep.edges[mef.edge].curve = Some(closing);
    // mef: the new face carries [first_he .. after_he) + the closing twin;
    // the parent keeps `mef.half_edge` (the closing's primary side).
    assign_cell_surface(brep, mef.face, carrier, gate)?;
    refresh_pcurves(brep, mef.face)?;
    Ok(mef.half_edge)
}

/// Assign the carrier surface to the cell and verify the whole cell
/// boundary lies on it.
fn assign_cell_surface(
    brep: &mut BRep,
    face: FaceHandle,
    carrier: &Carrier,
    gate: f64,
) -> Result<(), String> {
    let samples = loop_samples(brep, face)?;
    for sample in &samples {
        if !he_fits(sample, carrier, gate) {
            return Err(format!(
                "cell postcondition: boundary sample off carrier from {:?} at ({:.4},{:.4},{:.4})",
                carrier.source, sample.origin.x, sample.origin.y, sample.origin.z
            ));
        }
    }
    let orientation = brep
        .faces
        .get(carrier.source)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.orientation)
        .unwrap_or(true);
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: carrier.surface.clone(),
        orientation,
    });
    Ok(())
}

fn refresh_pcurves(brep: &mut BRep, face: FaceHandle) -> Result<(), String> {
    let mut refresh: std::collections::BTreeSet<FaceHandle> = std::collections::BTreeSet::new();
    refresh.insert(face);
    let samples = loop_samples(brep, face)?;
    for sample in samples {
        if let Some(twin_face) = brep
            .half_edges
            .get(sample.he)
            .and_then(|d| brep.half_edges.get(d.twin))
            .and_then(|t| brep.loops.get(t.loop_))
            .map(|l| l.face)
        {
            refresh.insert(twin_face);
        }
    }
    for face in refresh {
        if brep
            .faces
            .get(face)
            .is_none_or(|data| data.surface.is_none())
        {
            continue;
        }
        brep.recompute_face_pcurves(face)
            .map_err(|err| format!("pcurve authoring on {face:?} failed: {err}"))?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "side_support_cell_solve_tests.rs"]
mod tests;
