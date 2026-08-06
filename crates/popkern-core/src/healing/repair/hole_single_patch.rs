//! Single-patch G1 cap for 4-sided holes (PK-224).
//!
//! A 4-sided hole does not need a sector complex: one rectangular
//! degree-5 Gregory patch over the whole ring has NO internal seams
//! (the fan-crease class cannot exist) and meets each surviving
//! neighbour tangentially through its boundary ribbon rows. The
//! Gregory inner pairs exist precisely to let two adjacent sides'
//! ribbons disagree at a shared corner.
//!
//! # Corner reality
//!
//! At a hole corner the two neighbour walls meet at an angle; no cap
//! can be G1 to both at the corner point. Each rim side's quintic is
//! therefore built with ZERO end legs (`[S, S, A, B, E, E]` — the
//! curve geometry is unchanged, the parameterization has a cusp at
//! the corners), which kills the adjacent side's corner-leg
//! contribution to this side's cross field: the cross-derivative
//! along each open side is then ribbon + tangential terms only —
//! G1 to the neighbour by construction away from the corner points.
//! Design: `docs/architecture/hole-cap-single-patch-g1.md`.

use crate::error::CoreError;
use crate::geometry::gregory::{GregoryPatch, RectangularGregoryDeg5};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::topology::BRep;
use crate::topology::coherence::CoherenceValidationScope;

use super::hole_fill::DetectedHole;
use super::n_sided_fill_g1::{G1HoleFillError, boundary_edges_with_frames};

/// Mid-band fractions used for the measured cap↔neighbour G1 gate.
/// Corner bands are excluded by the structural corner argument (see
/// module doc); the open middle of each side is where a G1 cap must
/// actually be tangent to its wall.
const G1_SAMPLE_FRACTIONS: [f64; 5] = [0.25, 0.375, 0.5, 0.625, 0.75];

/// Measured cap↔neighbour G1 bar, `|sin|` of the normal angle.
const G1_GATE: f64 = 5.0e-2;

/// Why the single-patch cap declined. Stringly summarized into the
/// heal report like the other rungs.
#[derive(Debug, thiserror::Error)]
pub enum SinglePatchCapError {
    /// Only 4-sided holes are served by this rung.
    #[error("single-patch cap serves 4-sided holes; this hole has {sides}")]
    UnsupportedSideCount {
        /// Sides of the detected hole.
        sides: usize,
    },
    /// The boundary frames could not be built.
    #[error("single-patch cap: {0}")]
    Frames(#[from] G1HoleFillError),
    /// A side's quintic fit missed the sampled rim curve.
    #[error(
        "single-patch cap: side {side} quintic fit residual {residual:.3e} exceeds {budget:.3e}"
    )]
    FitResidual {
        /// Side index in ring order.
        side: usize,
        /// Worst sample distance to the fitted Bezier.
        residual: f64,
        /// Linear tolerance budget scaled for fitting.
        budget: f64,
    },
    /// A side's ribbon direction could not be oriented.
    #[error("single-patch cap: side {side} has no usable ribbon direction")]
    RibbonDegenerate {
        /// Side index in ring order.
        side: usize,
    },
    /// The measured cap↔neighbour G1 exceeded the gate mid-side.
    #[error(
        "single-patch cap: side {side} measured cap-neighbour G1 {worst:.3e} exceeds {gate:.1e}"
    )]
    G1Exceeded {
        /// Side index in ring order.
        side: usize,
        /// Worst |sin| over the mid-band samples.
        worst: f64,
        /// The gate.
        gate: f64,
    },
    /// The topology or pcurve layer declined.
    #[error("single-patch cap: {0}")]
    Topology(#[from] CoreError),
    /// The committed cap failed the coherence gate.
    #[error(
        "single-patch cap: coherence gate failed — {missing_curves} curve(s), {missing_pcurves} pcurve(s) missing, worst L2 {l2_worst:.3e}"
    )]
    CoherenceRejected {
        /// L1 violations on the cap face.
        missing_curves: usize,
        /// L3 violations on the cap face.
        missing_pcurves: usize,
        /// Worst sampled boundary disagreement.
        l2_worst: f64,
    },
    /// Post-commit validation failed.
    #[error("single-patch cap: post-commit validation failed: {0}")]
    PostValidationFailed(#[from] crate::error::ValidationError),
}

/// A committed single-patch cap.
#[derive(Debug, Clone)]
pub struct SinglePatchCapCommit {
    /// The new cap face.
    pub face: crate::topology::FaceHandle,
    /// Worst L2 disagreement from the coherence gate.
    pub worst_l2_gap: f64,
    /// Worst measured cap↔neighbour G1 (`|sin|`) over the mid-band
    /// samples of all four sides.
    pub worst_mid_band_g1: f64,
}

/// Build and transactionally commit the single-patch G1 cap for a
/// 4-sided hole. Any decline leaves `brep` untouched.
pub fn commit_single_patch_hole_cap(
    brep: &mut BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Result<SinglePatchCapCommit, SinglePatchCapError> {
    if hole.boundary_edges.len() != 4 {
        return Err(SinglePatchCapError::UnsupportedSideCount {
            sides: hole.boundary_edges.len(),
        });
    }
    let frames = boundary_edges_with_frames(brep, hole, tolerance)?;

    // Quintic side rows with zero end legs, chained CCW.
    let mut boundaries = [[Point3::new(0.0, 0.0, 0.0); 6]; 4];
    let mut chords = [0.0_f64; 4];
    for (side, frame) in frames.iter().enumerate() {
        boundaries[side] = fit_quintic_zero_leg(frame, side, tolerance)?;
        chords[side] = frame.start.distance_to(&frame.end);
    }

    // Per-side ribbon directions at the 1/3 and 2/3 stations,
    // oriented away from the neighbour face across the rim.
    let mut ribbon = [[Vec3::new(0.0, 0.0, 0.0); 2]; 4];
    for (side, frame) in frames.iter().enumerate() {
        let witness = neighbour_witness(brep, hole, side);
        for (station, f) in [(0usize, 1.0 / 3.0), (1usize, 2.0 / 3.0)] {
            let n =
                interpolate_frame_vector(frame.normal_start, frame.normal_mid, frame.normal_end, f);
            let t = interpolate_frame_vector(
                frame.tangent_start,
                frame.tangent_mid,
                frame.tangent_end,
                f,
            );
            let mut w = n.cross(&t);
            let w_len = w.length();
            if w_len < 1e-12 {
                return Err(SinglePatchCapError::RibbonDegenerate { side });
            }
            w = w * (1.0 / w_len);
            if let Some(wit) = witness
                && w.dot(&wit) < 0.0
            {
                w = w * -1.0;
            }
            ribbon[side][station] = w;
        }
    }

    let patch = assemble_patch(&boundaries, &chords, &ribbon, tolerance)
        .map_err(|e| SinglePatchCapError::Topology(CoreError::PreconditionFailed(e.to_string())))?;

    // Measured G1 gate BEFORE committing topology: sample each side's
    // mid band and compare the cap normal against the neighbour's.
    let cap_geometry =
        SurfaceGeometry::Gregory(GregoryPatch::RectangularDeg5(Box::new(patch.clone())));
    let mut worst_mid_band_g1 = 0.0_f64;
    for side in 0..4 {
        let mut worst = 0.0_f64;
        for &f in &G1_SAMPLE_FRACTIONS {
            let (u, v) = side_uv(side, f);
            let cap_n = cap_geometry.normal_at(u, v);
            let Some(cap_n) = cap_n.normalized() else {
                continue;
            };
            let rim_point = cap_geometry.point_at(u, v);
            let Some(neighbour) = neighbour_surface_geometry(brep, hole, side) else {
                continue;
            };
            let Some((nu, nv)) = neighbour.closest_point(&rim_point) else {
                continue;
            };
            let Some(wall_n) = neighbour.normal_at(nu, nv).normalized() else {
                continue;
            };
            let sin = cap_n.cross(&wall_n).length();
            worst = worst.max(sin);
        }
        worst_mid_band_g1 = worst_mid_band_g1.max(worst);
        if worst > G1_GATE {
            return Err(SinglePatchCapError::G1Exceeded {
                side,
                worst,
                gate: G1_GATE,
            });
        }
    }

    // Transactional commit: single-face ring cap + surface + pcurves
    // + coherence + validation on a candidate clone.
    let mut candidate = brep.clone();
    let ring = candidate.cap_boundary_ring(&hole.boundary_edges)?;
    if let Some(face) = candidate.faces.get_mut(ring.face) {
        face.surface = Some(SurfaceOnFace {
            geometry: cap_geometry.clone(),
            orientation: true,
        });
    }
    candidate.recompute_face_pcurves(ring.face)?;
    // The orientation flag is derived state (PK-244); reconcile it
    // from the committed winding rather than trusting the stamp.
    let _ = candidate.reconcile_orientation_flags();

    let scope = CoherenceValidationScope::Faces([ring.face].into_iter().collect());
    let report = candidate.validate_coherence(&scope);
    let worst_l2_gap = report
        .l2_face_boundary_disagreements
        .iter()
        .map(|d| d.worst_gap)
        .fold(0.0_f64, f64::max);
    if !report.is_empty() {
        return Err(SinglePatchCapError::CoherenceRejected {
            missing_curves: report.l1_missing_edge_curves.len(),
            missing_pcurves: report.l3_missing_pcurves.len(),
            l2_worst: worst_l2_gap,
        });
    }
    candidate.validate()?;

    *brep = candidate;
    Ok(SinglePatchCapCommit {
        face: ring.face,
        worst_l2_gap,
        worst_mid_band_g1,
    })
}

/// Quintic Bezier `[S, S, A, B, E, E]` for the side. Zero end legs
/// by construction (see module doc); `A`, `B` interpolate the rim at
/// the 1/3 and 2/3 stations, and a GEOMETRIC residual gate (sampled
/// point-to-curve distance — the cusped parameterization cannot and
/// need not match the rim's own parameter) verifies containment.
/// `orient_loop` flips a side's endpoints without flipping its dense
/// samples, so the station lookup re-orients them by endpoint
/// distance first.
fn fit_quintic_zero_leg(
    frame: &crate::geometry::FillBoundaryEdge,
    side: usize,
    tolerance: f64,
) -> Result<[Point3; 6], SinglePatchCapError> {
    let bern = |i: usize, t: f64| -> f64 {
        let s = 1.0 - t;
        let c = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0][i];
        c * t.powi(i as i32) * s.powi(5 - i as i32)
    };

    // Dense rim samples in THIS side's orientation.
    let mut samples: Vec<(f64, Point3)> = frame.extra_samples.clone();
    let reversed = samples
        .first()
        .map(|&(_, p)| p.distance_to(&frame.end) < p.distance_to(&frame.start))
        .unwrap_or(false);
    if reversed {
        for (f, _) in &mut samples {
            *f = 1.0 - *f;
        }
        samples.reverse();
    }
    samples.push((0.5, frame.mid));

    // Station targets at 1/3 and 2/3 (nearest dense samples).
    let station_point = |f0: f64| -> Point3 {
        samples
            .iter()
            .min_by(|a, b| (a.0 - f0).abs().total_cmp(&(b.0 - f0).abs()))
            .map(|&(_, p)| p)
            .unwrap_or(frame.mid)
    };
    let (t1, t2) = (1.0 / 3.0, 2.0 / 3.0);
    let p1 = station_point(t1);
    let p2 = station_point(t2);

    // Solve A, B so the curve interpolates p1 at t1 and p2 at t2.
    let a11 = bern(2, t1);
    let a12 = bern(3, t1);
    let a21 = bern(2, t2);
    let a22 = bern(3, t2);
    let det = a11 * a22 - a12 * a21;
    if det.abs() < 1e-18 {
        return Err(SinglePatchCapError::FitResidual {
            side,
            residual: f64::INFINITY,
            budget: tolerance,
        });
    }
    let fixed = |t: f64, sc: f64, ec: f64| -> f64 {
        (bern(0, t) + bern(1, t)) * sc + (bern(4, t) + bern(5, t)) * ec
    };
    let solve_coord = |ps: f64, pe: f64, v1: f64, v2: f64| -> (f64, f64) {
        let r1 = v1 - fixed(t1, ps, pe);
        let r2 = v2 - fixed(t2, ps, pe);
        ((a22 * r1 - a12 * r2) / det, (a11 * r2 - a21 * r1) / det)
    };
    let (ax, bx) = solve_coord(frame.start.x, frame.end.x, p1.x, p2.x);
    let (ay, by) = solve_coord(frame.start.y, frame.end.y, p1.y, p2.y);
    let (az, bz) = solve_coord(frame.start.z, frame.end.z, p1.z, p2.z);
    let cps = [
        frame.start,
        frame.start,
        Point3::new(ax, ay, az),
        Point3::new(bx, by, bz),
        frame.end,
        frame.end,
    ];

    // Geometric residual: every dense sample must lie on the curve
    // (min distance over a dense evaluation — parameterization-free).
    let eval = |t: f64| -> Point3 {
        let mut p = Point3::new(0.0, 0.0, 0.0);
        for (i, cp) in cps.iter().enumerate() {
            let w = bern(i, t);
            p = Point3::new(p.x + w * cp.x, p.y + w * cp.y, p.z + w * cp.z);
        }
        p
    };
    let budget = (tolerance * 10.0).max(1e-9);
    let mut residual = 0.0_f64;
    for &(_, p) in &samples {
        // Coarse scan + windowed refinement (parameterization-free
        // point-to-curve distance; a coarse polyline minimum alone
        // reads the sampling spacing, not the geometry).
        let mut t_min = 0.0_f64;
        let mut d_min = f64::INFINITY;
        for k in 0..=64 {
            let t = k as f64 / 64.0;
            let d = eval(t).distance_to(&p);
            if d < d_min {
                d_min = d;
                t_min = t;
            }
        }
        let mut window = 1.0 / 64.0;
        for _ in 0..12 {
            let lo = (t_min - window).max(0.0);
            let hi = (t_min + window).min(1.0);
            for k in 0..=8 {
                let t = lo + (hi - lo) * k as f64 / 8.0;
                let d = eval(t).distance_to(&p);
                if d < d_min {
                    d_min = d;
                    t_min = t;
                }
            }
            window *= 0.4;
        }
        residual = residual.max(d_min);
    }
    if residual > budget {
        return Err(SinglePatchCapError::FitResidual {
            side,
            residual,
            budget,
        });
    }
    Ok(cps)
}

/// Assemble the degree-5 Gregory patch: chained boundaries + ribbon
/// inner pairs per corner.
fn assemble_patch(
    boundaries: &[[Point3; 6]; 4],
    chords: &[f64; 4],
    ribbon: &[[Vec3; 2]; 4],
    tolerance: f64,
) -> Result<RectangularGregoryDeg5, crate::geometry::gregory::GregoryBuildError> {
    // Ribbon member for side s at its start (station 0) and end
    // (station 1): a point over the side's own row shifted by the
    // ribbon direction. The tangential part is irrelevant to G1 (it
    // lies in the neighbour plane); the ribbon part carries it.
    let member = |s: usize, station: usize| -> Point3 {
        let q = &boundaries[s];
        let base = if station == 0 {
            Point3::new(
                0.5 * (q[1].x + q[2].x),
                0.5 * (q[1].y + q[2].y),
                0.5 * (q[1].z + q[2].z),
            )
        } else {
            Point3::new(
                0.5 * (q[3].x + q[4].x),
                0.5 * (q[3].y + q[4].y),
                0.5 * (q[3].z + q[4].z),
            )
        };
        let h = chords[s] / 3.0;
        base + ribbon[s][station] * h
    };
    // Corner incident sides (CCW chaining b0..b3, corner k = start
    // of side k): side k starts at corner k, side (k+3)%4 ends at
    // corner k. Member-to-pair-slot assignment is pinned by the
    // calibration test `deg5_pair_member_convention` — on the u=0
    // edge the blend limits to the P_u member, so the side lying
    // along u=0 (boundaries[0]) owns the P_u slots of corners 0/1,
    // and so on around the ring.
    let inner_cps: [(Point3, Point3); 4] = [
        (member(0, 0), member(3, 1)),
        (member(0, 1), member(1, 0)),
        (member(2, 1), member(1, 1)),
        (member(2, 0), member(3, 0)),
    ];
    RectangularGregoryDeg5::from_boundary_beziers(*boundaries, inner_cps, tolerance.max(1e-9))
}

/// `(u, v)` on the patch for fraction `f` along ring side `s`,
/// following the `from_boundary_beziers` convention (side 0 = u=0
/// edge walking v, side 1 = v=1 edge walking u, side 2 = u=1 edge
/// walking v backwards, side 3 = v=0 edge walking u backwards).
fn side_uv(side: usize, f: f64) -> (f64, f64) {
    match side {
        0 => (0.0, f),
        1 => (f, 1.0),
        2 => (1.0, 1.0 - f),
        _ => (1.0 - f, 0.0),
    }
}

/// The neighbour face's surface beside ring side `s`.
fn neighbour_surface_geometry<'a>(
    brep: &'a BRep,
    hole: &DetectedHole,
    side: usize,
) -> Option<&'a SurfaceGeometry> {
    let edge = brep.edges.get(*hole.boundary_edges.get(side)?)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let face = brep.loops.get(he.loop_)?.face;
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// Witness direction pointing across the rim away from the neighbour
/// face: from the neighbour's outer-loop centroid toward the rim
/// side's midpoint.
fn neighbour_witness(brep: &BRep, hole: &DetectedHole, side: usize) -> Option<Vec3> {
    let edge = brep.edges.get(*hole.boundary_edges.get(side)?)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let face = brep.loops.get(he.loop_)?.face;
    let positions = brep.face_outer_loop_positions(face);
    if positions.is_empty() {
        return None;
    }
    let n = positions.len() as f64;
    let centroid = Point3::new(
        positions.iter().map(|p| p.x).sum::<f64>() / n,
        positions.iter().map(|p| p.y).sum::<f64>() / n,
        positions.iter().map(|p| p.z).sum::<f64>() / n,
    );
    let curve = edge.curve.as_ref()?;
    let mid = {
        use crate::geometry::curve::CurveEval;
        (&curve.geometry as &dyn CurveEval).point_at(0.5 * (curve.t_start + curve.t_end))
    };
    (mid - centroid).normalized()
}

/// Linear blend of the three sampled frame vectors at fraction `f`.
fn interpolate_frame_vector(start: Vec3, mid: Vec3, end: Vec3, f: f64) -> Vec3 {
    let v = if f <= 0.5 {
        let t = f * 2.0;
        start * (1.0 - t) + mid * t
    } else {
        let t = (f - 0.5) * 2.0;
        mid * (1.0 - t) + end * t
    };
    v.normalized().map_or(mid, |unit| unit)
}

#[cfg(test)]
#[path = "hole_single_patch_tests.rs"]
mod tests;
