//! The variational-fit cap on a carved setback region (PK-574).
//!
//! The carved blend face is bounded by the n setback sides. The core fit
//! (`popkern_core::geometry::variational_cap`) returns n NURBS sectors
//! meeting at a centre o; this module subdivides the blend face into
//! those sectors (spokes from o to each side's midpoint m_s), gives every
//! cap edge its exact iso-curve and every cap half-edge its exact Linear
//! iso pcurve, re-authors the neighbours' pcurves, and certifies the
//! cap's seams (see `setback_fit_certify`). Any failure is a cap-stage
//! decline; the surgery restores the snapshot. The corner-patched mode
//! (PK-619) builds several such pieces, see `setback_fit_patched`.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::variational_cap::{
    CapEnergy, VariationalCap, VariationalCapOptions, VariationalCapTargets, fit_variational_cap,
};

use super::*;
use crate::blend_continuity_evidence::ContinuityGates;
use crate::construct::kp20_paper_border_extract::Kp20SetbackArc;
use crate::construct::setback_commit_gate::fit_cap_region_self_intersection;
use crate::construct::setback_decline::CapCommitError;
use crate::construct::setback_fit_cap::fit_sides;
use crate::types::{SetbackFitCertified, SetbackFitMode};

#[path = "setback_fit_certify.rs"]
mod setback_fit_certify;
#[path = "setback_fit_patched.rs"]
mod setback_fit_patched;
#[path = "setback_fit_pcurves.rs"]
mod setback_fit_pcurves;
#[path = "setback_fit_sectors.rs"]
mod setback_fit_sectors;
use setback_fit_sectors::build_sector_faces;
#[path = "setback_fit_tilt.rs"]
mod setback_fit_tilt;
use setback_fit_pcurves::author_cap_pcurves;
use setback_fit_tilt::{tilt_rim_row, tilted_face};

/// A certified fit serve.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FitServe {
    pub mode: SetbackFitMode,
    pub certified: SetbackFitCertified,
    /// Per cap face (in cap-face order): which piece it belongs to and
    /// its sector there (PK-619).
    pub pieces: Vec<FacePiece>,
}

/// A cap face's place in its cap: the piece (0 for a star cap or the
/// main cap of a corner-patched one, k for its k-th corner patch), its
/// sector within that piece, and the piece's sector count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FacePiece {
    pub piece: usize,
    pub sector: usize,
    pub sectors: usize,
}

/// The star cap's faces: one piece.
fn star_pieces(n: usize) -> Vec<FacePiece> {
    (0..n)
        .map(|sector| FacePiece {
            piece: 0,
            sector,
            sectors: n,
        })
        .collect()
}

/// Fit, build and certify the cap over `blend_face`.
pub(super) fn commit_fit_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    arcs: &[Kp20SetbackArc],
    mode: SetbackFitMode,
    (r_min, rival): (f64, Option<f64>),
    tolerance: &Tolerance,
) -> Result<(Vec<FaceHandle>, FitServe), CapCommitError> {
    let gates = ContinuityGates::commercial();
    // Neighbour orientation from the carved topology: the face across
    // each side's middle from the blend face.
    let orientation: Vec<bool> = arcs
        .iter()
        .map(|arc| side_neighbour_orientation(brep, blend_face, arc))
        .collect::<Result<_, String>>()?;
    let outward = |k: usize, _: Point3, _: popkern_core::geometry::Vec3| orientation[k];
    let sides = fit_sides(arcs, fillet_surfaces, &outward)?;
    let (targets, options) = fit_settings(mode, (r_min, rival), tolerance);
    if mode == SetbackFitMode::TangentCornerPatched {
        return setback_fit_patched::commit_patched_fan(
            brep,
            blend_face,
            fillet_surfaces,
            &sides,
            (&targets, &options),
            tolerance,
        );
    }
    let mut cap = fit_variational_cap(&sides, &targets, &options, tolerance.linear())
        .map_err(|e| format!("variational fit: {e}"))?;
    if crate::dispatch_overrides::fit_cap_rim_tilted() {
        tilt_rim_row(&mut cap.sectors[0]);
        tilted_face(&mut cap.faces[0], &cap.sectors[0]);
    }
    if mode == SetbackFitMode::CurvatureWeighted {
        curvature_witness_gate(&cap, gates.g2_curvature)?;
    }
    let cap_faces = build_sector_faces(brep, blend_face, &cap, &sides)?;
    author_cap_pcurves(brep, &cap_faces, &cap)?;
    super::author_setback_pcurves(brep, &[], &cap_faces)?;
    author_neighbour_pcurves(brep, &cap_faces, tolerance)?;
    super::setback_boundary_tolerance::record_boundary_tolerances(brep, &cap_faces);
    let certified = setback_fit_certify::certify_cap(brep, &cap_faces, &gates, &cap.report)?;
    let required = mode.level();
    if level_rank(certified.level) < level_rank(required) {
        return Err(CapCommitError {
            reason: format!(
                "variational fit ({mode:?}) certified {:?}, needs {required:?}: rim (g0 {:.3e}, g1 {:.3e}, g2 {:.3e}), spoke (g0 {:.3e}, g1 {:.3e}, g2 {:.3e})",
                certified.level,
                certified.rim[0],
                certified.rim[1],
                certified.rim[2],
                certified.spoke[0],
                certified.spoke[1],
                certified.spoke[2],
            ),
            fit_certified: Some(certified),
        });
    }
    if let Err(intersection) = fit_cap_region_self_intersection(brep, &cap_faces, fillet_surfaces) {
        return Err(CapCommitError {
            reason: format!("fit cap self-intersection postcondition failed: {intersection}"),
            fit_certified: Some(certified),
        });
    }
    let pieces = star_pieces(cap_faces.len());
    Ok((
        cap_faces,
        FitServe {
            mode,
            certified,
            pieces,
        },
    ))
}

/// The fit's targets and options for `mode`, with the fairness and
/// thin-sector limits for a corner of smallest radius `r_min` and the
/// rival fan's reading.
fn fit_settings(
    mode: SetbackFitMode,
    (r_min, rival): (f64, Option<f64>),
    tolerance: &Tolerance,
) -> (VariationalCapTargets, VariationalCapOptions) {
    let (targets, options) = match mode {
        SetbackFitMode::CurvatureWeighted => (
            VariationalCapTargets::g2(
                tolerance.linear(),
                1e-4,
                ContinuityGates::commercial().g2_curvature,
            ),
            VariationalCapOptions {
                max_refinement_levels: 0,
                energy: CapEnergy::MeanCurvatureVariation,
                // PK-614: one freeze bent the blunt pentagon to 3.03 / r_min
                // (three: 1.31) and the hexagon h1.2 to 12.4 (4.39).
                reference_passes: 3,
                // PK-595: each rim owns its corner 2-jet (dual-jet cells, degree 9).
                corner_cells: true,
                degree: 9,
                ..Default::default()
            },
        ),
        SetbackFitMode::Tangent => (
            VariationalCapTargets::g1(tolerance.linear(), 1e-6),
            VariationalCapOptions {
                max_refinement_levels: 2,
                energy: CapEnergy::Bending,
                reference_passes: 3, // PK-605: one freeze left a rim curvature layer
                ..Default::default()
            },
        ),
        // PK-615: the plain tangent fit folds or creases on mixed-radius
        // G1 corners (pyramid-3 h0.4 wide: rim g1 0.99); with corner
        // cells it certifies G1 at 3.2 / r_min.
        // PK-619: the corner-patched cap's pieces are corner-cell tangent
        // fits.
        SetbackFitMode::TangentCornerCells | SetbackFitMode::TangentCornerPatched => (
            VariationalCapTargets::g1(tolerance.linear(), 1e-6),
            VariationalCapOptions {
                max_refinement_levels: 0,
                energy: CapEnergy::Bending,
                reference_passes: 3,
                corner_cells: true,
                degree: 9,
                ..Default::default()
            },
        ),
    };
    // PK-618: a layout bending past the fairness bound is refit with the
    // spokes' feet nearest the fitted centre; the fairer cap is kept. A
    // sector under the served regularity floor serves only under the
    // bound or the unfair fan the rung holds (a short strip between long
    // bridges gives its sector a thin edge, not a fold; PK-621: the
    // octagon h1.2 variable's thin fit reads 15.1 against the fan's 17.9).
    let bound = crate::construct::setback_cap_fairness::CAP_FAIRNESS_BOUND;
    let usable = r_min.is_finite() && r_min > 0.0;
    let fairness_limit = usable.then(|| bound / r_min);
    let thin_limit = usable.then(|| rival.map_or(bound, |fan| fan.max(bound)) / r_min);
    let options = VariationalCapOptions {
        layout_refit_curvature: fairness_limit,
        thin_sector_curvature: thin_limit,
        ..options
    };
    (targets, options)
}

/// Co-parameterise every neighbour half-edge along the rim with the
/// rim edge's new curve (the fit's iso-curve). The rim geometry was
/// replaced, so a neighbour's split or projected pcurve is stale or
/// missing — analytic neighbours (cylinder, torus strips) are not
/// re-authored by the carve's pcurve pass at all. The certifier reads
/// the seam under fraction correspondence, which needs exactly this
/// co-parameterisation. The verification gate is ten linear
/// tolerances: the rim curve lies off the neighbour by the fit's G0
/// residual, which the certifier then bounds.
fn author_neighbour_pcurves(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    tolerance: &Tolerance,
) -> Result<(), String> {
    let gate = 10.0 * tolerance.linear();
    for &face in cap_faces {
        for he in face_loop_half_edges(brep, face)? {
            let twin = brep.half_edges[he].twin;
            let neighbour = brep.loops[brep.half_edges[twin].loop_].face;
            if cap_faces.contains(&neighbour) {
                continue;
            }
            let edge = brep.half_edges[he].edge;
            let curve = brep.edges[edge]
                .curve
                .clone()
                .ok_or("rim edge has no curve")?;
            let authored = popkern_core::topology::edge_pcurve_reauthor::coparameterised_pcurves(
                brep, edge, &curve, gate, tolerance,
            )
            .ok_or_else(|| {
                format!("rim edge {edge:?}: neighbour pcurve does not verify within {gate:.1e}")
            })?;
            for (h, pcurve) in authored {
                if h == twin {
                    brep.half_edges[twin].pcurve = Some(pcurve);
                }
            }
        }
    }
    Ok(())
}

/// Decline a curvature-weighted fit before any topology when its own
/// rim samples already prove G2 cannot certify. A sample is an attained
/// value, so it bounds the seam's supremum from below; the fit reports
/// the spectral norm of the second-fundamental-form jump, which is at
/// most twice the largest component the certifier bounds (2×2
/// symmetric), hence the factor 2. Measured (E4): certification was
/// 70 % of the rung's time, mostly spent on these doomed caps.
fn curvature_witness_gate(cap: &VariationalCap, g2_gate: f64) -> Result<(), CapCommitError> {
    match cap.report.rim.iter().find(|r| r.g2 > 2.0 * g2_gate) {
        Some(r) => Err(CapCommitError::from(format!(
            "curvature-weighted fit: sampled rim curvature jump {:.3e} on side {} ({:?} half, t = {:.3}) exceeds twice the G2 gate {g2_gate:.1e}; G2 cannot certify",
            r.g2, r.side, r.half, r.parameter
        ))),
        None => Ok(()),
    }
}

fn level_rank(level: crate::types::ContinuityLevel) -> u8 {
    match level {
        crate::types::ContinuityLevel::G0 => 0,
        crate::types::ContinuityLevel::G1 => 1,
        crate::types::ContinuityLevel::G2 => 2,
    }
}

/// `SurfaceOnFace::orientation` of the face across the blend face from
/// the middle of `arc` (true: the neighbour's natural normal is its
/// outward normal).
fn side_neighbour_orientation(
    brep: &BRep,
    blend_face: FaceHandle,
    arc: &Kp20SetbackArc,
) -> Result<bool, String> {
    let curve = match &arc.bridge {
        Some(bridge) => bridge.curve.clone(),
        None => arc.section_curve(),
    };
    let (a, b) = curve.parameter_range();
    let mid = curve.point_at(0.5 * (a + b));
    let mut best: Option<(f64, HalfEdgeHandle)> = None;
    for he in face_loop_half_edges(brep, blend_face)? {
        let edge = brep.half_edges[he].edge;
        let Some(c) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let d = curve_distance(&c.geometry, c.t_start, c.t_end, mid);
        if best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, he));
        }
    }
    let (_, he) = best.ok_or("blend face has no curved boundary edge")?;
    let twin = brep.half_edges[he].twin;
    let face = brep.loops[brep.half_edges[twin].loop_].face;
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.orientation)
        .ok_or_else(|| "side neighbour face has no surface".to_string())
}

/// Parameter of the point of `curve` (on [0, 1]) nearest `p`: a dense
/// scan refined by ternary search, snapped to the ends within 1e-12.
fn param_on(curve: &popkern_core::geometry::NurbsCurve, p: Point3) -> f64 {
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..=256 {
        let t = i as f64 / 256.0;
        let d = curve.point_at(t).distance_to(&p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let (mut lo, mut hi) = (
        (best.1 - 1.0 / 256.0).max(0.0),
        (best.1 + 1.0 / 256.0).min(1.0),
    );
    for _ in 0..80 {
        let t1 = lo + (hi - lo) / 3.0;
        let t2 = hi - (hi - lo) / 3.0;
        if curve.point_at(t1).distance_to(&p) < curve.point_at(t2).distance_to(&p) {
            hi = t2;
        } else {
            lo = t1;
        }
    }
    let t = 0.5 * (lo + hi);
    if t < 1e-12 {
        0.0
    } else if t > 1.0 - 1e-12 {
        1.0
    } else {
        t
    }
}
