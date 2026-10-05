//! Post-commit strip-conflict TRIGGER for setback vertex blends
//! (PK-225; demoted from verdict to trigger by PK-518).
//!
//! PK-518 re-measured the class on today's substrate: the rho = 0.35
//! cube corner below now serves self-intersection clean
//! (`validate_all` and OpenCASCADE `BRepAlgoAPI_Check`), and every G2
//! strip — flatter than circular by construction (rho 0.4351 vs
//! 0.4443 on the r = 0.16 triangular-pyramid apex) — serves clean
//! too. The analytic class below therefore no longer declines on its
//! own: when it fires, the surgery measures the committed corner
//! region for self-intersection and declines only on a real one
//! (`setback_vertex_topology`, strip-conflict postcondition). The
//! history below records why the class was first a verdict.
//!
//! The N=3 setback extends the fillet strips to the setback arcs, and
//! FLAT conic sections at small radii genuinely interpenetrate near
//! the shared contact corner (measured r=0.1 rho=0.35: strip-strip
//! SelfIntersection on the served body; rho >= 0.4142 clean).
//!
//! Mechanism record (all measured this campaign, none survived):
//! - full-model / full-domain pair validation: minutes per NURBS
//!   strip pair;
//! - 2-D corner-windowed surface intersector: catches the dirty class
//!   in 0.45 s but proving a CLEAN corner costs ~70 s per
//!   near-tangent pair;
//! - contact-extent sampling: the crossing curve (~1e-2) slips
//!   between affordable grid samples;
//! - signed-straddle against full or kept-filtered neighbour
//!   surfaces: feet beyond the neighbour's material (or over-filtered
//!   kept bands) read arbitrary or no sides.
//!
//! The class boundary is ANALYTIC, and the measurements sit exactly
//! on it: a conic section with rho below the circular equivalent
//! `rho_c = cos(theta/2) / (1 + cos(theta/2))` sags INSIDE the
//! rolling-ball envelope, and the setback corner contract (tangent-
//! coherent contact corners on shared supports) holds on or outside
//! that envelope only. Measured: rho = 0.4142 (= rho_c for the 90°
//! cube corner) serves clean at r ∈ {0.1, 0.15}; rho = 0.35 self-
//! intersects at r ∈ {0.1, 0.12}. The gate measures each strip's
//! actual mid-edge section fullness against its own opening angle —
//! no profile metadata needed, circular strips sit exactly on the
//! boundary and pass.

use popkern_core::geometry::gregory::{Kp20PaperCapReport, Kp20PaperGuideBorder};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

use crate::fillet_surface::FilletSurfaceResult;

/// Paper-cap acceptance before any setback topology commits: the
/// assembly's C0, G1 and border residuals must close at 1e-6 of the
/// cap scale, floored at ten times the INPUT corner mismatch (skew
/// adjacent trims — octagon8_mixed: no exact trim intersection,
/// corners are closest-approach midpoints — bound how far any
/// derivation can close; measured g2 2.1e-7 against a 4.0e-8
/// pure-scale gate), and regularity must stay positive.
///
/// The paper's G2 rows are a construction residual, not the
/// commercial G2 instrument (PK-518): at the regular valence n = 4 on
/// real borders they stay open (g2 5.5e-4 against this gate) while
/// the committed fan's LIVE seams measure inside the 5e-2 commercial
/// G2 gate. G2 is therefore decided downstream by the fan-seam
/// postcondition, alongside the fold and regularity gates.
pub(in crate::construct) fn paper_cap_acceptance(
    border: &Kp20PaperGuideBorder,
    report: &Kp20PaperCapReport,
    star: Point3,
) -> Result<(), String> {
    let scale = border
        .sectors
        .iter()
        .flat_map(|sector| sector.iter())
        .flat_map(|piece| piece.rows.iter().flatten())
        .map(|cp| cp.distance_to(&star))
        .fold(0.0_f64, f64::max);
    let gate = (1.0e-6 * scale.max(f64::MIN_POSITIVE)).max(10.0 * report.input_corner_mismatch);
    if report.max_g1_residual <= gate
        && report.max_c0_residual <= gate
        && report.max_border_residual <= gate
        && report.min_regularity > 0.0
    {
        return Ok(());
    }
    Err(format!(
        "cap gates open: g1={:.3e} g2={:.3e} c0={:.3e} border={:.3e} min_reg={:.3e} gate={gate:.3e}",
        report.max_g1_residual,
        report.max_g2_residual,
        report.max_c0_residual,
        report.max_border_residual,
        report.min_regularity,
    ))
}

/// The fitted cap's self-intersection postcondition (PK-574, scoped
/// by PK-598): the cap's pieces against the whole candidate region
/// (cap plus edge-adjacent neighbours).
///
/// Pairs among the untouched neighbours are re-proved only when the
/// strips are in the flat-section interpenetration class — the one
/// class where two carved strips genuinely cross near the corner (the
/// rho 0.35 cube), exactly the trigger the KP20 fan's strip-conflict
/// postcondition uses. Otherwise the neighbour pairs are the strips
/// and supports the carve already committed, which the fit does not
/// change, and re-proving them cost 0.08-0.64 s per attempt (measured
/// 2026-09-29 on pentagon/hexagon/octagon apexes).
pub(in crate::construct) fn fit_cap_region_self_intersection(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(), String> {
    use crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::{
        candidate_region, validate_candidate_region_self_intersection,
    };
    if flat_section_class(fillet_surfaces).is_some() {
        return validate_candidate_region_self_intersection(brep, cap_faces);
    }
    let region = candidate_region(brep, cap_faces);
    brep.validate_self_intersection_touching(&region, cap_faces)
        .map_err(|error| format!("{error}"))
}

/// Decline reason when any strip's section is flatter than its
/// circular equivalent (the interpenetration class), else `None`.
pub(in crate::construct) fn flat_section_class(
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<String> {
    const FLATNESS_SLACK: f64 = 1.0e-3;
    for (k, fillet) in fillet_surfaces.iter().enumerate() {
        let Some((rho, rho_circular)) = worst_section_fullness(&fillet.surface) else {
            continue;
        };
        if rho < rho_circular - FLATNESS_SLACK {
            return Some(format!(
                "strip {k} section fullness rho={rho:.4} is below its circular equivalent \
                 {rho_circular:.4} — flat sections sag inside the rolling-ball envelope and \
                 the setback kept-domains interpenetrate at the corner"
            ));
        }
    }
    None
}

/// `(rho, rho_circular)` for the strip's least-full sampled section.
/// The section axis is detected per surface (sweeps carry the
/// cross-section on either parameter): the axis with the larger
/// turning angle is the section.
fn worst_section_fullness(surface: &popkern_core::geometry::SurfaceGeometry) -> Option<(f64, f64)> {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        return None;
    }
    let section = |along: f64, across: f64, v_is_section: bool| -> Point3 {
        if v_is_section {
            surface.point_at(u0 + (u1 - u0) * along, v0 + (v1 - v0) * across)
        } else {
            surface.point_at(u0 + (u1 - u0) * across, v0 + (v1 - v0) * along)
        }
    };
    let turning = |v_is_section: bool| -> f64 {
        let p0 = section(0.5, 0.0, v_is_section);
        let pm = section(0.5, 0.5, v_is_section);
        let p1 = section(0.5, 1.0, v_is_section);
        let (Some(a), Some(b)) = ((pm - p0).normalized(), (p1 - pm).normalized()) else {
            return 0.0;
        };
        a.dot(&b).clamp(-1.0, 1.0).pacos()
    };
    let v_is_section = turning(true) >= turning(false);
    let mut worst: Option<(f64, f64)> = None;
    for station in [0.25, 0.5, 0.75] {
        let Some(pair) = section_fullness(|t| section(station, t, v_is_section)) else {
            continue;
        };
        worst = Some(match worst {
            None => pair,
            Some(best) => {
                if pair.0 - pair.1 < best.0 - best.1 {
                    pair
                } else {
                    best
                }
            }
        });
    }
    worst
}

/// Fullness of one cross-section: `rho` = shoulder offset over apex
/// offset (the classic conic fullness), `rho_circular` from the
/// section's own opening angle.
fn section_fullness(section: impl Fn(f64) -> Point3) -> Option<(f64, f64)> {
    let p0 = section(0.0);
    let p1 = section(1.0);
    let shoulder = section(0.5);
    let h = 1.0e-4;
    let t0: Vec3 = ((section(h) - p0) * (1.0 / h)).normalized()?;
    let t1: Vec3 = ((section(1.0 - h) - p1) * (-1.0 / h)).normalized()?;
    // Apex: intersection of the two end tangent lines (closest
    // approach midpoint for robustness).
    let w = p1 - p0;
    let (a, b, c) = (t0.dot(&t0), t0.dot(&t1), t1.dot(&t1));
    let (d, e) = (t0.dot(&w), t1.dot(&w));
    let denom = a * c - b * b;
    if denom.abs() < 1e-14 {
        return None;
    }
    let s = (c * d - b * e) / denom;
    let t = (b * d - a * e) / denom;
    let apex_a = p0 + t0 * s;
    let apex_b = p1 + t1 * t;
    let apex = Point3::new(
        0.5 * (apex_a.x + apex_b.x),
        0.5 * (apex_a.y + apex_b.y),
        0.5 * (apex_b.z + apex_a.z),
    );
    let chord_mid = Point3::new(
        0.5 * (p0.x + p1.x),
        0.5 * (p0.y + p1.y),
        0.5 * (p0.z + p1.z),
    );
    let apex_offset = apex - chord_mid;
    let apex_len = apex_offset.length();
    if apex_len < 1e-12 {
        return None;
    }
    let rho = (shoulder - chord_mid).dot(&(apex_offset * (1.0 / apex_len))) / apex_len;
    // Opening angle from the end tangents: the section turns by
    // pi - theta between its ends.
    let cos_turn = t0.dot(&(t1 * -1.0)).clamp(-1.0, 1.0);
    let half = 0.5 * cos_turn.pacos();
    let c_half = half.pcos();
    let rho_circular = c_half / (1.0 + c_half);
    Some((rho, rho_circular))
}

#[cfg(test)]
#[path = "setback_commit_gate_tests.rs"]
mod tests;
