//! Construct-time spoke/strip footprint clearance (PK-261).
//!
//! At steep radius ratios an inner spoke's interior can dive through a
//! neighbouring fillet strip's live footprint while both of its ends
//! stay legal (measured pyramid3_extreme 1:10: the shared sector spoke
//! penetrates the middle r=0.10 strip to -4.71e-2 at t=0.787; the pole
//! itself is radially inside the strip's cylinder but axially OUTSIDE
//! its footprint, so a pole-side detour exists). The serve-time
//! footprint repair (PK-251) moves interior control points only and
//! provably stalls on such boundary-row crossings, so the clearance
//! must happen where the spoke is authored.
//!
//! Mechanism: after corner-compatible shaping, measure every spoke
//! against every corner strip's live footprint — the strip's analytic
//! implicit form, restricted to an arc-derived angular window and a
//! vertex-bounded axial window (a trim-blind signed distance against
//! the infinite primitive is meaningless; see `StripFootprint` for
//! why the angular window must come from the end arc rather than the
//! loop vertices). The legal side of each strip is the side its OWN
//! RIM occupies there: the prepared outer boundary curves are the
//! committed boundary contract — trusted legal — and the flanking
//! bridge boundaries measurably traverse the mid strip's footprint
//! window on the clear side (pyramid3_extreme: +5.4e-2 / +1.06e-1
//! against the offending strip). The majority implicit-form sign of
//! in-footprint rim samples beyond the band fixes the side once per
//! corner. (The serve-time repair's authority, the component material
//! sense, is unusable mid-candidate — the component walk panics on
//! transient construction handles — and a spoke's own approach sign
//! is unusable too: a spoke anchored on the strip's end arc enters
//! the footprint already inside material, which a self-referential
//! side rule would bless as legal. Both were measured, not guessed.)
//! Violating spokes are corrected by pushing their FREE
//! control points: the anchor handle only along its locked in-plane
//! direction (seam-end G1, PK-125) and the pole-side handle freely;
//! anchor and pole positions never move, degree and knots never
//! change (`spoke_in_basis` hostability), so spoke G0 reconciliation,
//! the pole meeting point, and cross-spoke pairing are unaffected.
//! Twin records of a shared spoke stay bit-identical because the
//! correction is a deterministic function of the (identical) curve
//! and the corner-global strip set.
//!
//! Correction-only containment boundary: a spoke this pass cannot
//! clear records an informational residual diagnostic and the
//! candidate proceeds unchanged — downstream the serving
//! self-intersection retry detects any surviving crossing and serves
//! the legacy dispatch, which is the pinned honest contract for the
//! class. Kill-switch: `POPKERN_SPOKE_CLEARANCE=0`.

use popkern_core::geometry::analytic_surfaces::analytic_signed_distance;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::corner::{GeneratedInnerSpoke, PreparedExactBoundaryInnerSpokes};
use super::sector::ExactBoundaryCapSectorInput;
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

use footprint::collect_strip_footprints;

#[path = "spoke_strip_clearance_footprint.rs"]
mod footprint;

/// Interior samples per spoke per measurement round (endpoints are
/// locked and excluded — the anchor legitimately sits ON the strip's
/// end arc).
const SPOKE_SAMPLES: usize = 64;

/// Push rounds per spoke before recording a residual.
const MAX_ROUNDS_PER_SPOKE: usize = 24;

/// Grazing band as a multiple of the model's linear tolerance:
/// penetration within the band is the validator's sub-tolerance
/// grazing class (PK-230, material-two-sided rule), not a crossing.
const GRAZING_BAND_FACTOR: f64 = 4.0;

/// A corner strip's live footprint: analytic surface plus the UV
/// window it actually occupies. The angular (u) window comes from the
/// strip's own END ARC — the prepared fillet boundary curve, exact
/// on-surface construction data — NOT from projecting the face's
/// outer-loop vertices: mid-candidate those loops are the Phase C
/// merge heuristic's incoherent state (measured pyramid3_extreme:
/// loop vertices 0.1063 off-surface, and the strip's far-end trim
/// curve sweeps the complementary side of the cylinder entirely —
/// the arc is the only trustworthy angular authority). The axial (v)
/// window keeps the vertex bound — it spans both strip ends and is
/// robust to the same vertex noise.
struct StripFootprint {
    face: FaceHandle,
    surface: SurfaceGeometry,
    u_periodic: bool,
    u_range: (f64, f64),
    v_range: (f64, f64),
}

/// Worst material-side violation of a spoke against the strip set.
#[derive(Clone)]
struct SpokeViolation {
    depth: f64,
    t: f64,
    out_dir: Vec3,
    strip_index: usize,
}

/// Measure every prepared inner spoke against the corner's strip
/// footprints and reshape the free control points of violating spokes
/// until they clear the grazing band. Returns informational
/// diagnostics only — this pass never declines the candidate.
pub(super) fn apply_spoke_strip_clearance(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    inner: &mut PreparedExactBoundaryInnerSpokes,
) -> Vec<ApexConstructionDiagnostic> {
    if std::env::var("POPKERN_SPOKE_CLEARANCE").as_deref() == Ok("0") {
        return Vec::new();
    }
    let band = brep.tolerance.linear() * GRAZING_BAND_FACTOR;
    let strips = collect_strip_footprints(brep, sectors, &inner.sectors);
    let trace = std::env::var("POPKERN_DIRECT_APEX_TRACE").as_deref() == Ok("1");
    if strips.is_empty() {
        return Vec::new();
    }
    let sides = rim_legal_sides(&inner.sectors, &strips, band);
    if trace {
        eprintln!(
            "[spoke-clearance] strips={:?} rim sides={sides:?} band={band:.3e}",
            strips.iter().map(|s| s.face).collect::<Vec<_>>()
        );
    }
    let mut diagnostics = Vec::new();
    for sector in inner.sectors.iter_mut() {
        let sector_name = sector.sector_name.clone();
        for (side, spoke) in [
            ("FilletEnd", &mut sector.fillet_inner_spoke),
            ("FaceEdge", &mut sector.face_inner_spoke),
        ] {
            let diagnostic = clear_one_spoke(spoke, &strips, &sides, band, &sector_name, side);
            if trace && let Some(diagnostic) = &diagnostic {
                eprintln!("[spoke-clearance] {sector_name}/{side}: {diagnostic:?}");
            }
            if let Some(diagnostic) = diagnostic {
                diagnostics.push(diagnostic);
            }
        }
    }
    diagnostics
}

/// Clear one spoke against the strip set. Returns an informational
/// diagnostic when the spoke was corrected or a residual remains;
/// `None` when the spoke was already clear.
fn clear_one_spoke(
    spoke: &mut GeneratedInnerSpoke,
    strips: &[StripFootprint],
    sides: &[Option<f64>],
    band: f64,
    sector_name: &str,
    side: &str,
) -> Option<ApexConstructionDiagnostic> {
    let initial = worst_violation(&spoke.curve, strips, sides, band)?;
    let mut curve = cubic_bezier_form(&spoke.curve)?;
    let mut worst = worst_violation(&curve, strips, sides, band).or(Some(initial.clone()));
    let mut previous_depth = f64::INFINITY;
    let mut rounds = 0;
    while let Some(violation) = worst {
        if rounds >= MAX_ROUNDS_PER_SPOKE || violation.depth >= previous_depth - 1.0e-9 {
            // Stalled or out of budget: record the residual and leave
            // the curve at its best measured state.
            spoke.curve = curve;
            return Some(clearance_diagnostic(
                "spoke_footprint_residual",
                sector_name,
                side,
                strips[violation.strip_index].face,
                violation.depth,
                violation.t,
                band,
            ));
        }
        previous_depth = violation.depth;
        if !push_free_control_points(&mut curve, &violation, band) {
            spoke.curve = curve;
            return Some(clearance_diagnostic(
                "spoke_footprint_endpoint_locked",
                sector_name,
                side,
                strips[violation.strip_index].face,
                violation.depth,
                violation.t,
                band,
            ));
        }
        rounds += 1;
        worst = worst_violation(&curve, strips, sides, band);
    }
    spoke.curve = curve;
    Some(clearance_diagnostic(
        "spoke_footprint_cleared",
        sector_name,
        side,
        strips[initial.strip_index].face,
        initial.depth,
        initial.t,
        band,
    ))
}

/// The legal side of each strip, read off the corner's own rim: the
/// majority implicit-form sign of prepared outer-boundary samples
/// within the strip's AXIAL window beyond the grazing band. The side
/// vote deliberately ignores the angular window — legal rim curves
/// clear the patch angularly (that is what makes them legal), but the
/// material side of a rolling-ball strip is a property of the
/// primitive along its whole length, so a legal on-body point
/// alongside the strip still fixes which side is clear (measured
/// pyramid3_extreme: the flanking bridge boundary rides the offending
/// strip's axial span at +5.4e-2). The strip's own end arc
/// contributes nothing (its samples sit within the band). `None` when
/// no rim sample resolves the side; such strips are skipped rather
/// than guessed at.
fn rim_legal_sides(
    sectors: &[super::corner::PreparedExactBoundarySectorBoundary],
    strips: &[StripFootprint],
    band: f64,
) -> Vec<Option<f64>> {
    let mut votes = vec![(0usize, 0usize); strips.len()];
    for sector in sectors {
        for boundary in [&sector.fillet_boundary, &sector.face_boundary] {
            for k in 0..=SPOKE_SAMPLES {
                let t = k as f64 / SPOKE_SAMPLES as f64;
                let point = boundary.curve.point_at(t);
                for (index, strip) in strips.iter().enumerate() {
                    if !strip.contains_axial_projection(&point) {
                        continue;
                    }
                    let Some(sd) = analytic_signed_distance(&strip.surface, &point) else {
                        continue;
                    };
                    if sd > band {
                        votes[index].0 += 1;
                    } else if sd < -band {
                        votes[index].1 += 1;
                    }
                }
            }
        }
    }
    votes
        .into_iter()
        .map(|(positive, negative)| match positive.cmp(&negative) {
            std::cmp::Ordering::Greater => Some(1.0),
            std::cmp::Ordering::Less => Some(-1.0),
            std::cmp::Ordering::Equal => None,
        })
        .collect()
}

/// Deepest wrong-side penetration of the spoke's interior samples
/// inside any strip's footprint window, beyond the grazing band.
/// "Wrong side" is relative to the strip's fixed approach side.
fn worst_violation(
    curve: &NurbsCurve,
    strips: &[StripFootprint],
    sides: &[Option<f64>],
    band: f64,
) -> Option<SpokeViolation> {
    let mut worst: Option<SpokeViolation> = None;
    for k in 1..SPOKE_SAMPLES {
        let t = k as f64 / SPOKE_SAMPLES as f64;
        let point = curve.point_at(t);
        for (strip_index, strip) in strips.iter().enumerate() {
            let Some(side) = sides[strip_index] else {
                continue;
            };
            if !strip.contains_footprint_projection(&point) {
                continue;
            }
            let Some((sd, gradient)) = implicit_outward(&strip.surface, &point) else {
                continue;
            };
            let depth = -(sd * side);
            if depth <= band {
                continue;
            }
            if worst.as_ref().is_none_or(|w| depth > w.depth) {
                worst = Some(SpokeViolation {
                    depth,
                    t,
                    out_dir: gradient * side,
                    strip_index,
                });
            }
        }
    }
    worst
}

/// Implicit signed distance and its normalized gradient (finite
/// differences — generic across the analytic types). The gradient
/// points toward the positive implicit side.
fn implicit_outward(surface: &SurfaceGeometry, p: &Point3) -> Option<(f64, Vec3)> {
    let sd = analytic_signed_distance(surface, p)?;
    const H: f64 = 1.0e-6;
    let gradient = Vec3::new(
        (analytic_signed_distance(surface, &Point3::new(p.x + H, p.y, p.z))? - sd) / H,
        (analytic_signed_distance(surface, &Point3::new(p.x, p.y + H, p.z))? - sd) / H,
        (analytic_signed_distance(surface, &Point3::new(p.x, p.y, p.z + H))? - sd) / H,
    );
    let gradient = gradient.normalized()?;
    Some((sd, gradient))
}

/// The spoke as a 4-control-point cubic Bézier: shaped spokes already
/// are one; straight generated spokes are elevated exactly. Any other
/// form is left alone (`None`).
fn cubic_bezier_form(curve: &NurbsCurve) -> Option<NurbsCurve> {
    if curve.degree == 3 && curve.control_points.len() == 4 {
        return Some(curve.clone());
    }
    if curve.degree == 1 && curve.control_points.len() == 2 {
        let a = curve.control_points[0];
        let b = curve.control_points[1];
        let third = (b - a) * (1.0 / 3.0);
        return Some(NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![a, a + third, b - third, b],
            weights: vec![1.0; 4],
        });
    }
    None
}

/// Push the free control points so the curve at the violation
/// parameter moves out of material past the grazing band: the anchor
/// handle only along its locked direction, the pole handle freely,
/// distributed by their Bernstein responses. Returns `false` when the
/// violation sits too close to the locked endpoints to respond.
fn push_free_control_points(curve: &mut NurbsCurve, violation: &SpokeViolation, band: f64) -> bool {
    let t = violation.t;
    let b1 = 3.0 * t * (1.0 - t) * (1.0 - t);
    let b2 = 3.0 * t * t * (1.0 - t);
    let anchor = curve.control_points[0];
    let start_dir = (curve.control_points[1] - anchor)
        .normalized()
        .or_else(|| (curve.control_points[3] - anchor).normalized());
    let Some(start_dir) = start_dir else {
        return false;
    };
    let along = violation.out_dir.dot(&start_dir);
    let w1 = b1 * along.max(0.0);
    let w2 = b2;
    let denominator = w1 * w1 + w2 * w2;
    if denominator < 1.0e-9 {
        return false;
    }
    let needed = violation.depth + band;
    curve.control_points[1] = curve.control_points[1] + start_dir * (needed * w1 / denominator);
    curve.control_points[2] =
        curve.control_points[2] + violation.out_dir * (needed * w2 / denominator);
    true
}

#[cfg(test)]
#[path = "spoke_strip_clearance_tests.rs"]
mod tests;

fn clearance_diagnostic(
    code: &str,
    sector_name: &str,
    side: &str,
    strip: FaceHandle,
    depth: f64,
    t: f64,
    band: f64,
) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
        stage: "spoke_strip_clearance".to_string(),
        code: code.to_string(),
        detail_code: None,
        sector_name: sector_name.to_string(),
        edge: None,
        support_half_edge: None,
        support_face: Some(strip),
        residual: Some(depth),
        outer_side: Some(side.to_string()),
        sample_parameter: Some(t),
        tolerance: Some(band),
        shared_fact_ids: Vec::new(),
    })
}
