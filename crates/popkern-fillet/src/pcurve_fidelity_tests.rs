//! Validation of the pcurve-fidelity instrument.
//!
//! An instrument that has only been run on the thing it is
//! investigating is not evidence. These pin a known-good reading and a
//! known-bad one, plus the orientation trap that produced a retracted
//! root cause during the PnS campaign.

use super::*;
use popkern_core::tolerance::Tolerance;

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("unit cube fixture")
}

fn all_faces(brep: &BRep) -> Vec<FaceHandle> {
    brep.faces.iter().map(|(f, _)| f).collect()
}

/// KNOWN GOOD: a unit cube's pcurves are exact, so the instrument must
/// read essentially zero. If this ever drifts, the instrument is lying
/// about everything else it reports.
#[test]
fn a_clean_body_reads_near_zero() {
    let brep = cube();
    let fidelity =
        measure_pcurve_fidelity(&brep, &all_faces(&brep)).expect("a cube has measurable pcurves");
    assert!(
        fidelity.measured >= 12,
        "a cube should measure at least a dozen half-edges, got {}",
        fidelity.measured
    );
    assert!(
        fidelity.worst < 1e-9,
        "a clean body must read near zero, got {:.3e} on {:?}",
        fidelity.worst,
        fidelity.worst_edge
    );
}

/// KNOWN BAD: displace one pcurve and the instrument must find it, and
/// find it on the face that was damaged. Without this the "near zero"
/// result above could equally mean the instrument measures nothing.
#[test]
fn a_damaged_pcurve_is_found_and_located() {
    let mut brep = cube();
    let (face, half_edge) = brep
        .faces
        .iter()
        .find_map(|(f, d)| {
            let loop_data = brep.loops.get(d.outer_loop)?;
            Some((f, loop_data.half_edge))
        })
        .expect("cube has a face with a loop");
    let displaced = {
        use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
        let pc = brep.half_edges[half_edge]
            .pcurve
            .as_ref()
            .expect("cube half-edge has a pcurve")
            .clone();
        // Replace it with a straight segment between its own endpoints
        // slid 0.25 in u. Variant-agnostic on purpose: the point is to
        // put the pcurve somewhere it demonstrably is not supposed to
        // be, whatever form it started in.
        let a = pc.point_at(pc.t_start);
        let b = pc.point_at(pc.t_end);
        Pcurve {
            geometry: PcurveGeometry::Linear {
                start: (a.x + 0.25, a.y),
                end: (b.x + 0.25, b.y),
            },
            t_start: pc.t_start,
            t_end: pc.t_end,
        }
    };
    brep.half_edges[half_edge].pcurve = Some(displaced);

    let fidelity =
        measure_pcurve_fidelity(&brep, &all_faces(&brep)).expect("damaged body still measures");
    assert!(
        fidelity.worst > 1e-3,
        "a displaced pcurve must be detected, got {:.3e}",
        fidelity.worst
    );
    assert_eq!(
        fidelity.worst_face, face,
        "the instrument must point at the face that was damaged"
    );
}

/// THE ORIENTATION TRAP. Reversing a pcurve's control points leaves it
/// tracing the same curve backwards — geometrically identical, and it
/// must still read near zero. A probe that samples both sides forward
/// reads this as ~0.9 * chord and reports a defect that is not there;
/// that error produced a retracted root cause (PK-202).
#[test]
fn a_reversed_pcurve_is_not_a_defect() {
    let mut brep = cube();
    let half_edge = brep
        .faces
        .iter()
        .find_map(|(_, d)| Some(brep.loops.get(d.outer_loop)?.half_edge))
        .expect("cube has a face with a loop");
    let reversed = {
        use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
        let pc = brep.half_edges[half_edge]
            .pcurve
            .as_ref()
            .expect("pcurve present")
            .clone();
        // Rebuild the same UV segment traversed the other way —
        // geometrically identical, parametrically reversed. Built from
        // the original's own endpoints so it is a true reversal
        // whatever variant the original was.
        let a = pc.point_at(pc.t_start);
        let b = pc.point_at(pc.t_end);
        Pcurve {
            geometry: PcurveGeometry::Linear {
                start: (b.x, b.y),
                end: (a.x, a.y),
            },
            t_start: pc.t_start,
            t_end: pc.t_end,
        }
    };
    brep.half_edges[half_edge].pcurve = Some(reversed);

    let fidelity =
        measure_pcurve_fidelity(&brep, &all_faces(&brep)).expect("reversed body still measures");
    assert!(
        fidelity.worst < 1e-9,
        "a REVERSED pcurve traces the same curve and must not read as a defect; \
         got {:.3e}. This is the exact error that produced a retracted root cause.",
        fidelity.worst
    );
}

/// No faces means no measurement — never a zero, which would read as
/// a pass.
#[test]
fn nothing_measurable_returns_none() {
    let brep = cube();
    assert!(measure_pcurve_fidelity(&brep, &[]).is_none());
}
