//! Serving-boundary pcurve image-honesty audit over the fastener
//! generator bodies (PK-490).
//!
//! A served body's half-edge pcurves are load-bearing state for every
//! chart-trusting consumer (trim loops, chart subdivision in later
//! booleans, STEP export). The mesher survives a wrong-branch pcurve
//! because its ring build works from canonical edge samples — so a
//! pcurve whose image runs hundreds of chart spans from its own edge
//! curve is LATENT: watertight today, torn or misclassified the
//! moment anything follows the pcurve (measured pre-PK-468: a
//! groove-flank use ended at u = −388.8 with image deviation 4.7e7 on
//! a watertight M5 countersunk serve).
//!
//! The audit: for every half-edge use carrying a pcurve, sample the
//! pcurve at 0, ¼, ½, ¾, 1, map through the face surface, and require
//! each image to lie within the audit band of the NEAREST point of
//! its own 3-D edge curve (nearest-point, not fraction-matched, so
//! legal reversed/reparameterised pcurves pass; branch lies still
//! fail by construction — a wrong preimage lands far from the whole
//! curve).

use popkern_core::geometry::{CurveEval, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_part::fastener::{CountersunkScrew, HexBolt, HexNut, SocketCapScrew, ThreadSpec};

/// Worst per-use image deviation over the whole body, with the
/// offending use named for the failure message.
fn worst_pcurve_image_deviation(brep: &popkern_core::topology::BRep) -> (f64, String) {
    let mut worst = 0.0_f64;
    let mut worst_at = String::from("<none>");
    for (fh, face) in brep.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        for lh in loops {
            let Some(ld) = brep.loops.get(lh) else {
                continue;
            };
            let start = ld.half_edge;
            let mut cur = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let Some(hd) = brep.half_edges.get(cur) else {
                    break;
                };
                if let (Some(pc), Some(c)) = (
                    hd.pcurve.as_ref(),
                    brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()),
                ) {
                    for k in 0..=4 {
                        let f = k as f64 / 4.0;
                        let q = pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * f);
                        let img = sof.geometry.point_at(q.x, q.y);
                        let mut nearest = f64::MAX;
                        for j in 0..=64 {
                            let t = c.t_start + (c.t_end - c.t_start) * j as f64 / 64.0;
                            nearest = nearest.min(c.geometry.point_at(t).distance_to(&img));
                        }
                        if nearest > worst {
                            worst = nearest;
                            worst_at = format!(
                                "face {fh:?} he {cur:?} edge {:?} pcurve f={f:.2} uv ({:.4},{:.4})",
                                hd.edge, q.x, q.y
                            );
                        }
                    }
                }
                cur = hd.next;
                if cur == start {
                    break;
                }
            }
        }
    }
    (worst, worst_at)
}

/// The audit band: chord sampling of the nearest-point search plus
/// fit slop. 64-segment sampling of the longest fastener arcs leaves
/// ~1e-3 chord sag; the historical defect class sits 5+ orders above.
const AUDIT_BAND: f64 = 5.0e-2;

fn assert_pcurves_honest(brep: &popkern_core::topology::BRep, label: &str) {
    let (worst, at) = worst_pcurve_image_deviation(brep);
    assert!(
        worst < AUDIT_BAND,
        "{label}: pcurve image deviation {worst:.6} at {at} exceeds the audit band {AUDIT_BAND}"
    );
}

#[test]
fn m5_countersunk_pcurve_uses_are_image_honest() {
    let spec = ThreadSpec::parse("M5").expect("M5");
    let screw = CountersunkScrew::new(spec, 10.0)
        .thread_length(3.0)
        .build(Tolerance::default())
        .expect("countersunk screw");
    assert_pcurves_honest(screw.brep(), "m5 countersunk");
}

#[test]
fn m4_socket_cap_pcurve_uses_are_image_honest() {
    let spec = ThreadSpec::parse("M4").expect("M4");
    let screw = SocketCapScrew::new(spec, 8.0)
        .thread_length(3.0)
        .build(Tolerance::default())
        .expect("socket cap screw");
    assert_pcurves_honest(screw.brep(), "m4 socket cap");
}

#[test]
fn m6_bolt_and_nut_pcurve_uses_are_image_honest() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let bolt = HexBolt::new(spec.clone(), 12.0)
        .thread_length(6.0)
        .build(Tolerance::default())
        .expect("hex bolt");
    assert_pcurves_honest(bolt.brep(), "m6 hex bolt");
    let nut = HexNut::new(spec)
        .clearance(0.15)
        .build(Tolerance::default())
        .expect("nut");
    assert_pcurves_honest(nut.brep(), "m6 nut");
}
