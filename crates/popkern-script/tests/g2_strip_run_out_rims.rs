//! PK-396/PK-397 — the served G2 strip run-out rims and support-face
//! pcurves of the rim+verticals box class are exact.
//!
//! Two producers left them wrong: the face-trim stage handed every trim
//! piece the whole rail's support-side pcurve (image spanning the
//! original rail — the PK-391 gaps), and the edge reconcile mis-read an
//! exact run-out rim as 2.88e-2 off the G2 loft (a wrong-minimum
//! projection), fitted a numeric intersection from that seed and wrote
//! it over the exact rim.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn boxed(w: f64, d: f64, h: f64, r: f64) -> Model {
    let mut m = Model::new();
    let plate = m
        .profile(Plane::XY)
        .rect(w, d)
        .sketch("plate")
        .expect("sketch")
        .extrude(&mut m, End::Blind(h.into()))
        .expect("extrude");
    plate
        .fillet_continuity(
            &mut m,
            edges()
                .straight()
                .at_height(h)
                .or(edges().straight().direction(Dir::Z)),
            r,
            ContinuityKind::G2,
            "rounds",
        )
        .expect("G2 fillet serves");
    m
}

#[test]
fn run_out_rims_lie_on_both_surfaces_and_rails_carry_no_gap() {
    for (w, d, h, r) in [
        (2.4, 1.6, 0.4, 0.08),
        (1.0, 1.0, 0.4, 0.08),
        (1.0, 1.0, 1.0, 0.1),
    ] {
        let mut m = boxed(w, d, h, r);
        let b = m.to_brep().expect("brep");
        let mut rims = 0;
        for (eh, e) in &b.edges {
            let he = &b.half_edges[e.half_edge];
            let tw = &b.half_edges[he.twin];
            let a = b.vertices[he.origin].position;
            let z = b.vertices[tw.origin].position;
            if a.z.abs() > 1e-9 || z.z.abs() > 1e-9 {
                continue;
            }
            let faces = [b.loops[he.loop_].face, b.loops[tw.loop_].face];
            let strip = faces.iter().any(|f| {
                matches!(
                    b.faces[*f].surface.as_ref().map(|s| &s.geometry),
                    Some(SurfaceGeometry::Nurbs(_))
                )
            });
            if !strip {
                continue;
            }
            rims += 1;
            let c = e.curve.as_ref().expect("rim curve");
            for f in faces {
                let s = b.faces[f].surface.as_ref().expect("surface");
                for i in 0..=16 {
                    let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 16.0;
                    let p = c.geometry.point_at(t);
                    let dist = s
                        .geometry
                        .closest_point(&p)
                        .map(|(u, v)| s.geometry.point_at(u, v).distance_to(&p))
                        .expect("closest point");
                    assert!(
                        dist < 1e-9,
                        "box {w}x{d}x{h} r={r}: rim {eh:?} is {dist:.3e} off face {f:?}"
                    );
                }
            }
        }
        assert_eq!(
            rims, 8,
            "box {w}x{d}x{h}: four strips × two half-section rims"
        );
        // The r-proportional rail artifact (0.64–0.73·r) is gone: no
        // seam carries a gap within an order of magnitude of it. What
        // remains (≤ 1.1e-4 with sub-gate tangent jumps on the unit
        // cube's corner fans) is the PK-392 corner-class residual.
        let report = m
            .evidence(&m.feature("rounds").expect("label"))
            .expect("evidence");
        let c = report
            .blend
            .as_ref()
            .expect("blend")
            .continuity
            .as_ref()
            .expect("continuity");
        assert!(
            c.worst_g0_gap < 1e-3,
            "box {w}x{d}x{h} r={r}: worst seam gap {:.3e} (rail artifact would be {:.3e})",
            c.worst_g0_gap,
            0.64 * r
        );
    }
}
