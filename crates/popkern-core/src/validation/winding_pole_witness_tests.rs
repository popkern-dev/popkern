//! Pole-side witness pins (PK-882).

use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn line(u0: f64, v0: f64, u1: f64, v1: f64) -> crate::geometry::pcurve::Pcurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
        weights: vec![1.0, 1.0],
    }
    .into()
}

/// The cone as `BRep::cone` built it before PK-881: the cap canonical on
/// the base circle, the lateral pcurves laid for it, the cap chart
/// (x, y) — both base uses ran the rim backwards.
fn pre_pk881_cone() -> BRep {
    let mut cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    let tau = std::f64::consts::TAU;
    let h = 2.0;
    let base = cone
        .edges
        .keys()
        .find(|&e| {
            let he = cone.edges[e].half_edge;
            cone.half_edges[he].origin == cone.half_edges[cone.half_edges[he].twin].origin
        })
        .unwrap();
    let lateral = cone.edges[base].half_edge;
    let cap = cone.half_edges[lateral].twin;
    cone.edges[base].half_edge = cap;
    cone.half_edges[lateral].pcurve = Some(line(0.0, h, tau, h));
    let up = cone.half_edges[lateral].next;
    let down = cone.half_edges[up].next;
    cone.half_edges[up].pcurve = Some(line(tau, h, tau, 0.0));
    cone.half_edges[down].pcurve = Some(line(0.0, 0.0, 0.0, h));
    let cap_face = cone.loops[cone.half_edges[cap].loop_].face;
    if let Some(sof) = cone.faces[cap_face].surface.as_mut()
        && let SurfaceGeometry::Plane { v_axis, .. } = &mut sof.geometry
    {
        *v_axis = crate::geometry::Vec3::new(0.0, 1.0, 0.0);
    }
    cone
}

fn lateral_verdict(cone: &BRep) -> (Option<bool>, &'static str) {
    let audit = cone.audit_winding_honesty();
    let lateral = audit
        .iter()
        .find(|v| {
            matches!(
                cone.faces[v.face].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cone { .. })
            )
        })
        .unwrap();
    (lateral.honest, lateral.reason)
}

#[test]
fn the_pole_side_witness_reads_the_cone_lateral() {
    let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
    assert_eq!(lateral_verdict(&cone), (Some(true), "pole-side-witness"));
    assert!(cone.validate_winding_honesty_all().is_empty());
}

#[test]
fn the_pole_side_witness_flags_the_pre_pk881_cone() {
    let cone = pre_pk881_cone();
    assert_eq!(lateral_verdict(&cone), (Some(false), "pole-side-witness"));
    assert!(!cone.validate_winding_honesty_all().is_empty());
}
