use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A unit-square sheet with authored geometry (plane + line curves),
/// the state the history sketch path hands to extrude.
fn square_sheet() -> (BRep, crate::topology::FaceHandle) {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    brep.attach_geometry().unwrap();
    brep.attach_planar_pcurves().unwrap();
    (brep, mef.face)
}

/// PK-316: extrude authors its own geometry — the result validates
/// WITHOUT any post-pass, every face carries a distinct correct
/// surface, and every edge carries a curve. (Before: every mef face
/// inherited the base plane and new edges had no curves — the history
/// replay path shipped validate-red bodies.)
#[test]
fn extrude_authors_geometry_and_validates() {
    let (mut brep, face) = square_sheet();
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();

    brep.validate()
        .expect("extruded body validates with no post-pass");

    for (eh, ed) in &brep.edges {
        assert!(ed.curve.is_some(), "edge {eh:?} must carry a curve");
    }
    // Side faces: four distinct planes, none of them z-normal.
    assert_eq!(result.side_faces.len(), 4);
    for &f in &result.side_faces {
        let sof = brep.faces[f].surface.as_ref().expect("side surface");
        let SurfaceGeometry::Plane { normal, .. } = &sof.geometry else {
            panic!("straight-edge side face must be planar");
        };
        assert!(
            normal.z.abs() < 1e-9,
            "side face normal {normal:?} must be lateral, not the base plane's"
        );
    }
    // Cap: the base plane translated to z = 1.
    let cap = brep.faces[result.cap_face].surface.as_ref().unwrap();
    let SurfaceGeometry::Plane { origin, normal, .. } = &cap.geometry else {
        panic!("cap must be planar");
    };
    assert!(
        (origin.z - 1.0).abs() < 1e-9,
        "cap origin at z=1, got {origin:?}"
    );
    assert!(normal.z.abs() > 0.99, "cap normal along z, got {normal:?}");
}

/// The full-validation screen is clean too (self-intersection etc.).
#[test]
fn extruded_sheet_passes_full_validation() {
    let (mut brep, face) = square_sheet();
    brep.extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    let findings = brep.validate_all();
    assert!(
        findings.is_empty(),
        "full validation clean, got {findings:?}"
    );
}

/// The unit square with its `y = 0` side bulging outward as a
/// circular arc (centre `(0.5, 0.5, 0)`, radius `√0.5`, through
/// `(0.5, 0.5 − √0.5, 0)`), lying in the sheet plane unless `tilt`
/// rotates the arc plane about the chord (an arc that still meets the
/// same two vertices but whose normal is no longer the sheet normal).
fn square_sheet_with_arc(tilt: f64) -> (BRep, crate::topology::FaceHandle) {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    let (mut brep, face) = square_sheet();
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let (eh, forward) = brep
        .edges
        .iter()
        .find_map(|(eh, ed)| {
            let he = &brep.half_edges[ed.half_edge];
            let p = brep.vertices[he.origin].position;
            let q = brep.vertices[brep.half_edges[he.twin].origin].position;
            if p.distance_to(&a) < 1e-12 && q.distance_to(&b) < 1e-12 {
                Some((eh, true))
            } else if p.distance_to(&b) < 1e-12 && q.distance_to(&a) < 1e-12 {
                Some((eh, false))
            } else {
                None
            }
        })
        .expect("bottom edge");
    // Arc plane: spanned by the chord (x) and a direction rotated by
    // `tilt` from +y towards +z; centre sits `0.5` along that
    // direction from the chord midpoint, radius √0.5.
    let up = Vec3::new(0.0, tilt.cos(), tilt.sin());
    let normal = Vec3::new(1.0, 0.0, 0.0).cross(&up);
    let center = Point3::new(0.5, 0.0, 0.0) + up * 0.5;
    let radius = 0.5_f64.sqrt();
    // Angles are measured in the plane's own frame; probe the frame
    // the kernel uses by evaluating the chord endpoints.
    let probe = CurveGeometry::Arc {
        center,
        normal,
        radius,
        start_angle: 0.0,
        end_angle: std::f64::consts::TAU,
    };
    let angle_of = |p: Point3| {
        let mut best = (f64::MAX, 0.0);
        for i in 0..3600 {
            let t = i as f64 * std::f64::consts::TAU / 3600.0;
            let d = probe.point_at(t).distance_to(&p);
            if d < best.0 {
                best = (d, t);
            }
        }
        best.1
    };
    let (t_a, t_b) = (angle_of(a), angle_of(b));
    // Sweep from a to b through the side AWAY from the square
    // (through `centre − up·r`): pick the direction whose midpoint
    // lands there.
    let outward = center - up * radius;
    let ccw = if t_b > t_a {
        t_b - t_a
    } else {
        t_b + std::f64::consts::TAU - t_a
    };
    let mid_ccw = probe.point_at(t_a + ccw / 2.0);
    let (start_angle, end_angle) = if mid_ccw.distance_to(&outward) < 1e-6 {
        (t_a, t_a + ccw)
    } else {
        (t_a, t_a - (std::f64::consts::TAU - ccw))
    };
    let (start_angle, end_angle) = if forward {
        (start_angle, end_angle)
    } else {
        (end_angle, start_angle)
    };
    brep.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    });
    for fh in brep.faces.keys().collect::<Vec<_>>() {
        for heh in brep.half_edges.keys().collect::<Vec<_>>() {
            brep.half_edges[heh].pcurve = None;
        }
        brep.recompute_face_pcurves(fh).unwrap();
    }
    (brep, face)
}

/// PK-377: an arc base edge extruded along its plane normal authors a
/// TRUE `Cylinder` side face (centre, axis, radius exact), not a ruled
/// NURBS — and the body validates.
#[test]
fn arc_edge_extruded_along_its_normal_promotes_to_cylinder() {
    let (mut brep, face) = square_sheet_with_arc(0.0);
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    brep.validate().expect("arc-sided extrusion validates");
    let mut cylinders = 0;
    for &f in &result.side_faces {
        if let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } = &brep.faces[f].surface.as_ref().unwrap().geometry
        {
            cylinders += 1;
            assert!(origin.distance_to(&Point3::new(0.5, 0.5, 0.0)) < 1e-12);
            assert!((axis.z.abs() - 1.0).abs() < 1e-12);
            assert!((radius - 0.5_f64.sqrt()).abs() < 1e-12);
        }
    }
    assert_eq!(cylinders, 1, "exactly the arc side promotes");
    // The cap's arc edge is the base arc translated, still an Arc.
    let arcs = brep
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(crate::geometry::curve::CurveGeometry::Arc { .. })
            )
        })
        .count();
    assert_eq!(arcs, 2);
    for e in brep.edges.values() {
        let c = e.curve.as_ref().unwrap();
        let _ = c.geometry.point_at(c.t_start);
    }
}

/// PK-377 fallback boundary: an arc whose plane is tilted off the
/// displacement sweeps an ELLIPTIC cylinder; it keeps the exact ruled
/// NURBS rather than a wrong analytic promotion.
#[test]
fn obliquely_extruded_arc_keeps_the_ruled_nurbs() {
    let (mut brep, face) = square_sheet_with_arc(0.3);
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    let nurbs = result
        .side_faces
        .iter()
        .filter(|&&f| {
            matches!(
                brep.faces[f].surface.as_ref().unwrap().geometry,
                SurfaceGeometry::Nurbs(_)
            )
        })
        .count();
    let cylinders = result
        .side_faces
        .iter()
        .filter(|&&f| {
            matches!(
                brep.faces[f].surface.as_ref().unwrap().geometry,
                SurfaceGeometry::Cylinder { .. }
            )
        })
        .count();
    assert_eq!((nurbs, cylinders), (1, 0));
}

/// A bigon lamina — one straight edge and one arc between the SAME two
/// vertices (the line-plus-arc "D"; a split circle is the same class)
/// — extrudes with each top edge carrying ITS OWN base curve. Pairing
/// top edges to base edges by vertex pair handed both the arc: the
/// straight wall's rim lay 0.5 off its plane and the body measured
/// wrong while the geometric validator (a full cylinder contains both
/// arcs) stayed green.
#[test]
fn bigon_lamina_top_edges_carry_their_own_base_curves() {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    let mut brep = BRep::new(Tolerance::default());
    let a = Point3::new(0.0, 0.0, 0.0);
    let b = Point3::new(1.0, 0.0, 0.0);
    let mvfs = brep.mvfs(a).unwrap();
    let mev = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, b).unwrap();
    let mef = brep.mef(mev.half_edge_twin, mev.half_edge).unwrap();
    // Straight edge (mev) and a semicircle (mef) bulging to +y.
    brep.edges[mev.edge].curve = Some(crate::constructors::line_between(a, b));
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let (u, v) = crate::geometry::analytic_curves::circle_basis(normal);
    let center = Point3::new(0.5, 0.0, 0.0);
    let d = b - center;
    let t_b = d.dot(&v).atan2(d.dot(&u));
    // From b counter-clockwise (through +y) to a.
    brep.edges[mef.edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal,
            radius: 0.5,
            start_angle: t_b,
            end_angle: t_b + std::f64::consts::PI,
        },
        t_start: t_b,
        t_end: t_b + std::f64::consts::PI,
    });
    // The face whose loop runs a → b along the line (then back along
    // the arc through +y) winds counter-clockwise about +z.
    let ccw_face = |brep: &BRep, fh: crate::topology::FaceHandle| {
        let he = brep.loops[brep.faces[fh].outer_loop].half_edge;
        let hed = &brep.half_edges[he];
        let along_line = hed.edge == mev.edge;
        let from_a = brep.vertices[hed.origin].position.distance_to(&a) < 1e-12;
        along_line == from_a
    };
    for fh in brep.faces.keys().collect::<Vec<_>>() {
        let plane = crate::geometry::surface::SurfaceGeometry::Plane {
            origin: a,
            normal,
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let orientation = ccw_face(&brep, fh);
        brep.faces[fh].surface = Some(crate::geometry::surface::SurfaceOnFace {
            geometry: plane,
            orientation,
        });
        brep.recompute_face_pcurves(fh).unwrap();
    }
    assert!(
        brep.validate_all().is_empty(),
        "lamina: {:?}",
        brep.validate_all()
    );
    let face = brep.faces.keys().find(|&fh| ccw_face(&brep, fh)).unwrap();
    let result = brep
        .extrude(face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
        .unwrap();
    assert_eq!(result.top_edge_of_base.len(), 2);
    for &(base, top) in &result.top_edge_of_base {
        let kind = |e: crate::topology::EdgeHandle| {
            std::mem::discriminant(&brep.edges[e].curve.as_ref().unwrap().geometry)
        };
        assert_eq!(kind(base), kind(top), "top edge mirrors its own base edge");
    }
    let kinds: std::collections::HashSet<_> = result
        .top_edge_of_base
        .iter()
        .map(|&(_, top)| std::mem::discriminant(&brep.edges[top].curve.as_ref().unwrap().geometry))
        .collect();
    assert_eq!(kinds.len(), 2, "one straight rim and one arc rim");
    assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
}
