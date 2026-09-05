use crate::constructors::builder::BRepBuilder;
use crate::geometry::Point3;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::tolerance::Tolerance;

/// A geometrically u-closed NURBS chart: a degree-1 ring of 4 distinct
/// columns whose first and last control rows coincide.
fn closed_wall() -> SurfaceGeometry {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let ring = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)];
    // Layout: index = i_u * num_v + i_v (u-major).
    let mut cps = Vec::new();
    for (x, y) in ring {
        for z in [0.0, 1.0] {
            cps.push(Point3::new(x, y, z));
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        num_u: 5,
        num_v: 2,
        control_points: cps,
        weights: vec![1.0; 10],
    })
}

#[test]
fn aliased_seam_loop_normalizes_to_a_continuous_contour() {
    // Build a minimal one-face brep whose rectangle loop carries BOTH
    // seam sides at u = 0 (the aliased authorship the boolean
    // assembly's projection recompute produces).
    let mut b = BRepBuilder::new(Tolerance::default());
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let v0 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 1.0));
    let (_, up, down) = b.add_edge(v0, v1);
    let (_, bottom_f, _bottom_r) = b.add_edge(v0, v0);
    let (_, top_f, _top_r) = b.add_edge(v1, v1);
    let face = b.add_face(shell);
    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: closed_wall(),
            orientation: true,
        },
    );
    b.add_loop(&[bottom_f, up, top_f, down], face).unwrap();
    let mut brep = b.build_unchecked();
    // Aliased authorship: both seam halves at u = 0.
    brep.half_edges[bottom_f].pcurve = Some(Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0));
    brep.half_edges[up].pcurve = Some(Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0));
    brep.half_edges[top_f].pcurve = Some(Pcurve::linear((1.0, 1.0), (0.0, 1.0), 0.0, 1.0));
    brep.half_edges[down].pcurve = Some(Pcurve::linear((0.0, 1.0), (0.0, 0.0), 0.0, 1.0));

    brep.normalize_loop_chart_continuity(&[face]);

    // The up-seam must land on u = 1 and every junction must meet.
    let mut he = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut prev_end: Option<Point3> = None;
    for _ in 0..4 {
        let h = &brep.half_edges[he];
        let pc = h.pcurve.as_ref().expect("pcurve");
        let s = pc.point_at(pc.t_start);
        let e = pc.point_at(pc.t_end);
        if let Some(pe) = prev_end {
            assert!(
                (pe.x - s.x).abs() < 1e-9 && (pe.y - s.y).abs() < 1e-9,
                "loop discontinuity: prev end ({},{}) vs start ({},{})",
                pe.x,
                pe.y,
                s.x,
                s.y
            );
        }
        prev_end = Some(e);
        he = h.next;
    }
    let up_pc = brep.half_edges[up].pcurve.as_ref().unwrap();
    assert!(
        (up_pc.point_at(up_pc.t_start).x - 1.0).abs() < 1e-9,
        "up-seam side must sit on u = 1 after normalization"
    );
}

#[test]
fn a_continuous_loop_is_untouched() {
    let mut b = BRepBuilder::new(Tolerance::default());
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let v0 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 1.0));
    let (_, up, down) = b.add_edge(v0, v1);
    let (_, bottom_f, _) = b.add_edge(v0, v0);
    let (_, top_f, _) = b.add_edge(v1, v1);
    let face = b.add_face(shell);
    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: closed_wall(),
            orientation: true,
        },
    );
    b.add_loop(&[bottom_f, up, top_f, down], face).unwrap();
    let mut brep = b.build_unchecked();
    brep.half_edges[bottom_f].pcurve = Some(Pcurve::linear((0.0, 0.0), (1.0, 0.0), 0.0, 1.0));
    brep.half_edges[up].pcurve = Some(Pcurve::linear((1.0, 0.0), (1.0, 1.0), 0.0, 1.0));
    brep.half_edges[top_f].pcurve = Some(Pcurve::linear((1.0, 1.0), (0.0, 1.0), 0.0, 1.0));
    brep.half_edges[down].pcurve = Some(Pcurve::linear((0.0, 1.0), (0.0, 0.0), 0.0, 1.0));

    brep.normalize_loop_chart_continuity(&[face]);

    let up_pc = brep.half_edges[up].pcurve.as_ref().unwrap();
    assert!((up_pc.point_at(up_pc.t_start).x - 1.0).abs() < 1e-12);
    let down_pc = brep.half_edges[down].pcurve.as_ref().unwrap();
    assert!(down_pc.point_at(down_pc.t_start).x.abs() < 1e-12);
}
