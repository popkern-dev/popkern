use super::*;

fn build_flat_hex() -> (BRep, FaceHandle, [Point3; 6]) {
    let (brep, face, pts) = build_flat_ngon(6);
    let hex_points: [Point3; 6] = pts.try_into().expect("6 points");
    (brep, face, hex_points)
}

fn build_flat_ngon(n_he: usize) -> (BRep, FaceHandle, Vec<Point3>) {
    assert!(n_he >= 3, "need at least a triangle");
    let mut brep = BRep::new(Tolerance::default());
    let points: Vec<Point3> = (0..n_he)
        .map(|k| {
            let a = 2.0 * std::f64::consts::PI * (k as f64) / (n_he as f64);
            pt(a.cos(), a.sin(), 0.0)
        })
        .collect();
    let mvfs = brep.mvfs(points[0]).expect("mvfs");
    let first = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, points[1])
        .expect("mev first");
    let mut prev_twin = first.half_edge_twin;
    for p in points.iter().skip(2) {
        prev_twin = brep.mev(prev_twin, *p).expect("mev").half_edge_twin;
    }
    let mef = brep.mef(prev_twin, first.half_edge).expect("mef");
    (brep, mef.face, points)
}

#[test]
fn insert_pns_cap_creates_n_subfaces_with_nurbs_surfaces() {
    use super::super::pns_promote::insert_pns_cap_into_blend_face;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let (mut brep, hex_face, hex_points) = build_flat_hex();
    let (boundary, fillets) = flat_boundary_and_fillets(&hex_points, 3);
    let f_before = brep.faces.len();
    let v_before = brep.vertices.len();

    let result = insert_pns_cap_into_blend_face(&mut brep, hex_face, &boundary, &fillets)
        .expect("insert should succeed");

    assert_eq!(result.sub_faces.len(), 3);
    assert_eq!(brep.vertices.len(), v_before + 1);
    assert_eq!(brep.faces.len(), f_before + 2);

    for sub_face in &result.sub_faces {
        let surf = brep.faces[*sub_face]
            .surface
            .as_ref()
            .expect("sub-face should have a surface assigned")
            .geometry
            .clone();
        match surf {
            SurfaceGeometry::Nurbs(ns) => {
                assert_eq!(ns.degree_u, 4);
                assert_eq!(ns.degree_v, 4);
                assert_eq!(ns.num_u, 5);
                assert_eq!(ns.num_v, 5);
            }
            _ => panic!("expected NurbsSurface on sub-face"),
        }
    }

    let cv = brep
        .vertices
        .get(result.central_vertex)
        .expect("central vertex exists");
    use super::super::pns_adapter::blend_boundary_to_pns_ring_with_surfaces;
    use popkern_core::geometry::gregory::pole_from_boundary;
    let ring = blend_boundary_to_pns_ring_with_surfaces(&boundary, &fillets, &brep).unwrap();
    let expected_pole = pole_from_boundary(&ring);
    assert!(cv.position.distance_to(&expected_pole) < 1e-10);
}

#[test]
fn try_pns_dispatch_succeeds_on_hex_blend() {
    use super::super::pns_promote::try_pns_dispatch;
    use crate::types::CornerMethod;
    use crate::vertex_classifier::VertexBlendClassification;

    let (mut brep, hex_face, hex_points) = build_flat_hex();
    let vertex = VertexHandle::default();
    let (boundary, fillets) = flat_boundary_and_fillets(&hex_points, 3);

    let patches = try_pns_dispatch(
        &mut brep,
        hex_face,
        &boundary,
        &fillets,
        vertex,
        VertexBlendClassification::SameSideConvex,
        Vec::new(),
    )
    .expect("PnS dispatch should succeed on hex blend");

    assert_eq!(patches.len(), 3, "3 sub-faces for a 6-HE hex");
    for patch in &patches {
        assert_eq!(patch.method, CornerMethod::PolyhedralNetSpline);
        assert_eq!(patch.n_sides, 3);
    }
}

#[test]
fn try_pns_dispatch_succeeds_at_higher_valences() {
    use super::super::pns_promote::try_pns_dispatch;
    use crate::types::CornerMethod;
    use crate::vertex_classifier::VertexBlendClassification;

    for &n_fillets in &[5usize, 6] {
        let n_he = 2 * n_fillets;
        let (mut brep, face, ring_points) = build_flat_ngon(n_he);
        let vertex = VertexHandle::default();
        let (boundary, fillets) = flat_boundary_and_fillets_vec(&ring_points, n_fillets);

        let patches = try_pns_dispatch(
            &mut brep,
            face,
            &boundary,
            &fillets,
            vertex,
            VertexBlendClassification::SameSideConvex,
            Vec::new(),
        )
        .unwrap_or_else(|| panic!("N={n_fillets}: PnS dispatch should succeed"));

        assert_eq!(patches.len(), n_fillets);
        for (k, patch) in patches.iter().enumerate() {
            assert_eq!(patch.method, CornerMethod::PolyhedralNetSpline);
            assert_eq!(patch.n_sides, n_fillets, "N={n_fillets} k={k}");
        }
    }
}

fn flat_boundary_and_fillets(
    points: &[Point3; 6],
    n_fillets: usize,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
) {
    flat_boundary_and_fillets_vec(points, n_fillets)
}

fn flat_boundary_and_fillets_vec(
    points: &[Point3],
    n_fillets: usize,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let plane = SurfaceGeometry::Plane {
        origin: pt(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let n_he = 2 * n_fillets;
    let mut segments = Vec::with_capacity(n_he);
    let mut fillets = Vec::new();
    for k in 0..n_he {
        let p0 = points[k];
        let p1 = points[(k + 1) % n_he];
        let curve = line_curve(p0, p1);
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: k,
            curve: curve.clone(),
        });
        fillets.push(crate::fillet_surface::FilletSurfaceResult {
            surface: plane.clone(),
            trim_curve_1: curve.clone(),
            trim_curve_2: curve,
            continuity_downgrade: None,
        });
    }
    let boundary = BlendBoundary {
        vertex: VertexHandle::default(),
        vertex_position: pt(0.0, 0.0, 0.3),
        segments,
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    (boundary, fillets)
}
