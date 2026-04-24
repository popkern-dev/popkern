//! Tests for `pns_promote` — end-to-end BlendBoundary → surfaces.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::pns_promote::promote_blend_boundary_to_pns_surfaces;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn line_curve(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Synthesize a `2·n_fillets`-segment BlendBoundary with all
/// segments as FilletEnd on a z=0 plane fillet surface.
fn synth_planar_blend_for_promote(
    radius: f64,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
    BRep,
) {
    synth_planar_blend_for_promote_n(radius, 3)
}

fn synth_planar_blend_for_promote_n(
    radius: f64,
    n_fillets: usize,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
    BRep,
) {
    let mut segments = Vec::with_capacity(2 * n_fillets);
    let vertex = VertexHandle::default();

    // Plane z = 0.
    let plane = SurfaceGeometry::Plane {
        origin: pt(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let mut fillets = Vec::new();
    let n_seg = 2 * n_fillets;
    for k in 0..n_seg {
        let a0 = 2.0 * std::f64::consts::PI * (k as f64) / (n_seg as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((k + 1) as f64) / (n_seg as f64);
        let p0 = pt(radius * a0.cos(), radius * a0.sin(), 0.0);
        let p1 = pt(radius * a1.cos(), radius * a1.sin(), 0.0);
        let curve = line_curve(p0, p1);
        // All FilletEnd for simplicity — we need a fillet surface
        // index per segment.
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: k,
            curve: curve.clone(),
        });
        let dummy = line_curve(p0, p1);
        fillets.push(crate::fillet_surface::FilletSurfaceResult {
            surface: plane.clone(),
            trim_curve_1: dummy.clone(),
            trim_curve_2: dummy,
            continuity_downgrade: None,
        });
    }

    let boundary = BlendBoundary {
        vertex,
        vertex_position: pt(0.0, 0.0, 0.4),
        segments,
        bridge_fallbacks: Vec::new(),
    };
    let brep = BRep::new(Tolerance::default());
    (boundary, fillets, brep)
}

#[test]
fn promote_produces_n_surfaces_matching_pairs() {
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep)
        .expect("promote should succeed");
    // Each PnS sector covers a (arc, bridge) pair, so surface
    // count = segments.len() / 2.
    assert_eq!(surfaces.len(), bb.segments.len() / 2);
    for (k, surf) in surfaces.iter().enumerate() {
        match surf {
            SurfaceGeometry::Nurbs(ns) => {
                assert_eq!(ns.degree_u, 4, "sector {k} degree_u");
                assert_eq!(ns.degree_v, 4, "sector {k} degree_v");
                assert_eq!(ns.num_u, 5);
                assert_eq!(ns.num_v, 5);
            }
            _ => panic!("sector {k}: expected NurbsSurface, got other variant"),
        }
    }
}

#[test]
fn promoted_surfaces_interpolate_pole_at_canonical_origin() {
    // Sanity: each surface's NURBS point at (u=0, v=0) = canonical
    // patch pole = the auto-derived central point.
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();

    // Derive the auto-pole the same way promote does, for the
    // expected value.
    use super::pns_adapter::blend_boundary_to_pns_ring_with_surfaces;
    use popkern_core::geometry::gregory::pole_from_boundary;
    let ring = blend_boundary_to_pns_ring_with_surfaces(&bb, &fillets, &brep).unwrap();
    let expected_pole = pole_from_boundary(&ring);

    for (k, surf) in surfaces.iter().enumerate() {
        let SurfaceGeometry::Nurbs(ns) = surf else {
            panic!("expected NurbsSurface");
        };
        let p00 = ns.point_at(0.0, 0.0);
        let diff = (p00.x - expected_pole.x).abs()
            + (p00.y - expected_pole.y).abs()
            + (p00.z - expected_pole.z).abs();
        assert!(
            diff < 1e-10,
            "sector {k} pole at (0, 0) = {p00:?}, expected {expected_pole:?}"
        );
    }
}

#[test]
fn promoted_surfaces_boundary_matches_paired_segment_endpoints() {
    // Each surface covers a (arc, bridge) pair. The surface's
    // canonical v_c=4 boundary (v=1 in NURBS params) should go
    // from arc_start (u=0) to bridge_end (u=1).
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();
    for (k, surf) in surfaces.iter().enumerate() {
        let SurfaceGeometry::Nurbs(ns) = surf else {
            panic!("expected NurbsSurface");
        };
        let seg_arc = &bb.segments[2 * k];
        let seg_bridge = &bb.segments[2 * k + 1];
        let arc_curve = match seg_arc {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let bridge_curve = match seg_bridge {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let (ta0, _) = arc_curve.parameter_range();
        let (_, tb1) = bridge_curve.parameter_range();
        let expected_start = arc_curve.point_at(ta0);
        let expected_end = bridge_curve.point_at(tb1);
        let nurbs_u0v1 = ns.point_at(0.0, 1.0);
        let nurbs_u1v1 = ns.point_at(1.0, 1.0);
        assert!(
            (nurbs_u0v1.x - expected_start.x).abs() < 1e-10
                && (nurbs_u0v1.y - expected_start.y).abs() < 1e-10
                && (nurbs_u0v1.z - expected_start.z).abs() < 1e-10,
            "sector {k} boundary start: nurbs (0, 1) = {nurbs_u0v1:?}, expected {expected_start:?}"
        );
        assert!(
            (nurbs_u1v1.x - expected_end.x).abs() < 1e-10
                && (nurbs_u1v1.y - expected_end.y).abs() < 1e-10
                && (nurbs_u1v1.z - expected_end.z).abs() < 1e-10,
            "sector {k} boundary end: nurbs (1, 1) = {nurbs_u1v1:?}, expected {expected_end:?}"
        );
    }
}

/// Build a hexagonal BRep face with 6 half-edges — mimics the
/// shape of a 3-fillet blend boundary on a unit-circle outline.
fn build_flat_hex() -> (BRep, FaceHandle, [Point3; 6]) {
    let (brep, face, pts) = build_flat_ngon(6);
    let hex_points: [Point3; 6] = pts.try_into().expect("6 points");
    (brep, face, hex_points)
}

/// Build a flat regular `n_he`-sided BRep face on a unit circle
/// in z=0. Used to exercise the PnS poke path at valences
/// `N = n_he / 2 ∈ {3, 5, 6, …}`.
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
    use popkern_core::geometry::surface::SurfaceGeometry;

    use super::pns_promote::insert_pns_cap_into_blend_face;

    let (mut brep, hex_face, hex_points) = build_flat_hex();
    let n_fillets = 3;

    // Plane z = 0 serves as the supporting surface for all 6
    // segments.
    let plane = SurfaceGeometry::Plane {
        origin: pt(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let mut segments = Vec::with_capacity(2 * n_fillets);
    let mut fillets = Vec::new();
    for k in 0..6 {
        let p0 = hex_points[k];
        let p1 = hex_points[(k + 1) % 6];
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
    };

    let f_before = brep.faces.len();
    let v_before = brep.vertices.len();

    let result = insert_pns_cap_into_blend_face(&mut brep, hex_face, &boundary, &fillets)
        .expect("insert should succeed");

    // 3 sub-faces from the 2N=6 HE hex → N=3 sub-faces via poke.
    assert_eq!(result.sub_faces.len(), n_fillets);

    // Topology: +1 vertex (central), +2 faces (1 → 3).
    assert_eq!(brep.vertices.len(), v_before + 1);
    assert_eq!(brep.faces.len(), f_before + 2);

    // Each sub-face has a NurbsSurface assigned.
    for sub_face in &result.sub_faces {
        let face_data = brep.faces.get(*sub_face).expect("sub-face exists");
        let surf = face_data
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

    // Central vertex sits at the derived pole.
    let cv = brep
        .vertices
        .get(result.central_vertex)
        .expect("central vertex exists");
    use super::pns_adapter::blend_boundary_to_pns_ring_with_surfaces;
    use popkern_core::geometry::gregory::pole_from_boundary;
    let ring = blend_boundary_to_pns_ring_with_surfaces(&boundary, &fillets, &brep).unwrap();
    let expected_pole = pole_from_boundary(&ring);
    assert!(
        (cv.position.x - expected_pole.x).abs() < 1e-10
            && (cv.position.y - expected_pole.y).abs() < 1e-10
            && (cv.position.z - expected_pole.z).abs() < 1e-10,
        "central vertex {:?} should be at pole {:?}",
        cv.position,
        expected_pole
    );
}

#[test]
fn try_pns_dispatch_succeeds_on_hex_blend() {
    use super::pns_promote::try_pns_dispatch;
    use crate::types::CornerMethod;
    use crate::vertex_classifier::VertexBlendClassification;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let (mut brep, hex_face, hex_points) = build_flat_hex();
    let vertex = VertexHandle::default();

    let plane = SurfaceGeometry::Plane {
        origin: pt(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let mut segments = Vec::with_capacity(6);
    let mut fillets = Vec::new();
    for k in 0..6 {
        let p0 = hex_points[k];
        let p1 = hex_points[(k + 1) % 6];
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
        vertex,
        vertex_position: pt(0.0, 0.0, 0.3),
        segments,
        bridge_fallbacks: Vec::new(),
    };

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
        assert_eq!(patch.method, CornerMethod::Gregory);
        assert_eq!(patch.n_sides, 3);
    }
}

#[test]
fn try_pns_dispatch_succeeds_at_higher_valences() {
    // N=5 and N=6 exercise the paper's general eq-(2) path (N=3
    // is a special-case branch inside `assemble_pns_bi4_cap`).
    // This confirms the full BRep topology pipeline — poke into
    // N sub-faces + per-sub-face NURBS assignment — works at
    // valences beyond the cube-corner trio.
    use super::pns_promote::try_pns_dispatch;
    use crate::types::CornerMethod;
    use crate::vertex_classifier::VertexBlendClassification;
    use popkern_core::geometry::surface::SurfaceGeometry;

    for &n_fillets in &[5usize, 6] {
        let n_he = 2 * n_fillets;
        let (mut brep, face, ring_points) = build_flat_ngon(n_he);
        let vertex = VertexHandle::default();

        let plane = SurfaceGeometry::Plane {
            origin: pt(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let mut segments = Vec::with_capacity(n_he);
        let mut fillets = Vec::new();
        for k in 0..n_he {
            let p0 = ring_points[k];
            let p1 = ring_points[(k + 1) % n_he];
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
            vertex,
            vertex_position: pt(0.0, 0.0, 0.3),
            segments,
            bridge_fallbacks: Vec::new(),
        };

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

        assert_eq!(
            patches.len(),
            n_fillets,
            "N={n_fillets}: expected {n_fillets} sub-faces from poke, got {}",
            patches.len()
        );
        for (k, patch) in patches.iter().enumerate() {
            assert_eq!(patch.method, CornerMethod::Gregory, "N={n_fillets} k={k}");
            assert_eq!(patch.n_sides, n_fillets, "N={n_fillets} k={k}");
        }
    }
}

#[test]
fn cap_boundary_g1_with_planar_support() {
    // When every sector's supporting surface is the same z=0 plane
    // (the synth_planar setup), the bi-4 cap is expected to be
    // coplanar with z=0 — its outer boundary (v=1) must lie in
    // the plane and its surface normal there must be ±z. This
    // verifies the paper's eq (2) G¹ construction produces a cap
    // whose tangent plane matches the support at the boundary.
    //
    // Tested for N ∈ {3, 5, 6, 7}. N=4 is rejected by the cap
    // assembler (σ singularity); N≥5 exercises the general eq-(2)
    // path.
    for n_fillets in [3usize, 5, 6, 7] {
        let (bb, fillets, brep) = synth_planar_blend_for_promote_n(1.5, n_fillets);
        let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();
        assert_eq!(
            surfaces.len(),
            n_fillets,
            "N={n_fillets}: expected {n_fillets} sub-face surfaces"
        );
        let sample_us = [0.0_f64, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0];
        for (k, surf) in surfaces.iter().enumerate() {
            let SurfaceGeometry::Nurbs(ns) = surf else {
                panic!("N={n_fillets} sector {k}: expected NurbsSurface");
            };
            for &u in &sample_us {
                let p = ns.point_at(u, 1.0);
                assert!(
                    p.z.abs() < 1e-9,
                    "N={n_fillets} sector {k} u={u}: boundary point z={} should be 0",
                    p.z
                );
                let n = ns.normal_at(u, 1.0);
                let nl = n.length();
                assert!(
                    nl > 1e-9,
                    "N={n_fillets} sector {k} u={u}: degenerate normal"
                );
                let nz = (n.z / nl).abs();
                assert!(
                    (nz - 1.0).abs() < 1e-6,
                    "N={n_fillets} sector {k} u={u}: normal z-component {nz} should be ±1; \
                     n=({}, {}, {})",
                    n.x,
                    n.y,
                    n.z
                );
            }
        }
    }
}

#[test]
fn cap_boundary_g1_with_tilted_planar_support() {
    // Same shape as the planar test, but rotate the supporting
    // plane + boundary curves so the support normal is NOT aligned
    // with any axis. This catches regressions where the adapter
    // silently emits a z-axis normal regardless of the actual
    // supporting surface.
    let (bb_flat, fillets_flat, brep) = synth_planar_blend_for_promote(1.5);

    // Rotate everything about the x-axis by 30° (unit vectors):
    //   y' = cos θ · y − sin θ · z
    //   z' = sin θ · y + cos θ · z
    let theta = 30_f64.to_radians();
    let (c, s) = (theta.cos(), theta.sin());
    let rot = |p: Point3| Point3::new(p.x, c * p.y - s * p.z, s * p.y + c * p.z);
    let rot_v = |v: Vec3| Vec3::new(v.x, c * v.y - s * v.z, s * v.y + c * v.z);
    let rotate_curve = |curve: &NurbsCurve| NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: curve.control_points.iter().copied().map(rot).collect(),
        weights: curve.weights.clone(),
    };

    let expected_normal = rot_v(Vec3::new(0.0, 0.0, 1.0));

    let tilted_plane = SurfaceGeometry::Plane {
        origin: rot(pt(0.0, 0.0, 0.0)),
        normal: expected_normal,
        u_axis: rot_v(Vec3::new(1.0, 0.0, 0.0)),
        v_axis: rot_v(Vec3::new(0.0, 1.0, 0.0)),
    };

    let segments: Vec<_> = bb_flat
        .segments
        .iter()
        .map(|s| match s {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => BoundarySegment::FilletEnd {
                fillet_index: *fillet_index,
                curve: rotate_curve(curve),
            },
            BoundarySegment::FaceEdge { face, curve } => BoundarySegment::FaceEdge {
                face: *face,
                curve: rotate_curve(curve),
            },
        })
        .collect();
    let bb = BlendBoundary {
        vertex: bb_flat.vertex,
        vertex_position: rot(bb_flat.vertex_position),
        segments,
        bridge_fallbacks: bb_flat.bridge_fallbacks.clone(),
    };
    let fillets: Vec<_> = fillets_flat
        .iter()
        .map(|f| crate::fillet_surface::FilletSurfaceResult {
            surface: tilted_plane.clone(),
            trim_curve_1: rotate_curve(&f.trim_curve_1),
            trim_curve_2: rotate_curve(&f.trim_curve_2),
            continuity_downgrade: f.continuity_downgrade.clone(),
        })
        .collect();

    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();

    let sample_us = [0.1_f64, 0.3, 0.5, 0.7, 0.9];
    for (k, surf) in surfaces.iter().enumerate() {
        let SurfaceGeometry::Nurbs(ns) = surf else {
            panic!("expected NurbsSurface");
        };
        for &u in &sample_us {
            let n = ns.normal_at(u, 1.0);
            let nl = n.length();
            assert!(nl > 1e-9, "sector {k} u={u}: degenerate normal |n|={nl}");
            let unit = Vec3::new(n.x / nl, n.y / nl, n.z / nl);
            let dot = unit.x * expected_normal.x
                + unit.y * expected_normal.y
                + unit.z * expected_normal.z;
            assert!(
                dot.abs() > 1.0 - 1e-6,
                "sector {k} u={u}: cap normal {unit:?} not aligned with plane normal {expected_normal:?} (|dot|={})",
                dot.abs()
            );
        }
    }
}

#[test]
fn cap_boundary_g1_with_cube_corner_cylinder_supports() {
    // Real cube-corner geometry: three r-radius cylinders as fillet
    // ribbons on the +x, +y, +z edges meeting at the origin. Each
    // sector is supported by one cylinder; the three arcs meet
    // pairwise at (0,r,0), (0,0,r), (r,0,0). Splitting each quarter
    // arc in half gives n=3 sectors × 2 segments = 6 segments, all
    // consistent with their cylinder support. At the strip sample
    // points (v=1, u ∈ {0, 0.5, 1} which Bernstein-interpolate to
    // strip.pos[0..=2]) the cap normal must match the cylinder's
    // outer normal — this exercises G¹ with a *curved* support.
    use popkern_core::geometry::surface::SurfaceGeometry;
    let r = 0.3_f64;

    // Cubic Bézier from P_0 to P_3 with prescribed endpoint tangents.
    let cubic_with_tangents = |p0: Point3, t0: Vec3, p3: Point3, t3: Vec3| -> NurbsCurve {
        let p1 = Point3::new(p0.x + t0.x / 3.0, p0.y + t0.y / 3.0, p0.z + t0.z / 3.0);
        let p2 = Point3::new(p3.x - t3.x / 3.0, p3.y - t3.y / 3.0, p3.z - t3.z / 3.0);
        NurbsCurve {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![p0, p1, p2, p3],
            weights: vec![1.0, 1.0, 1.0, 1.0],
        }
    };

    // Arc point + tangent on cylinder with axis through `axis_pt`
    // along `axis_dir`, radius r, at cross-section parameter θ
    // (angle in the plane perpendicular to axis). `basis` gives
    // (perp1, perp2) — orthogonal to axis, forming a right-handed
    // basis. The arc lies in a single perpendicular cross-section.
    let arc_sample = |axis_pt: Point3, perp1: Vec3, perp2: Vec3, theta: f64, dtheta_dt: f64| {
        let cos = theta.cos();
        let sin = theta.sin();
        let p = Point3::new(
            axis_pt.x + r * cos * perp1.x + r * sin * perp2.x,
            axis_pt.y + r * cos * perp1.y + r * sin * perp2.y,
            axis_pt.z + r * cos * perp1.z + r * sin * perp2.z,
        );
        // dP/dθ = -r·sin · perp1 + r·cos · perp2; chain dθ/dt.
        let t = Vec3::new(
            (-r * sin * perp1.x + r * cos * perp2.x) * dtheta_dt,
            (-r * sin * perp1.y + r * cos * perp2.y) * dtheta_dt,
            (-r * sin * perp1.z + r * cos * perp2.z) * dtheta_dt,
        );
        (p, t)
    };

    // Three cylinders (axes along +x, +y, +z, shifted by r into
    // the solid). For each arc, use perp basis (perp1, perp2) so
    // that θ=3π/2 lands on a face-tangent point and θ=π on the
    // neighbouring face-tangent point, tracing the quarter arc on
    // the cube-corner side of each cylinder.
    //
    //  +x fillet: axis (0,r,r)+t·x̂ ; perp1 = ŷ, perp2 = ẑ.
    //      θ=3π/2 → (0, r,  0)  (z=0 tangent)
    //      θ=π    → (0, 0,  r)  (y=0 tangent)
    //  +y fillet: axis (r,0,r)+t·ŷ ; perp1 = ẑ, perp2 = x̂.
    //      θ=3π/2 → (0, 0, r)   (x=0 tangent)
    //      θ=π    → (r, 0, 0)   (z=0 tangent)
    //  +z fillet: axis (r,r,0)+t·ẑ ; perp1 = x̂, perp2 = ŷ.
    //      θ=3π/2 → (r, 0, 0)   (y=0 tangent)
    //      θ=π    → (0, r, 0)   (x=0 tangent)
    let cyl_axes: [(Point3, Vec3, Vec3); 3] = [
        (
            pt(0.0, r, r),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ),
        (
            pt(r, 0.0, r),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
        ),
        (
            pt(r, r, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
    ];
    let cyl_axis_dir: [Vec3; 3] = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];

    // Build the 6 half-arc segments (arc k splits into seg 2k + 2k+1).
    // θ decreases from 3π/2 → 5π/4 (half 1) → π (half 2); dθ/dt = −π/4.
    let theta_start = 1.5 * std::f64::consts::PI;
    let theta_mid = 1.25 * std::f64::consts::PI;
    let theta_end = std::f64::consts::PI;
    let dtheta_dt = -std::f64::consts::FRAC_PI_4;

    let mut segments = Vec::with_capacity(6);
    let mut fillets: Vec<crate::fillet_surface::FilletSurfaceResult> = Vec::with_capacity(6);
    for k in 0..3 {
        let (axis_pt, perp1, perp2) = cyl_axes[k];
        let axis_dir = cyl_axis_dir[k];
        let (p0, t0) = arc_sample(axis_pt, perp1, perp2, theta_start, dtheta_dt);
        let (p1, tm) = arc_sample(axis_pt, perp1, perp2, theta_mid, dtheta_dt);
        let (p2, t2) = arc_sample(axis_pt, perp1, perp2, theta_end, dtheta_dt);
        let half_a = cubic_with_tangents(p0, t0, p1, tm);
        let half_b = cubic_with_tangents(p1, tm, p2, t2);
        let cylinder = SurfaceGeometry::Cylinder {
            origin: axis_pt,
            axis: axis_dir,
            radius: r,
        };
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: 2 * k,
            curve: half_a.clone(),
        });
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: 2 * k + 1,
            curve: half_b.clone(),
        });
        fillets.push(crate::fillet_surface::FilletSurfaceResult {
            surface: cylinder.clone(),
            trim_curve_1: half_a.clone(),
            trim_curve_2: half_a,
            continuity_downgrade: None,
        });
        fillets.push(crate::fillet_surface::FilletSurfaceResult {
            surface: cylinder,
            trim_curve_1: half_b.clone(),
            trim_curve_2: half_b,
            continuity_downgrade: None,
        });
    }

    let boundary = BlendBoundary {
        vertex: VertexHandle::default(),
        vertex_position: pt(0.0, 0.0, 0.0),
        segments,
        bridge_fallbacks: Vec::new(),
    };
    let brep = BRep::new(Tolerance::default());
    let surfaces = promote_blend_boundary_to_pns_surfaces(&boundary, &fillets, &brep)
        .expect("cube-corner promote");
    assert_eq!(surfaces.len(), 3);

    // Sample at u=0.5 — the arc/bridge junction, midway between
    // the two ring-corner spokes. Away from sector corners where
    // twist couples the patch's ∂u and ∂v in ways dependent on
    // neighbours, the cap normal here is purely a function of the
    // sector's own strip data. G¹ with a bi-4 polynomial cap fit
    // against a true circular cylinder still incurs a small
    // polynomial-approximation residual — tolerate up to ~3°.
    for (k, surf) in surfaces.iter().enumerate() {
        let SurfaceGeometry::Nurbs(ns) = surf else {
            panic!("expected NurbsSurface");
        };
        let (axis_pt, _, _) = cyl_axes[k];
        let axis_dir = cyl_axis_dir[k];
        let p = ns.point_at(0.5, 1.0);
        let to_p = Vec3::new(p.x - axis_pt.x, p.y - axis_pt.y, p.z - axis_pt.z);
        let along = to_p.dot(&axis_dir);
        let radial = Vec3::new(
            to_p.x - along * axis_dir.x,
            to_p.y - along * axis_dir.y,
            to_p.z - along * axis_dir.z,
        );
        let rl = radial.length();
        assert!(rl > 1e-9, "sector {k}: p on axis");
        let support_n = Vec3::new(radial.x / rl, radial.y / rl, radial.z / rl);

        let cap_n = ns.normal_at(0.5, 1.0);
        let nl = cap_n.length();
        assert!(nl > 1e-9, "sector {k}: degenerate cap normal");
        let cap_unit = Vec3::new(cap_n.x / nl, cap_n.y / nl, cap_n.z / nl);
        let dot = cap_unit.x * support_n.x + cap_unit.y * support_n.y + cap_unit.z * support_n.z;
        let angle_deg = dot.abs().acos().to_degrees();
        assert!(
            dot.abs() > 0.998,
            "sector {k} u=0.5: cap normal {cap_unit:?} should track \
             cylinder normal {support_n:?} (|dot|={}, angle~{:.3}°)",
            dot.abs(),
            angle_deg,
        );
    }
}

#[test]
fn pns_env_enabled_reads_env_var() {
    use super::pns_promote::pns_env_enabled;
    // Sanity — depends on test-run env. Just check the function
    // doesn't panic and returns a bool.
    let _ = pns_env_enabled();
}

#[test]
fn _unused_faces_stub() {
    let _ = FaceHandle::default();
}
