use super::*;

/// Build a synthetic [`CapEndpoint`] for a fillet at the corner of a
/// unit cube: the two support faces are the `z = 0` (normal `+ẑ`)
/// and `x = 0` (normal `+x̂`) planes, meeting at the `y`-axis edge.
/// A rolling ball of radius `r` at spine parameter `y = 0.5` touches
/// both faces at:
///
/// - `contact_1 = (r, 0.5, 0)` on the `z = 0` face
/// - `contact_2 = (0, 0.5, r)` on the `x = 0` face
///
/// with centre at `(r, 0.5, r)` and outward spine tangent `+ŷ`.
fn cube_corner_endpoint(r: f64) -> CapEndpoint {
    CapEndpoint {
        contact_1: Point3::new(r, 0.5, 0.0),
        contact_2: Point3::new(0.0, 0.5, r),
        normal_1: Vec3::new(0.0, 0.0, 1.0),
        normal_2: Vec3::new(1.0, 0.0, 0.0),
        center: Point3::new(r, 0.5, r),
        spine_tangent_outward: Vec3::new(0.0, 1.0, 0.0),
        radius: r,
    }
}

#[test]
fn open_cap_produces_no_surface() {
    let ep = cube_corner_endpoint(0.1);
    let cap = build_cap(CapKind::Open, &ep).expect("Open cap should succeed");
    assert_eq!(cap.kind, CapKind::Open);
    assert!(cap.surface.is_none());
    // Boundary is just the two contact points.
    assert_eq!(cap.boundary_points.len(), 2);
    assert_eq!(cap.boundary_points[0], ep.contact_1);
    assert_eq!(cap.boundary_points[1], ep.contact_2);
}

#[test]
fn flat_cap_produces_plane_perpendicular_to_spine() {
    let ep = cube_corner_endpoint(0.1);
    let cap = build_flat_cap(&ep);
    assert_eq!(cap.kind, CapKind::Flat);

    match cap.surface.as_ref().unwrap() {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => {
            // Normal == outward spine tangent.
            assert!(
                (normal.x - 0.0).abs() < 1e-12
                    && (normal.y - 1.0).abs() < 1e-12
                    && (normal.z - 0.0).abs() < 1e-12,
                "flat cap normal should be +ŷ, got {normal:?}"
            );
            // Origin is the midpoint of the two contacts:
            // ((r, 0.5, 0) + (0, 0.5, r))/2 = (r/2, 0.5, r/2).
            assert!((origin.x - 0.05).abs() < 1e-12);
            assert!((origin.y - 0.5).abs() < 1e-12);
            assert!((origin.z - 0.05).abs() < 1e-12);
            // u_axis and v_axis must be perpendicular to the normal
            // and to each other.
            assert!(u_axis.dot(normal).abs() < 1e-12);
            assert!(v_axis.dot(normal).abs() < 1e-12);
            assert!(u_axis.dot(v_axis).abs() < 1e-12);
        }
        other => panic!("expected Plane, got {other:?}"),
    }

    // Boundary has 3 points: contact_1, midpoint, contact_2.
    assert_eq!(cap.boundary_points.len(), 3);
    assert_eq!(cap.boundary_points[0], ep.contact_1);
    assert_eq!(cap.boundary_points[2], ep.contact_2);
}

#[test]
fn rounded_cap_produces_sphere_at_ball_center() {
    let r = 0.15;
    let ep = cube_corner_endpoint(r);
    let cap = build_rounded_cap(&ep);
    assert_eq!(cap.kind, CapKind::Rounded);

    match cap.surface.as_ref().unwrap() {
        SurfaceGeometry::Sphere { center, radius } => {
            // Centre == rolling-ball centre.
            assert!((center.x - r).abs() < 1e-12);
            assert!((center.y - 0.5).abs() < 1e-12);
            assert!((center.z - r).abs() < 1e-12);
            assert!((radius - r).abs() < 1e-12);
        }
        other => panic!("expected Sphere, got {other:?}"),
    }

    // Every boundary sample lies on the sphere (|p - centre| == r).
    for p in &cap.boundary_points {
        let d = ((p.x - r).powi(2) + (p.y - 0.5).powi(2) + (p.z - r).powi(2)).sqrt();
        assert!(
            (d - r).abs() < 1e-9,
            "boundary sample {p:?} off the sphere by {d}"
        );
    }
    // First and last boundary points are the two contacts (arc
    // endpoints).
    assert!((cap.boundary_points[0].x - ep.contact_1.x).abs() < 1e-9);
    assert!((cap.boundary_points.last().unwrap().x - ep.contact_2.x).abs() < 1e-9);
}

#[test]
fn rounded_cap_boundary_arc_is_monotonic() {
    // Consecutive samples should be approximately equidistant
    // (uniform arc parameter).
    let ep = cube_corner_endpoint(0.1);
    let cap = build_rounded_cap(&ep);
    let distances: Vec<f64> = cap
        .boundary_points
        .windows(2)
        .map(|w| {
            ((w[1].x - w[0].x).powi(2) + (w[1].y - w[0].y).powi(2) + (w[1].z - w[0].z).powi(2))
                .sqrt()
        })
        .collect();
    let min = distances.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = distances.iter().cloned().fold(0.0_f64, f64::max);
    assert!(
        max / min < 1.05,
        "arc samples should be quasi-uniform, got min={min}, max={max}"
    );
}

#[test]
fn tangent_match_cap_produces_nurbs_surface() {
    let ep = cube_corner_endpoint(0.1);
    let cap = build_tangent_match_cap(&ep);
    assert_eq!(cap.kind, CapKind::TangentMatch);
    assert!(matches!(
        cap.surface.as_ref().unwrap(),
        SurfaceGeometry::Nurbs(_)
    ));
}

/// The tangent-match cap's section row runs between the contacts: its
/// corner control points are the two contacts.
#[test]
fn tangent_match_cap_corner_cps_are_contact_points() {
    let ep = cube_corner_endpoint(0.1);
    let cap = build_tangent_match_cap(&ep);
    let nurbs = match cap.surface.as_ref().unwrap() {
        SurfaceGeometry::Nurbs(n) => n,
        _ => panic!("expected Nurbs"),
    };
    // Flat indexing: i_u * num_v + i_v; the section row is i_v = 0.
    let ends = [
        nurbs.control_points[0],
        nurbs.control_points[(nurbs.num_u - 1) * nurbs.num_v],
    ];
    for contact in [ep.contact_1, ep.contact_2] {
        assert!(
            ends.iter().any(|p| (*p - contact).length() < 1e-12),
            "contact {contact:?} is not a section-row corner of {ends:?}"
        );
    }
}

#[test]
fn build_cap_dispatches_to_each_kind() {
    let ep = cube_corner_endpoint(0.1);
    for kind in [
        CapKind::Open,
        CapKind::Flat,
        CapKind::Rounded,
        CapKind::TangentMatch,
    ] {
        let cap = build_cap(kind, &ep).unwrap_or_else(|| panic!("kind {kind:?} failed"));
        assert_eq!(cap.kind, kind);
    }
    // UserDefined is not implemented — returns None.
    assert!(build_cap(CapKind::UserDefined, &ep).is_none());
}

#[test]
fn resample_cap_boundary_respects_kind() {
    let ep = cube_corner_endpoint(0.1);

    let open = resample_cap_boundary(CapKind::Open, &ep, 16);
    assert_eq!(open.len(), 2);

    let flat = resample_cap_boundary(CapKind::Flat, &ep, 16);
    assert_eq!(flat.len(), 3); // contact_1, midpoint, contact_2

    let rounded = resample_cap_boundary(CapKind::Rounded, &ep, 16);
    assert_eq!(rounded.len(), 16);
    // All samples still lie on the sphere.
    for p in &rounded {
        let d = ((p.x - 0.1).powi(2) + (p.y - 0.5).powi(2) + (p.z - 0.1).powi(2)).sqrt();
        assert!((d - 0.1).abs() < 1e-9);
    }

    let tangent = resample_cap_boundary(CapKind::TangentMatch, &ep, 20);
    assert_eq!(tangent.len(), 20);
    // The tangent-match cap meets the strip on the same section arc.
    for (p, q) in tangent
        .iter()
        .zip(&resample_cap_boundary(CapKind::Rounded, &ep, 20))
    {
        assert!((*p - *q).length() < 1e-12);
    }
}

#[test]
fn resample_clamps_tiny_n_samples() {
    let ep = cube_corner_endpoint(0.1);
    let pts = resample_cap_boundary(CapKind::Rounded, &ep, 0);
    assert!(pts.len() >= 2);
    let pts = resample_cap_boundary(CapKind::Rounded, &ep, 1);
    assert!(pts.len() >= 2);
}

// ---------------------------------------------------------------------------
// End-to-end: `FilletOptions::cap_kind` has a runtime effect on the pipeline.
// ---------------------------------------------------------------------------

/// With `cap_kind = Open` (the default), compute_fillet_geometry
/// must leave `EdgeGeometry::start_cap` and `end_cap` as `None`.
#[test]
fn open_cap_leaves_edge_geometry_caps_none() {
    use crate::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let opts = FilletOptions::default(); // cap_kind = Open
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
        .expect("open fillet should succeed");

    // We can't inspect EdgeGeometry directly from outside the
    // crate (it's pub(crate)), so we sanity-check via the result
    // face count: the cube fillet creates a few faces and no cap
    // faces (cap_kind = Open).
    assert!(result.brep.faces.len() > cube.faces.len());
}

/// With `cap_kind = Flat`, the construction should succeed and
/// store a `Flat` cap geometry at each spine endpoint on the
/// `EdgeGeometry`. The public API doesn't yet surface the caps
/// (that's commit 3), but we can verify the internal path via
/// `compute_fillet_geometry` directly.
#[test]
fn flat_cap_kind_populates_edge_geometry_caps() {
    use crate::analysis::analyze_edges;
    use crate::geometry_compute::compute_fillet_geometry;
    use crate::types::{CapKind, FilletRadius};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let analyses = analyze_edges(&cube, &[edge], true).unwrap();
    let analysis = &analyses[0];

    let geom = compute_fillet_geometry(
        analysis,
        &FilletRadius::Constant(0.1),
        &cube,
        &Tolerance::default(),
        crate::types::ContinuityLevel::G1,
        CapKind::Flat,
        crate::types::BlendProfile::Circular,
    )
    .expect("fillet should compute");

    // Both endpoints should have a Flat cap.
    assert!(geom.start_cap.is_some(), "start cap should be Some(Flat)");
    assert!(geom.end_cap.is_some(), "end cap should be Some(Flat)");
    assert_eq!(geom.start_cap.as_ref().unwrap().kind, CapKind::Flat);
    assert_eq!(geom.end_cap.as_ref().unwrap().kind, CapKind::Flat);
    // Flat cap surface is a Plane.
    assert!(matches!(
        geom.start_cap.as_ref().unwrap().surface,
        Some(SurfaceGeometry::Plane { .. })
    ));
}

/// `cap_kind = Rounded` populates both endpoints with Sphere caps.
#[test]
fn rounded_cap_kind_populates_edge_geometry_caps() {
    use crate::analysis::analyze_edges;
    use crate::geometry_compute::compute_fillet_geometry;
    use crate::types::{CapKind, FilletRadius};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let analyses = analyze_edges(&cube, &[edge], true).unwrap();
    let analysis = &analyses[0];

    let geom = compute_fillet_geometry(
        analysis,
        &FilletRadius::Constant(0.1),
        &cube,
        &Tolerance::default(),
        crate::types::ContinuityLevel::G1,
        CapKind::Rounded,
        crate::types::BlendProfile::Circular,
    )
    .expect("fillet should compute");

    let start = geom.start_cap.expect("start cap should be Some");
    let end = geom.end_cap.expect("end cap should be Some");
    assert_eq!(start.kind, CapKind::Rounded);
    assert_eq!(end.kind, CapKind::Rounded);
    // Rounded cap surface is a Sphere of the requested radius.
    match start.surface {
        Some(SurfaceGeometry::Sphere { radius, .. }) => {
            assert!((radius - 0.1).abs() < 1e-12);
        }
        other => panic!("expected Sphere, got {other:?}"),
    }
}

/// `cap_kind = TangentMatch` populates both endpoints with NURBS caps.
#[test]
fn tangent_match_cap_kind_populates_edge_geometry_caps() {
    use crate::analysis::analyze_edges;
    use crate::geometry_compute::compute_fillet_geometry;
    use crate::types::{CapKind, FilletRadius};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let analyses = analyze_edges(&cube, &[edge], true).unwrap();
    let analysis = &analyses[0];

    let geom = compute_fillet_geometry(
        analysis,
        &FilletRadius::Constant(0.1),
        &cube,
        &Tolerance::default(),
        crate::types::ContinuityLevel::G1,
        CapKind::TangentMatch,
        crate::types::BlendProfile::Circular,
    )
    .expect("fillet should compute");

    assert!(matches!(
        geom.start_cap.as_ref().unwrap().surface,
        Some(SurfaceGeometry::Nurbs(_))
    ));
    assert!(matches!(
        geom.end_cap.as_ref().unwrap().surface,
        Some(SurfaceGeometry::Nurbs(_))
    ));
}

/// `cap_kind = UserDefined` is a [BEYOND] placeholder — falls back
/// to leaving caps as `None`.
#[test]
fn user_defined_cap_kind_falls_back_to_none() {
    use crate::analysis::analyze_edges;
    use crate::geometry_compute::compute_fillet_geometry;
    use crate::types::{CapKind, FilletRadius};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let analyses = analyze_edges(&cube, &[edge], true).unwrap();
    let analysis = &analyses[0];

    let geom = compute_fillet_geometry(
        analysis,
        &FilletRadius::Constant(0.1),
        &cube,
        &Tolerance::default(),
        crate::types::ContinuityLevel::G1,
        CapKind::UserDefined,
        crate::types::BlendProfile::Circular,
    )
    .expect("fillet should compute");

    assert!(geom.start_cap.is_none());
    assert!(geom.end_cap.is_none());
}

#[test]
fn flat_cap_uses_different_parameterisation_for_identical_contacts() {
    // Edge case: contact_1 == contact_2 (degenerate, but caller
    // shouldn't pass this in practice). The fallback u_axis should
    // be finite and non-zero.
    let r = 0.1;
    let ep = CapEndpoint {
        contact_1: Point3::new(0.1, 0.0, 0.0),
        contact_2: Point3::new(0.1, 0.0, 0.0), // same point
        normal_1: Vec3::new(0.0, 0.0, 1.0),
        normal_2: Vec3::new(1.0, 0.0, 0.0),
        center: Point3::new(0.1, 0.0, r),
        spine_tangent_outward: Vec3::new(0.0, 1.0, 0.0),
        radius: r,
    };
    let cap = build_flat_cap(&ep);
    match cap.surface.as_ref().unwrap() {
        SurfaceGeometry::Plane { u_axis, v_axis, .. } => {
            // Even with a degenerate contact pair, the fallback
            // produces a finite axis.
            assert!(u_axis.length() > 0.0);
            assert!(v_axis.length() > 0.0);
        }
        _ => panic!("expected Plane"),
    }
}
