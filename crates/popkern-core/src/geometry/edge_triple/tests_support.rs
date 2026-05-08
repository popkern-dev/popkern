use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::EdgeTriple;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

// ─── Fixtures ─────────────────────────────────────────────────────

pub(super) fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

pub(super) fn xz_plane_at(z0: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

pub(super) fn z_axis_cyl(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

pub(super) fn unit_sphere_at(c: Point3) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: c,
        radius: 1.0,
    }
}

pub(super) fn z_axis_cone(half_angle: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle,
    }
}

pub(super) fn nurbs_torus_dummy() -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    }
}

// ─── Hard-guarantee assertions ────────────────────────────────────

pub(super) fn assert_guarantees(
    triple: &EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors_start: Point3,
    anchors_end: Point3,
    tol: &Tolerance,
) {
    // G3 / G5: t_start ≤ t_end and curve.point_at(t_start) is the
    // start anchor.
    assert!(
        triple.t_start <= triple.t_end,
        "G3/G5: t_start={} t_end={}",
        triple.t_start,
        triple.t_end
    );
    let p_start = triple.curve_3d.point_at(triple.t_start);
    let p_end = triple.curve_3d.point_at(triple.t_end);
    assert!(
        p_start.distance_to(&anchors_start) <= tol.linear(),
        "G5 start: {:?} vs {:?}",
        p_start,
        anchors_start
    );
    assert!(
        p_end.distance_to(&anchors_end) <= tol.linear(),
        "G5 end: {:?} vs {:?}",
        p_end,
        anchors_end
    );

    // G1: endpoint CP equality. For NURBS variants, check first/last
    // CPs explicitly. For analytic primitives the endpoint comes from
    // the surface evaluation and is exact-by-construction; the
    // distance check above already covers it.
    if let crate::geometry::curve::CurveGeometry::Nurbs(n) = &triple.curve_3d {
        let cps = &n.control_points;
        assert_eq!(*cps.first().unwrap(), anchors_start, "G1 first CP");
        assert_eq!(*cps.last().unwrap(), anchors_end, "G1 last CP");
    }

    // G2: per-sample positional invariant. Verify already runs 32
    // samples — call again to assert the achieved residual fits.
    let q = triple.verify(surf_a, surf_b, tol, 32).expect("G2 verify");
    assert!(
        q.max_positional_residual <= tol.linear(),
        "G2 residual {:?}",
        q
    );

    // G4: UV-box containment for bounded surfaces is checked
    // inside verify(); if we got here without an error, we're fine.

    // G6: determinism — second invocation of verify is bit-identical.
    let q2 = triple.verify(surf_a, surf_b, tol, 32).expect("G6 verify");
    assert_eq!(
        q.max_positional_residual.to_bits(),
        q2.max_positional_residual.to_bits()
    );
}

// ─── Stage 1: still-passing baseline tests ────────────────────────
