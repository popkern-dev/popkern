#![cfg(test)]
use super::*;

#[test]
fn cap_boundary_g1_with_cube_corner_cylinder_supports() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let r = 0.3_f64;

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

    let arc_sample = |axis_pt: Point3, perp1: Vec3, perp2: Vec3, theta: f64, dtheta_dt: f64| {
        let cos = theta.cos();
        let sin = theta.sin();
        let p = Point3::new(
            axis_pt.x + r * cos * perp1.x + r * sin * perp2.x,
            axis_pt.y + r * cos * perp1.y + r * sin * perp2.y,
            axis_pt.z + r * cos * perp1.z + r * sin * perp2.z,
        );
        let t = Vec3::new(
            (-r * sin * perp1.x + r * cos * perp2.x) * dtheta_dt,
            (-r * sin * perp1.y + r * cos * perp2.y) * dtheta_dt,
            (-r * sin * perp1.z + r * cos * perp2.z) * dtheta_dt,
        );
        (p, t)
    };

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
        for curve in [half_a, half_b] {
            fillets.push(crate::fillet_surface::FilletSurfaceResult {
                surface: cylinder.clone(),
                trim_curve_1: curve.clone(),
                trim_curve_2: curve,
                continuity_downgrade: None,
            });
        }
    }

    let boundary = BlendBoundary {
        vertex: VertexHandle::default(),
        vertex_position: pt(0.0, 0.0, 0.0),
        segments,
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    let brep = BRep::new(Tolerance::default());
    let surfaces = promote_blend_boundary_to_pns_surfaces(&boundary, &fillets, &brep)
        .expect("cube-corner promote");
    assert_eq!(surfaces.len(), 3);

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
             cylinder normal {support_n:?} (|dot|={}, angle~{:.3})",
            dot.abs(),
            angle_deg,
        );
    }
}
