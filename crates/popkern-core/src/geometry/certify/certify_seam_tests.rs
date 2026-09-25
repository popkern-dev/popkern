//! Certified seams against analytic references: an exact seam
//! certifies at roundoff, a rolling-ball rail brackets its `1/r`
//! curvature jump, the PK-498 perturbation controls certify above
//! their gates, and every certified bound dominates dense samples.

use super::*;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

fn options() -> SeamCertifyOptions {
    SeamCertifyOptions::for_model(1.0, 1e-6)
}

fn side(surface: &SurfaceGeometry, pcurve: &Pcurve, reversed: bool) -> SeamSideSource {
    SeamSideSource::new(surface, pcurve, reversed).unwrap_or_else(|e| panic!("{e}"))
}

/// Dense-sample maxima of (gap, sine, curvature) from point states.
fn dense(a: &mut SeamSideSource, b: &mut SeamSideSource) -> [f64; 3] {
    let mut out = [0.0_f64; 3];
    for k in 0..=2000 {
        let s = Iv::point(k as f64 / 2000.0);
        let (Ok((sa, _)), Ok((sb, _))) = (a.state(s), b.state(s)) else {
            continue;
        };
        let flip = seam_state::seam_quantities(&sa, &sb, false)
            .normal_dot
            .mid()
            < 0.0;
        let q = seam_state::seam_quantities(&sa, &sb, flip);
        let norm =
            |c: [seam_state::Series; 3]| c.iter().map(|x| x.c[0].mid().powi(2)).sum::<f64>().sqrt();
        out[0] = out[0].max(norm(q.gap));
        out[1] = out[1].max(norm(q.normal_cross));
        out[2] = out[2].max(
            q.curvature
                .iter()
                .map(|x| x.c[0].mid().abs())
                .fold(0.0, f64::max),
        );
    }
    out
}

fn assert_dominates(cert: &SeamCertificate, sampled: [f64; 3], slack: [f64; 3]) {
    let bounds = [cert.g0_upper, cert.g1_upper, cert.g2_upper];
    for k in 0..3 {
        assert!(
            bounds[k] >= sampled[k],
            "bound {k} = {:e} below dense sample {:e}",
            bounds[k],
            sampled[k]
        );
        assert!(
            bounds[k] <= sampled[k] + slack[k],
            "bound {k} = {:e} not tight (dense {:e}, slack {:e})",
            bounds[k],
            sampled[k],
            slack[k]
        );
    }
}

#[test]
fn a_plane_cylinder_rolling_ball_rail_brackets_its_curvature_jump() {
    let r = 0.4;
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, r),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: r,
    };
    let angle = -std::f64::consts::FRAC_PI_2;
    let mut a = side(
        &plane,
        &Pcurve::linear((0.0, 0.0), (2.0, 0.0), 0.0, 1.0),
        false,
    );
    let mut b = side(
        &cylinder,
        &Pcurve::linear((angle, 2.0), (angle, 0.0), 0.0, 1.0),
        true,
    );
    let cert = certify_seam(&mut a, &mut b, &options());
    assert!(cert.is_certified(), "{cert:?}");
    assert!(cert.g0_upper < 1e-13 && cert.g1_upper < 1e-13, "{cert:?}");
    // The transverse curvature jumps by exactly 1/r.
    assert!(
        cert.g2_upper >= 1.0 / r && cert.g2_upper <= 1.0 / r + 1e-3,
        "{cert:?}"
    );
    assert!(
        cert.intervals <= 16,
        "exact rail took {} intervals",
        cert.intervals
    );
    let sampled = dense(&mut a, &mut b);
    assert_dominates(&cert, sampled, [1e-13, 1e-13, 1e-3]);
}

fn rational_patch() -> NurbsSurface {
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..5 {
        for j in 0..4 {
            let (x, y) = (i as f64 * 0.25, j as f64 / 3.0);
            control_points.push(Point3::new(x, y, 0.3 * (x * x - y * y) + 0.1 * x * y));
            weights.push(1.0 + 0.1 * (i as f64 * 0.8).sin() * (j as f64 * 0.6).cos());
        }
    }
    NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points,
        num_u: 5,
        num_v: 4,
        weights,
    }
}

/// The two halves of one rational surface split at `u = 0.35`: an
/// exact seam of every order.
fn split_pair() -> (NurbsSurface, NurbsSurface, Pcurve, Pcurve) {
    let (left, right) = rational_patch()
        .split_at_u(0.35)
        .unwrap_or_else(|e| panic!("{e}"));
    let seam = Pcurve::linear((0.35, 0.0), (0.35, 1.0), 0.0, 1.0);
    (left, right, seam.clone(), seam)
}

#[test]
fn an_exact_rational_seam_certifies_at_roundoff() {
    let (left, right, pa, pb) = split_pair();
    let (sa, sb) = (SurfaceGeometry::Nurbs(left), SurfaceGeometry::Nurbs(right));
    let mut a = side(&sa, &pa, false);
    let mut b = side(&sb, &pb, false);
    let cert = certify_seam(&mut a, &mut b, &options());
    assert!(cert.is_certified(), "{cert:?}");
    // Bounds land within the refinement targets of the true zero.
    let targets = options().targets;
    assert!(cert.g0_upper <= 2.0 * targets.gap, "{cert:?}");
    assert!(cert.g1_upper <= 2.0 * targets.normal_sine, "{cert:?}");
    assert!(cert.g2_upper <= 2.0 * targets.curvature, "{cert:?}");
    assert!(cert.jacobian_lower.iter().all(|j| *j > 0.0));
    // E1 cost gate: an exact seam certifies in few intervals.
    assert!(cert.intervals <= 64, "{} intervals", cert.intervals);
}

/// Move the interior controls of the right half's row `row` (row 0 is
/// the seam) along z.
fn perturbed(surface: &NurbsSurface, row: usize, offset: f64) -> NurbsSurface {
    let mut out = surface.clone();
    for j in 1..out.num_v - 1 {
        out.control_points[row * out.num_v + j].z += offset;
    }
    out
}

#[test]
fn perturbation_controls_certify_above_their_gates_and_dominate_samples() {
    let (left, right, pa, pb) = split_pair();
    let sa = SurfaceGeometry::Nurbs(left);
    // Row 1 tilts the normal; row 2 changes only the curvature.
    for (row, offset) in [(1, 1e-2), (2, 5e-3)] {
        let sb = SurfaceGeometry::Nurbs(perturbed(&right, row, offset));
        let mut a = side(&sa, &pa, false);
        let mut b = side(&sb, &pb, false);
        let cert = certify_seam(&mut a, &mut b, &options());
        assert!(cert.is_certified(), "row {row}: {cert:?}");
        let sampled = dense(&mut a, &mut b);
        let t = options().targets;
        let slack = [
            2.0 * t.gap,
            2.0 * (t.normal_sine + t.relative * sampled[1]),
            2.0 * (t.curvature + t.relative * sampled[2]),
        ];
        assert_dominates(&cert, sampled, slack);
        assert!(
            cert.intervals <= 64,
            "row {row}: {} intervals",
            cert.intervals
        );
        if row == 1 {
            assert!(
                cert.g1_upper > 1e-2,
                "normal tilt certified G1 {:e}",
                cert.g1_upper
            );
        } else {
            let target = options().targets.normal_sine;
            assert!(
                cert.g1_upper <= 2.0 * target,
                "curvature-only control tilted: {cert:?}"
            );
            assert!(
                cert.g2_upper > 1e-2,
                "curvature control certified G2 {:e}",
                cert.g2_upper
            );
        }
    }
}

#[test]
fn a_starved_budget_reports_inconclusive_with_sound_bounds() {
    let (left, right, pa, pb) = split_pair();
    let (sa, sb) = (
        SurfaceGeometry::Nurbs(left),
        SurfaceGeometry::Nurbs(perturbed(&right, 2, 5e-3)),
    );
    let mut a = side(&sa, &pa, false);
    let mut b = side(&sb, &pb, false);
    let mut starved = options();
    starved.budget = 4;
    starved.targets.curvature = 1e-9;
    let cert = certify_seam(&mut a, &mut b, &starved);
    assert!(matches!(
        cert.status,
        CertificateStatus::Inconclusive {
            reason: InconclusiveReason::Budget,
            ..
        }
    ));
    let sampled = dense(&mut a, &mut b);
    assert!(cert.g2_upper >= sampled[2]);
    // Budget-stopped with nothing uncertified: the bounds still cover
    // the seam, only looser than the targets.
    assert!(cert.bounds_cover_seam(), "{cert:?}");
}

#[test]
fn a_seam_through_a_cone_apex_is_inconclusive_at_the_apex() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.5,
    };
    let pcurve = Pcurve::linear((0.3, -1.0), (0.3, 1.0), 0.0, 1.0);
    let mut a = side(&cone, &pcurve, false);
    let mut b = side(&cone, &pcurve, false);
    let cert = certify_seam(&mut a, &mut b, &options());
    let CertificateStatus::Inconclusive { reason, at } = cert.status else {
        panic!("apex seam certified: {cert:?}");
    };
    assert_eq!(reason, InconclusiveReason::Jacobian);
    assert!(at.0 <= 0.5 && 0.5 <= at.1, "apex not located: {at:?}");
    assert!(cert.uncertified_length > 0.0 && cert.uncertified_length < 1e-6);
    assert!(!cert.bounds_cover_seam(), "an apex never covers the seam");
}

/// A seam pcurve given as a degree-1 NURBS polyline with 16 spans
/// (the corpus cylinder rails carry 64): every pre-cut interval still
/// straddles a knot by an ulp in parameter space. United per-branch
/// Taylor forms must certify it like the single-span seam; the hulled
/// natural enclosure alone ran the corpus rails out of budget.
#[test]
fn a_many_span_polyline_pcurve_certifies_through_branch_union() {
    use crate::geometry::nurbs_curve::NurbsCurve;
    let (left, right, pa, _) = split_pair();
    let n = 16;
    let control_points: Vec<Point3> = (0..=n)
        .map(|k| Point3::new(0.35, k as f64 / n as f64, 0.0))
        .collect();
    let mut knots = vec![0.0];
    knots.extend((0..=n).map(|k| k as f64 / n as f64));
    knots.push(1.0);
    let polyline = Pcurve::from_nurbs(NurbsCurve {
        degree: 1,
        knots,
        control_points,
        weights: vec![1.0; n + 1],
    });
    let (sa, sb) = (SurfaceGeometry::Nurbs(left), SurfaceGeometry::Nurbs(right));
    let mut a = side(&sa, &pa, false);
    let mut b = side(&sb, &polyline, false);
    let cert = certify_seam(&mut a, &mut b, &options());
    assert!(cert.is_certified(), "{cert:?}");
    assert_eq!(cert.correspondence, SeamCorrespondence::Fraction);
    assert!(cert.intervals <= 256, "{} intervals", cert.intervals);
}
