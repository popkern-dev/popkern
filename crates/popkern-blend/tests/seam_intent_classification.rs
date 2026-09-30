//! PK-592: structural seam intent in the blend continuity evidence.
//!
//! A blend-to-support seam is intended smooth when the blend runs
//! tangent into the support somewhere along it (every rail does by
//! construction); a strip end that never touches tangency is a designed
//! crease. Before PK-592 a 30° angle bound decided it, and every flat
//! pyramid (strip ends meeting the base at 14–30°) read G0 although every
//! corner and strip seam was G2.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    BlendContinuityEvidence, BlendSeamKind, ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn evidence(
    brep: &BRep,
    pre: &BTreeSet<FaceHandle>,
    level: ContinuityLevel,
) -> BlendContinuityEvidence {
    let blend: Vec<FaceHandle> = brep.faces.keys().filter(|f| !pre.contains(f)).collect();
    measure_blend_continuity_evidence(brep, &blend, level, &ContinuityGates::commercial(), 33)
}

fn census(ev: &BlendContinuityEvidence) -> (usize, usize) {
    let gated = ev.seams.iter().filter(|s| s.gated).count();
    (gated, ev.seams.len() - gated)
}

#[test]
fn flat_pyramid_strip_ends_are_designed_creases() {
    // Measured 2026-09-30 (before: every row read G0, g1 0.26–0.45 on
    // the ends). The gated/sharp census is pinned beside the verdict: a
    // rail collapsing into "sharp" trips it.
    for (n, gated, sharp) in [(3usize, 36usize, 5usize), (5, 60, 9)] {
        let tol = Tolerance::default();
        let body = BRep::n_gonal_pyramid(n, 1.0, 0.25, tol).expect("pyramid");
        let apex = body
            .vertices
            .iter()
            .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
            .map(|(h, _)| h)
            .expect("apex");
        let mut edges: Vec<_> = body
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .collect();
        edges.sort();
        edges.dedup();
        let specs: Vec<FilletSpec> = edges
            .iter()
            .map(|&edge| FilletSpec {
                edge,
                radius: FilletRadius::Constant(0.16),
            })
            .collect();
        let pre: BTreeSet<FaceHandle> = body.faces.keys().collect();
        let options = FilletOptions {
            continuity: ContinuityLevel::G2,
            ..Default::default()
        };
        let result = fillet_edges(&body, &specs, &options, tol, None).expect("fillet");
        let ev = evidence(&result.brep, &pre, ContinuityLevel::G2);
        assert_eq!(ev.achieved, ContinuityLevel::G2, "n={n}: {ev:?}");
        assert_eq!(census(&ev), (gated, sharp), "n={n}: gated/sharp census");
        for seam in ev.seams.iter().filter(|s| !s.gated) {
            if seam.kind == BlendSeamKind::BlendToSupport {
                assert!(
                    seam.best_g1_normal >= 0.2,
                    "n={n}: a designed end never runs tangent (best {:.3})",
                    seam.best_g1_normal
                );
            }
        }
    }
}

/// A cube-edge fillet at `continuity`, its pre-fillet face set, and the
/// strip face.
fn filleted_cube_edge(continuity: ContinuityLevel) -> (BRep, BTreeSet<FaceHandle>, FaceHandle) {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let edge = cube.edges.keys().next().expect("edge");
    let pre: BTreeSet<FaceHandle> = cube.faces.keys().collect();
    let specs = [FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.2),
    }];
    let options = FilletOptions {
        continuity,
        ..Default::default()
    };
    let result = fillet_edges(&cube, &specs, &options, tol, None).expect("fillet");
    let strip = result
        .brep
        .faces
        .keys()
        .find(|f| !pre.contains(f))
        .expect("strip");
    (result.brep, pre, strip)
}

/// Negative control: a rail that degrades along its length (tangent at
/// its ends, creased at mid-length) stays gated and refuses the G1 claim;
/// the structural rule does not excuse it as a designed crease. The G2
/// strip is NURBS, so moving the control points one row in from each
/// boundary at mid-length leaves the rail curves in place and tilts the
/// strip's tangent plane there only.
#[test]
fn a_rail_degrading_along_its_length_stays_gated_and_fails_g1() {
    let (mut brep, pre, strip) = filleted_cube_edge(ContinuityLevel::G2);
    let before = evidence(&brep, &pre, ContinuityLevel::G1);
    assert!(before.achieved >= ContinuityLevel::G1, "{before:?}");
    let (gated_before, _) = census(&before);
    let Some(SurfaceGeometry::Nurbs(net)) =
        brep.faces[strip].surface.as_mut().map(|s| &mut s.geometry)
    else {
        panic!("a G2 strip is NURBS");
    };
    let (nu, nv) = (net.num_u, net.num_v);
    let at = |i: usize, j: usize| i * nv + j;
    let mut targets = vec![(1, nv / 2), (nu - 2, nv / 2), (nu / 2, 1), (nu / 2, nv - 2)];
    targets.dedup();
    for (i, j) in targets {
        let p = |a: usize, b: usize| {
            let q = net.control_points[at(a, b)];
            Vec3::new(q.x, q.y, q.z)
        };
        let normal = (p(i + 1, j) - p(i - 1, j))
            .cross(&(p(i, j + 1) - p(i, j - 1)))
            .normalized()
            .expect("regular net");
        let q = net.control_points[at(i, j)];
        net.control_points[at(i, j)] = Point3::new(q.x, q.y, q.z) + normal * 0.02;
    }
    let after = evidence(&brep, &pre, ContinuityLevel::G1);
    assert!(
        census(&after).0 >= gated_before,
        "the degraded rails stay gated: {after:?}"
    );
    assert!(after.achieved < ContinuityLevel::G1, "{after:?}");
}

/// Negative control: a rail creased uniformly never runs tangent, but its
/// constant sub-bound crease on a curved blend is not a strip end (an end
/// crease varies as the strip's normal sweeps its profile), so it stays
/// gated and fails G1. The support plane is turned rigidly (frame
/// included, so its pcurves keep their meaning) about the rail line.
#[test]
fn a_uniformly_creased_rail_stays_gated_and_fails_g1() {
    let (mut brep, pre, strip) = filleted_cube_edge(ContinuityLevel::G1);
    let (gated_before, sharp_before) = census(&evidence(&brep, &pre, ContinuityLevel::G1));
    // A rail: a strip edge whose other face is a support plane.
    let (support, a, b) = brep
        .half_edges
        .values()
        .filter(|he| brep.loops.get(he.loop_).map(|l| l.face) == Some(strip))
        .find_map(|he| {
            let twin = brep.half_edges.get(he.twin)?;
            let face = brep.loops.get(twin.loop_)?.face;
            let is_plane = matches!(
                brep.faces[face].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            );
            let (a, b) = (
                brep.vertices[he.origin].position,
                brep.vertices[twin.origin].position,
            );
            (is_plane && pre.contains(&face) && a.distance_to(&b) > 0.5).then_some((face, a, b))
        })
        .expect("a rail against a support plane");
    let axis = (b - a).normalized().expect("rail direction");
    let turn = |v: Vec3, angle: f64| {
        // Rodrigues rotation about `axis`.
        v * angle.cos() + axis.cross(&v) * angle.sin() + axis * (axis.dot(&v) * (1.0 - angle.cos()))
    };
    if let Some(SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis,
        v_axis,
    }) = brep.faces[support]
        .surface
        .as_mut()
        .map(|s| &mut s.geometry)
    {
        let angle = 0.2;
        let rel = Vec3::new(origin.x - a.x, origin.y - a.y, origin.z - a.z);
        let moved = turn(rel, angle);
        *origin = Point3::new(a.x + moved.x, a.y + moved.y, a.z + moved.z);
        *normal = turn(*normal, angle);
        *u_axis = turn(*u_axis, angle);
        *v_axis = turn(*v_axis, angle);
    }
    let after = evidence(&brep, &pre, ContinuityLevel::G1);
    assert_eq!(
        census(&after),
        (gated_before, sharp_before),
        "the creased rail stays gated: {after:?}"
    );
    assert!(after.achieved < ContinuityLevel::G1, "{after:?}");
}
