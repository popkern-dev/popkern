//! The period rule (plan §4.3a gate): points either side of a seam and
//! points beside a pole are one vertex; two points of a face folding
//! back on itself, close in space, stay two.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::{Incompatibility, ResolutionContext, resolve};
use super::tests::{TOL, cube, options};

/// A point of `face` (operand 0) at chart point `uv`.
fn at_uv(brep: &BRep, face: FaceHandle, uv: (f64, f64), key: u64) -> CandidateVertex {
    let surface = &brep.faces[face].surface.as_ref().expect("surface").geometry;
    CandidateVertex {
        position: surface.point_at(uv.0, uv.1),
        radius: TOL,
        zone: TOL,
        source: SourceKind::EdgeFaceHit,
        incidence: vec![Incidence::Face {
            operand: 0,
            face,
            uv: Some(uv),
        }],
        key: (SourceKind::EdgeFaceHit, 1, key, 0, 0),
    }
}

fn face_with(brep: &BRep, kind: impl Fn(&SurfaceGeometry) -> bool) -> FaceHandle {
    brep.faces
        .iter()
        .find(|(_, f)| f.surface.as_ref().is_some_and(|s| kind(&s.geometry)))
        .map(|(h, _)| h)
        .expect("face")
}

fn groups(brep: &BRep, cands: &[CandidateVertex]) -> super::resolve::Resolution {
    let ctx = ResolutionContext {
        operands: [brep, brep],
        records: &[],
    };
    resolve(cands, &ctx, options(10.0 * TOL)).expect("resolves")
}

/// Either side of a cylinder's seam, a whole period apart in u: one
/// vertex.
#[test]
fn points_either_side_of_a_seam_are_one_vertex() {
    let c = BRep::cylinder(2.0, 3.0, Tolerance::default()).expect("cylinder");
    let wall = face_with(&c, |g| matches!(g, SurfaceGeometry::Cylinder { .. }));
    let tau = 2.0 * std::f64::consts::PI;
    let cands = [
        at_uv(&c, wall, (1e-8, 1.5), 0),
        at_uv(&c, wall, (tau - 1e-8, 1.5), 1),
    ];
    assert!(cands[0].position.distance_to(&cands[1].position) < TOL);
    assert_eq!(groups(&c, &cands).groups.len(), 1);
}

/// Beside a sphere's pole, half a turn apart in u: one vertex.
#[test]
fn points_beside_a_pole_are_one_vertex() {
    let s = BRep::sphere(2.0, Tolerance::default()).expect("sphere");
    let face = face_with(&s, |g| matches!(g, SurfaceGeometry::Sphere { .. }));
    let geom = &s.faces[face].surface.as_ref().expect("surface").geometry;
    let pi = std::f64::consts::PI;
    // The chart's pole: the v where every u meets.
    let pole = [-0.5 * pi, 0.5 * pi, 0.0, pi]
        .into_iter()
        .find(|&v| geom.point_at(0.0, v).distance_to(&geom.point_at(pi, v)) < 1e-12)
        .expect("a pole");
    let off = 0.1 * TOL / 2.0;
    let v = pole - off * pole.signum();
    let cands = [at_uv(&s, face, (0.0, v), 0), at_uv(&s, face, (pi, v), 1)];
    assert!(cands[0].position.distance_to(&cands[1].position) < TOL);
    assert_eq!(groups(&s, &cands).groups.len(), 1);
}

/// A face folding back on itself: its u = 0 and u = 1 edges lie half a
/// tolerance apart in space and 2 apart on the face. Points there are
/// two vertices, as two points of an edge curling back are.
#[test]
fn points_of_a_face_folding_back_on_itself_stay_two() {
    let mut c = cube();
    let top = c
        .faces
        .keys()
        .find(|&f| {
            let ps = c.face_outer_loop_positions(f);
            !ps.is_empty() && ps.iter().all(|p| (p.z - 1.0).abs() < 1e-9)
        })
        .expect("top face");
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let folded = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.5, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            p(0.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(1.0, 1.0, 0.0),
            p(0.0, 0.0, 0.5 * TOL),
            p(0.0, 1.0, 0.5 * TOL),
        ],
        num_u: 3,
        num_v: 2,
        weights: vec![1.0; 6],
    };
    if let Some(s) = c.faces.get_mut(top).and_then(|f| f.surface.as_mut()) {
        s.geometry = SurfaceGeometry::Nurbs(folded);
    }
    let cands = [at_uv(&c, top, (0.0, 0.5), 0), at_uv(&c, top, (1.0, 0.5), 1)];
    assert!(cands[0].position.distance_to(&cands[1].position) < TOL);
    let res = groups(&c, &cands);
    assert_eq!(res.groups.len(), 2, "{res:?}");
    assert!(matches!(
        res.separations.as_slice(),
        [(_, _, Incompatibility::SeparatedAlongFace { length })] if (*length - 2.0).abs() < 1e-6
    ));
}
