//! Candidate resolution (plan §4.3a gate): separation rules hold as
//! hard constraints, groups anchor at their lowest-dimension member,
//! chains grow within the cap or decline, and the result does not
//! depend on the order the candidates arrive in.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};
use slotmap::Key;

use super::candidates::{CandidateVertex, Incidence, SourceKind, collect_candidates};
use super::resolve::{
    Incompatibility, ResolutionContext, ResolutionDecline, ResolutionOptions, resolve,
};
use crate::intersect::{EndpointBoundaryType, IntersectionRecord};

pub(super) const TOL: f64 = 1e-6;

pub(super) fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("cube")
}

pub(super) fn vertex(operand: usize, vertex: VertexHandle, at: Point3) -> CandidateVertex {
    CandidateVertex {
        position: at,
        radius: TOL,
        zone: TOL,
        source: SourceKind::OperandVertex,
        incidence: vec![Incidence::Vertex { operand, vertex }],
        key: (
            SourceKind::OperandVertex,
            operand,
            vertex.data().as_ffi(),
            0,
            0,
        ),
    }
}

pub(super) fn record_end(record: usize, s: f64, end: bool, at: Point3) -> CandidateVertex {
    CandidateVertex {
        position: at,
        radius: TOL,
        zone: TOL,
        source: SourceKind::RecordEndpoint,
        incidence: vec![Incidence::Record {
            record,
            s,
            end: Some(end),
        }],
        key: (
            SourceKind::RecordEndpoint,
            0,
            record as u64,
            u64::from(end),
            0,
        ),
    }
}

pub(super) fn options(cap: f64) -> ResolutionOptions {
    ResolutionOptions { growth_cap: cap }
}

pub(super) fn quadratic(points: [Point3; 3]) -> NurbsCurve {
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: points.to_vec(),
        weights: vec![1.0; 3],
    }
}

pub(super) fn record(curve_3d: NurbsCurve) -> IntersectionRecord {
    let uv = quadratic([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.5, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ]);
    IntersectionRecord {
        face_a: Default::default(),
        face_b: Default::default(),
        curve: CurveOnSurface {
            curve_3d,
            curve_uv1: uv.clone(),
            curve_uv2: uv,
            is_closed: false,
        },
        is_tangent: false,
        is_coincident: false,
        start_boundary: EndpointBoundaryType::OnEdgeA,
        end_boundary: EndpointBoundaryType::OnEdgeA,
        residue_owner_a: None,
    }
}

/// Two vertices of one operand closer than their tolerances stay two:
/// a valid operand's vertices are distinct.
#[test]
fn two_vertices_of_one_operand_stay_two() {
    let c = cube();
    let vs: Vec<VertexHandle> = c.vertices.keys().take(2).collect();
    let cands = vec![
        vertex(0, vs[0], Point3::new(0.0, 0.0, 0.0)),
        vertex(0, vs[1], Point3::new(0.5 * TOL, 0.0, 0.0)),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &[],
    };
    let res = resolve(&cands, &ctx, options(20.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 2);
    assert_eq!(res.components, 1);
    assert!(matches!(
        res.separations.as_slice(),
        [(0, 1, Incompatibility::DistinctOperandVertices)]
    ));
}

/// Negative control for the rule above: the same two points from
/// different operands are one vertex.
#[test]
fn vertices_of_different_operands_within_tolerance_are_one() {
    let c = cube();
    let vs: Vec<VertexHandle> = c.vertices.keys().take(2).collect();
    let cands = vec![
        vertex(0, vs[0], Point3::new(0.0, 0.0, 0.0)),
        vertex(1, vs[1], Point3::new(0.5 * TOL, 0.0, 0.0)),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &[],
    };
    let res = resolve(&cands, &ctx, options(20.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 1);
    assert_eq!(res.groups[0].anchor, 0);
}

/// A record curling back to within a tolerance of its own start keeps
/// both ends: they are 2 apart along it.
#[test]
fn a_record_curling_back_on_itself_keeps_both_ends() {
    let hairpin = record(quadratic([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 0.5 * TOL, 0.0),
    ]));
    let c = cube();
    let records = [hairpin];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let cands = vec![
        record_end(0, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
        record_end(0, 1.0, true, Point3::new(0.0, 0.5 * TOL, 0.0)),
    ];
    let res = resolve(&cands, &ctx, options(20.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 2);
    assert!(matches!(
        res.separations.as_slice(),
        [(0, 1, Incompatibility::SeparatedAlongRecord { arc })] if *arc > 1.0
    ));
}

/// A record end beside an operand vertex anchors at the vertex, which
/// precedes it, whichever order they arrive in.
#[test]
fn a_record_end_beside_a_vertex_anchors_at_the_vertex() {
    let c = cube();
    let v = c.vertices.keys().next().expect("vertex");
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &[],
    };
    let mut cands = vec![
        record_end(0, 0.0, false, Point3::new(0.6 * TOL, 0.0, 0.0)),
        vertex(0, v, Point3::new(0.0, 0.0, 0.0)),
    ];
    cands.sort_by_key(|c| c.key);
    let res = resolve(&cands, &ctx, options(20.0 * TOL)).expect("resolves");
    assert_eq!(res.groups.len(), 1);
    let anchor = &cands[res.groups[0].anchor];
    assert_eq!(anchor.source, SourceKind::OperandVertex);
    assert_eq!(anchor.position, Point3::new(0.0, 0.0, 0.0));
    // The vertex tolerance covers the record end's ball.
    assert!((res.groups[0].tolerance - 1.6 * TOL).abs() < 1e-12);
}

/// A chain of candidates each 1.8 tolerances from the next: every
/// neighbour pair overlaps. Short, it grows one recorded vertex; long,
/// past the cap, it declines instead of welding silently.
#[test]
fn a_chain_grows_within_the_cap_and_declines_past_it() {
    let c = cube();
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &[],
    };
    let chain = |n: usize| -> Vec<CandidateVertex> {
        (0..n)
            .map(|k| record_end(k, 0.0, false, Point3::new(1.8 * TOL * k as f64, 0.0, 0.0)))
            .collect()
    };
    let short = resolve(&chain(3), &ctx, options(10.0 * TOL)).expect("grows");
    assert_eq!(short.groups.len(), 1);
    assert!((short.groups[0].tolerance - 4.6 * TOL).abs() < 1e-12);
    match resolve(&chain(8), &ctx, options(10.0 * TOL)) {
        Err(ResolutionDecline::ToleranceGrowthExceeded { needed, cap, .. }) => {
            assert!(needed > cap);
        }
        other => panic!("a long chain must decline: {other:?}"),
    }
}

/// The resolution of a real boolean's candidates does not depend on
/// the order they arrive in: shuffled, re-sorted by key, the vertex
/// table is the same.
#[test]
fn resolution_is_independent_of_candidate_order() {
    let tol = Tolerance::default();
    let a = BRep::unit_cube(tol).expect("a");
    let b = BRep::unit_cube(tol)
        .expect("b")
        .transform(&popkern_core::geometry::transform::Transform3::translation(
            0.5, 0.25, 0.0,
        ))
        .expect("b placed");
    let ix = crate::intersect::find_intersections(&a, &b, &tol).expect("intersections");
    let cands = collect_candidates(&a, &b, &ix, &ix.records, tol.linear(), 20.0 * tol.linear());
    let ctx = ResolutionContext {
        operands: [&a, &b],
        records: &ix.records,
    };
    let base = resolve(&cands, &ctx, options(20.0 * tol.linear())).expect("resolves");
    let mut shuffled = cands.clone();
    shuffled.reverse();
    shuffled.rotate_left(cands.len() / 3);
    shuffled.sort_by_key(|c| c.key);
    let again = resolve(&shuffled, &ctx, options(20.0 * tol.linear())).expect("resolves");
    let anchors = |r: &super::resolve::Resolution, c: &[CandidateVertex]| -> Vec<_> {
        r.groups
            .iter()
            .map(|g| (c[g.anchor].key, g.members.len()))
            .collect()
    };
    assert_eq!(anchors(&base, &cands), anchors(&again, &shuffled));
    // Flush faces: the coincident corners resolve onto A's vertices.
    assert!(base.components > 0);
}

pub(super) fn line(a: Point3, b: Point3) -> NurbsCurve {
    quadratic([
        a,
        Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z)),
        b,
    ])
}

/// A record ending mid-span of another record lies on its interior: a
/// T-junction, the pave today's split pass inserts after the fact. A
/// record ending a thousandth short of it does not.
#[test]
fn a_record_ending_on_another_record_is_a_t_junction() {
    // Operands far from the records: only the records interfere.
    let c = cube()
        .transform(&popkern_core::geometry::transform::Transform3::translation(
            100.0, 0.0, 0.0,
        ))
        .expect("cube away");
    // Two face pairs: one pair's records are branches of one
    // intersection, which the branch rule keeps apart.
    let mut faces = c.faces.keys();
    let (fa, fb) = (faces.next().expect("face"), faces.next().expect("face"));
    for (gap, junction) in [(0.0, true), (0.5 * TOL, true), (1e-3, false)] {
        let mut records = [
            record(line(Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0))),
            record(line(Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, gap, 0.0))),
        ];
        records[0].face_a = fa;
        records[1].face_a = fb;
        let mut cands = vec![
            record_end(0, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
            record_end(0, 1.0, true, Point3::new(2.0, 0.0, 0.0)),
            record_end(1, 0.0, false, Point3::new(1.0, 1.0, 0.0)),
            record_end(1, 1.0, true, Point3::new(1.0, gap, 0.0)),
        ];
        super::incidence::add_interior_incidences(&mut cands, [&c, &c], &records, TOL);
        let found: Vec<bool> = cands.iter().map(super::incidence::is_t_junction).collect();
        assert_eq!(found, [false, false, false, junction], "gap {gap}");
        if junction {
            assert!(cands[3].incidence.iter().any(|i| matches!(
                i,
                Incidence::Record { record: 0, s, end: None } if (s - 0.5).abs() < 1e-6
            )));
        }
    }
}

/// Two computations of one crossing (an edge-face hit and a record
/// end) 1.6e-5 apart along an edge: one vertex where the edge crosses
/// the other operand's face at 5°, two where it crosses at 60°. The
/// zone is the radius over sin θ along the edge, so a shallow crossing
/// links them and a steep one does not.
#[test]
fn a_shallow_crossing_links_its_two_computations_along_the_edge() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let tol = Tolerance::default();
    let a = cube();
    let (edge, _) = a
        .edges
        .iter()
        .find(|(_, e)| {
            e.curve.as_ref().is_some_and(|c| {
                let (p, q) = (c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end));
                p.y.abs() + p.z.abs() + q.y.abs() + q.z.abs() < 1e-12
            })
        })
        .expect("the x edge at y = z = 0");
    for (theta_deg, one) in [(5.0_f64, true), (60.0, false)] {
        // B's face normal makes 90° − θ with the edge, so the edge
        // crosses its plane at θ.
        let phi = (90.0 - theta_deg).to_radians();
        let b = BRep::unit_cube(tol)
            .expect("b")
            .transform(&popkern_core::geometry::transform::Transform3::rotation_z(
                phi,
            ))
            .expect("b turned");
        let want = (phi.cos(), phi.sin());
        let (face, _) = b
            .faces
            .iter()
            .find(|(_, f)| {
                f.surface.as_ref().is_some_and(|s| {
                    let n = s.geometry.normal_at(0.0, 0.0);
                    (n.x - want.0).abs() + (n.y - want.1).abs() + n.z.abs() < 1e-9
                })
            })
            .expect("the turned x face");
        let on = |x: f64, source: SourceKind, extra: Incidence| CandidateVertex {
            position: Point3::new(x, 0.0, 0.0),
            radius: 2.0 * TOL,
            zone: 2.0 * TOL,
            source,
            incidence: vec![
                Incidence::Edge {
                    operand: 0,
                    edge,
                    t: x,
                },
                Incidence::Face {
                    operand: 1,
                    face,
                    uv: None,
                },
                extra,
            ],
            key: (source, 0, x.to_bits(), 0, 0),
        };
        let mut cands = vec![
            on(
                0.5,
                SourceKind::EdgeFaceHit,
                Incidence::Face {
                    operand: 1,
                    face,
                    uv: None,
                },
            ),
            on(
                0.5 + 1.6e-5,
                SourceKind::RecordEndpoint,
                Incidence::Record {
                    record: 0,
                    s: 0.0,
                    end: Some(false),
                },
            ),
        ];
        super::incidence::add_crossing_zones(&mut cands, [&a, &b], 20.0 * TOL);
        let ctx = ResolutionContext {
            operands: [&a, &b],
            records: &[],
        };
        let res = resolve(&cands, &ctx, options(20.0 * TOL)).expect("resolves");
        assert_eq!(
            res.groups.len() == 1,
            one,
            "θ = {theta_deg}°: {:?}",
            res.groups
        );
    }
}

/// Endpoint closure: a record end alone in the face interior is an
/// open intersection loop; beside an operand vertex it closes, and a
/// closed record's two ends close each other.
#[test]
fn a_record_end_off_every_boundary_is_open() {
    let c = cube();
    let v = c.vertices.keys().next().expect("vertex");
    let records = [record(line(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
    ))];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let open = |cands: &[CandidateVertex]| -> Vec<bool> {
        let res = resolve(cands, &ctx, options(20.0 * TOL)).expect("resolves");
        (0..cands.len())
            .filter(|&i| cands[i].source == SourceKind::RecordEndpoint)
            .map(|i| !super::shadow::end_closes(cands, &res, i))
            .collect()
    };
    let lone = vec![
        record_end(0, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
        record_end(0, 1.0, true, Point3::new(1.0, 0.0, 0.0)),
    ];
    assert_eq!(open(&lone), [true, true]);
    let mut at_vertex = vec![vertex(0, v, Point3::new(0.0, 0.0, 0.0))];
    at_vertex.extend(lone.iter().cloned());
    assert_eq!(open(&at_vertex), [false, true]);
    let loop_records = [record(quadratic([
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    ]))];
    let loop_ctx = ResolutionContext {
        operands: [&c, &c],
        records: &loop_records,
    };
    let closed = vec![
        record_end(0, 0.0, false, Point3::new(0.0, 0.0, 0.0)),
        record_end(0, 1.0, true, Point3::new(0.0, 0.0, 0.0)),
    ];
    let res = resolve(&closed, &loop_ctx, options(20.0 * TOL)).expect("resolves");
    assert!((0..2).all(|i| super::shadow::end_closes(&closed, &res, i)));
}
