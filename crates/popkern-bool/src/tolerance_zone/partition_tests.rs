//! The exact partition (plan §4.3a gate): a component the
//! canonical-order pass cannot prove minimal is searched for its fewest
//! vertices, or declined when too large to search.

use popkern_core::geometry::Point3;
use popkern_core::topology::VertexHandle;

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::{ResolutionContext, ResolutionDecline, ResolutionOptions, resolve};
use super::tests::{TOL, cube, line, record, vertex};

fn options(exact_limit: usize) -> ResolutionOptions {
    ResolutionOptions {
        growth_cap: 10.0 * TOL,
        exact_limit,
    }
}

/// Two vertices a, b of one operand 1.8 tol apart; a point c between
/// them touching both; and a vertex d of the other operand that a, also
/// coinciding with another vertex of that operand, may not join,
/// touching only c. The canonical-order pass puts c with a, so d starts
/// a third vertex; the fewest is two, {a} and {b, c, d}. (Built on the
/// retired branch rule until PK-1130; the incompatibility is now two
/// vertices of one operand.)
fn early_join() -> Vec<CandidateVertex> {
    let c = cube();
    let vs: Vec<VertexHandle> = c.vertices.keys().take(4).collect();
    let mut a = vertex(0, vs[0], Point3::new(0.0, 0.0, 0.0));
    a.incidence.push(Incidence::Vertex {
        operand: 1,
        vertex: vs[2],
    });
    let b = vertex(0, vs[1], Point3::new(1.8 * TOL, 0.0, 0.0));
    let mid = CandidateVertex {
        position: Point3::new(0.9 * TOL, 0.0, 0.0),
        radius: TOL,
        zone: TOL,
        source: SourceKind::EdgeFaceHit,
        incidence: vec![],
        key: (SourceKind::EdgeFaceHit, 0, 0, 0, 0),
    };
    let d = vertex(1, vs[3], Point3::new(-0.5 * TOL, 1.2 * TOL, 0.0));
    vec![a, b, mid, d]
}

fn records() -> [crate::intersect::IntersectionRecord; 2] {
    [
        record(line(
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        )),
        record(line(
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(-0.5 * TOL, 1.2 * TOL, 0.0),
        )),
    ]
}

#[test]
fn an_early_join_costing_a_vertex_is_searched_out() {
    let c = cube();
    let records = records();
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let cands = early_join();
    let res = resolve(&cands, &ctx, options(12)).expect("resolves");
    assert_eq!(res.groups.len(), 2, "{res:?}");
    assert_eq!(res.exact_partitions, 1);
    assert_eq!(res.groups[res.group_of[1]].anchor, 1, "{res:?}");
    assert_eq!(res.group_of[2], res.group_of[1], "{res:?}");
    assert_eq!(res.group_of[3], res.group_of[1], "{res:?}");
    assert_ne!(res.group_of[0], res.group_of[1], "{res:?}");
}

/// The same component over the search limit: the pass made three
/// vertices against a lower bound of two, so it declines rather than
/// guess.
#[test]
fn a_component_too_large_to_search_declines_when_not_provably_minimal() {
    let c = cube();
    let records = records();
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    match resolve(&early_join(), &ctx, options(0)) {
        Err(ResolutionDecline::CandidateResolutionAmbiguous {
            component: 0,
            members,
            vertices: 3,
            lower_bound: 2,
        }) => assert_eq!(members, vec![0, 1, 2, 3]),
        other => panic!("{other:?}"),
    }
}

/// Over the search limit, a pass that reaches the lower bound is
/// minimal by proof: two vertices of one operand stay two, no search.
#[test]
fn a_pass_reaching_the_lower_bound_is_minimal_without_search() {
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
    let res = resolve(&cands, &ctx, options(0)).expect("minimal by its bound");
    assert_eq!(res.groups.len(), 2);
    assert_eq!(res.exact_partitions, 0);
}
