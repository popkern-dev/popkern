//! Which way a pair runs along its span, open and closed.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::EdgeHandle;

use super::runs_same_way;
use crate::merge::graph::{EdgeSource, IntersectionGraph, NodeSource, Operand};

/// The circle of radius 12 about z at height -9, from `t0` to `t1`.
fn circle(t0: f64, t1: f64) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, -9.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 12.0,
        },
        t_start: t0,
        t_end: t1,
    }
}

/// A ring through one node, run each way: its ends cannot tell (both
/// are the node), its leaving tangent can (PK-1133).
#[test]
fn a_ring_run_backwards_is_told_apart_by_its_tangent() {
    let tau = std::f64::consts::TAU;
    let mut g = IntersectionGraph::new();
    let n = g.add_node(Point3::new(12.0, 0.0, -9.0), NodeSource::Intersection);
    let src = || EdgeSource::Original(Operand::A, EdgeHandle::default());
    let (cone, _) = g.add_half_edge_pair(n, n, src(), Some(circle(0.0, tau)));
    let (ball, _) = g.add_half_edge_pair(n, n, src(), Some(circle(tau, 0.0)));
    let (again, _) = g.add_half_edge_pair(n, n, src(), Some(circle(0.0, tau)));
    assert!(!runs_same_way(&g, ball, cone, 1e-6));
    assert!(runs_same_way(&g, again, cone, 1e-6));
    // Negative control: by its ends alone the reversed ring reads the
    // same way.
    let at = |k: usize| g.nodes[g.half_edges[k].origin].position;
    assert!(at(ball).distance_to(&at(cone)) <= 1e-6);
}

/// An open span is decided by its ends, whatever its curve says.
#[test]
fn an_open_span_is_decided_by_its_ends() {
    let mut g = IntersectionGraph::new();
    let a = g.add_node(Point3::new(12.0, 0.0, -9.0), NodeSource::Intersection);
    let b = g.add_node(Point3::new(0.0, 12.0, -9.0), NodeSource::Intersection);
    let src = || EdgeSource::Original(Operand::A, EdgeHandle::default());
    let quarter = std::f64::consts::FRAC_PI_2;
    let (ab, _) = g.add_half_edge_pair(a, b, src(), Some(circle(0.0, quarter)));
    let (ba, _) = g.add_half_edge_pair(b, a, src(), Some(circle(quarter, 0.0)));
    let (ab2, _) = g.add_half_edge_pair(a, b, src(), None);
    assert!(!runs_same_way(&g, ba, ab, 1e-6));
    assert!(runs_same_way(&g, ab2, ab, 1e-6));
}
