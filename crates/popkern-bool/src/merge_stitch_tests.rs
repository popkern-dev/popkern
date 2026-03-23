//! Tests for stitch.rs coverage: merge_coincident_vertices, stitch_twin_pair,
//! reconcile_edge_curves, and refit_curve_to_vertices.
#[cfg(test)]
use std::collections::HashSet;

use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

use crate::merge::merge_faces;
use crate::merge::stitch;
use crate::select::SelectionResult;
use crate::split::{FragmentMap, SplitResult};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Build a SplitResult with no actual splitting — both operands are
/// passed through unchanged with identity fragment maps.
fn identity_split(a: BRep, b: BRep) -> SplitResult {
    let fragments_a = FragmentMap {
        origin: a.faces.keys().map(|fh| (fh, fh)).collect(),
        intersection_edges: HashSet::new(),
    };
    let fragments_b = FragmentMap {
        origin: b.faces.keys().map(|fh| (fh, fh)).collect(),
        intersection_edges: HashSet::new(),
    };
    SplitResult {
        a,
        b,
        fragments_a,
        fragments_b,
    }
}

/// Build a SelectionResult that keeps only faces from operand A.
fn select_a_only(split: &SplitResult) -> SelectionResult {
    let keep_a: HashSet<_> = split.a.faces.keys().collect();
    SelectionResult {
        keep_a,
        keep_b: HashSet::new(),
        flip_a: HashSet::new(),
        flip_b: HashSet::new(),
    }
}

// ── stitch.rs coverage: merge_coincident_vertices remap path ──────────────────

/// Build a BRep with two triangular faces using **duplicate** vertex handles
/// for the same positions, then call stitch_boundary_edges directly.
/// This triggers `merge_coincident_vertices` to remap half_edges and remove
/// the duplicate vertices (lines 85, 89, 92, 99, 101-114).
#[test]
fn stitch_boundary_edges_merges_coincident_vertices() {
    use popkern_core::geometry::Point3;
    use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);

    let solid_h = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell_h = brep.shells.insert(Shell {
        face: Default::default(),
        solid: solid_h,
    });
    brep.solids[solid_h].shell = shell_h;

    // Face 1 vertices: v0=(0,0,0), v1=(1,0,0), v2=(0,1,0).
    let v0 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
    });
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 1.0, 0.0),
        outgoing: Default::default(),
    });
    // Face 2 duplicate positions: v3=(0,0,0), v4=(1,0,0), v5=(0,1,0).
    let v3 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
    });
    let v4 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
    });
    let v5 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 1.0, 0.0),
        outgoing: Default::default(),
    });

    let e01 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e12 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e20 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e34 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e45 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e53 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });

    let loop1 = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: Default::default(),
    });
    let he0 = brep.half_edges.insert(HalfEdge {
        origin: v0,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e01,
        loop_: loop1,
        pcurve: None,
    });
    let he1 = brep.half_edges.insert(HalfEdge {
        origin: v1,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e12,
        loop_: loop1,
        pcurve: None,
    });
    let he2 = brep.half_edges.insert(HalfEdge {
        origin: v2,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e20,
        loop_: loop1,
        pcurve: None,
    });
    brep.half_edges[he0].twin = he0;
    brep.half_edges[he1].twin = he1;
    brep.half_edges[he2].twin = he2;
    brep.half_edges[he0].next = he1;
    brep.half_edges[he1].next = he2;
    brep.half_edges[he2].next = he0;
    brep.half_edges[he0].prev = he2;
    brep.half_edges[he1].prev = he0;
    brep.half_edges[he2].prev = he1;
    brep.edges[e01].half_edge = he0;
    brep.edges[e12].half_edge = he1;
    brep.edges[e20].half_edge = he2;
    brep.loops[loop1].half_edge = he0;
    brep.vertices[v0].outgoing = he0;
    brep.vertices[v1].outgoing = he1;
    brep.vertices[v2].outgoing = he2;
    let face1 = brep.faces.insert(Face {
        outer_loop: loop1,
        inner_loops: vec![],
        shell: shell_h,
        surface: None,
    });
    brep.loops[loop1].face = face1;
    brep.shells[shell_h].face = face1;

    let loop2 = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: Default::default(),
    });
    // Face 2: reverse winding so half-edges can be stitched by position.
    // he3: v4→v3 (reverse of v0→v1 after dedup), he4: v5→v4, he5: v3→v5.
    let he3 = brep.half_edges.insert(HalfEdge {
        origin: v4,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e34,
        loop_: loop2,
        pcurve: None,
    });
    let he4 = brep.half_edges.insert(HalfEdge {
        origin: v5,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e45,
        loop_: loop2,
        pcurve: None,
    });
    let he5 = brep.half_edges.insert(HalfEdge {
        origin: v3,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e53,
        loop_: loop2,
        pcurve: None,
    });
    brep.half_edges[he3].twin = he3;
    brep.half_edges[he4].twin = he4;
    brep.half_edges[he5].twin = he5;
    brep.half_edges[he3].next = he4;
    brep.half_edges[he4].next = he5;
    brep.half_edges[he5].next = he3;
    brep.half_edges[he3].prev = he5;
    brep.half_edges[he4].prev = he3;
    brep.half_edges[he5].prev = he4;
    brep.edges[e34].half_edge = he3;
    brep.edges[e45].half_edge = he4;
    brep.edges[e53].half_edge = he5;
    brep.loops[loop2].half_edge = he3;
    brep.vertices[v3].outgoing = he5;
    brep.vertices[v4].outgoing = he3;
    brep.vertices[v5].outgoing = he4;
    let face2 = brep.faces.insert(Face {
        outer_loop: loop2,
        inner_loops: vec![],
        shell: shell_h,
        surface: None,
    });
    brep.loops[loop2].face = face2;

    // 6 vertices before stitching (3 unique positions, each duplicated).
    assert_eq!(brep.vertices.len(), 6);

    stitch::stitch_boundary_edges(&mut brep, &tol);

    // After merging coincident vertices: only 3 unique positions remain.
    assert_eq!(
        brep.vertices.len(),
        3,
        "coincident vertices should be merged"
    );
}

// ── stitch.rs coverage: stitch_twin_pair copies curve from discard_edge ───────

/// When stitch_twin_pair merges two edges and the keep_edge has no curve but
/// the discard_edge does, the curve must be copied to keep_edge (lines 140-143).
#[test]
fn stitch_twin_pair_copies_curve_from_discard_edge() {
    use popkern_core::geometry::{
        Point3,
        curve::{CurveGeometry, CurveOnEdge},
    };
    use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);

    let solid_h = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell_h = brep.shells.insert(Shell {
        face: Default::default(),
        solid: solid_h,
    });
    brep.solids[solid_h].shell = shell_h;

    let v0 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: Default::default(),
    });
    let v1 = brep.vertices.insert(Vertex {
        position: Point3::new(1.0, 0.0, 0.0),
        outgoing: Default::default(),
    });
    let v2 = brep.vertices.insert(Vertex {
        position: Point3::new(0.5, 1.0, 0.0),
        outgoing: Default::default(),
    });

    // keep_edge (from face1) has no curve.
    let e01 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    // discard_edge (from face2) has a line curve.
    let line_curve = CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    };
    let e10 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: Some(line_curve),
    });
    let e12 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e20 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e21 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });
    let e02 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
    });

    let loop1 = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: Default::default(),
    });
    let he_01 = brep.half_edges.insert(HalfEdge {
        origin: v0,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e01,
        loop_: loop1,
        pcurve: None,
    });
    let he_12 = brep.half_edges.insert(HalfEdge {
        origin: v1,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e12,
        loop_: loop1,
        pcurve: None,
    });
    let he_20 = brep.half_edges.insert(HalfEdge {
        origin: v2,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e20,
        loop_: loop1,
        pcurve: None,
    });
    brep.half_edges[he_01].twin = he_01;
    brep.half_edges[he_12].twin = he_12;
    brep.half_edges[he_20].twin = he_20;
    brep.half_edges[he_01].next = he_12;
    brep.half_edges[he_12].next = he_20;
    brep.half_edges[he_20].next = he_01;
    brep.half_edges[he_01].prev = he_20;
    brep.half_edges[he_12].prev = he_01;
    brep.half_edges[he_20].prev = he_12;
    brep.edges[e01].half_edge = he_01;
    brep.edges[e12].half_edge = he_12;
    brep.edges[e20].half_edge = he_20;
    brep.loops[loop1].half_edge = he_01;
    brep.vertices[v0].outgoing = he_01;
    brep.vertices[v1].outgoing = he_12;
    brep.vertices[v2].outgoing = he_20;
    let face1 = brep.faces.insert(Face {
        outer_loop: loop1,
        inner_loops: vec![],
        shell: shell_h,
        surface: None,
    });
    brep.loops[loop1].face = face1;
    brep.shells[shell_h].face = face1;

    // Face 2 uses exact same vertex handles (no duplicate positions needed).
    let loop2 = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: Default::default(),
    });
    // he_10: v1→v0 (reverse of he_01).
    let he_10 = brep.half_edges.insert(HalfEdge {
        origin: v1,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e10,
        loop_: loop2,
        pcurve: None,
    });
    let he_02 = brep.half_edges.insert(HalfEdge {
        origin: v0,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e02,
        loop_: loop2,
        pcurve: None,
    });
    let he_21 = brep.half_edges.insert(HalfEdge {
        origin: v2,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: e21,
        loop_: loop2,
        pcurve: None,
    });
    brep.half_edges[he_10].twin = he_10;
    brep.half_edges[he_02].twin = he_02;
    brep.half_edges[he_21].twin = he_21;
    brep.half_edges[he_10].next = he_02;
    brep.half_edges[he_02].next = he_21;
    brep.half_edges[he_21].next = he_10;
    brep.half_edges[he_10].prev = he_21;
    brep.half_edges[he_02].prev = he_10;
    brep.half_edges[he_21].prev = he_02;
    brep.edges[e10].half_edge = he_10;
    brep.edges[e02].half_edge = he_02;
    brep.edges[e21].half_edge = he_21;
    brep.loops[loop2].half_edge = he_10;
    let face2 = brep.faces.insert(Face {
        outer_loop: loop2,
        inner_loops: vec![],
        shell: shell_h,
        surface: None,
    });
    brep.loops[loop2].face = face2;

    stitch::stitch_boundary_edges(&mut brep, &tol);

    // After stitching he_01 ↔ he_10, the keep_edge (e01) should have the curve
    // copied from the discard_edge (e10).
    let edge_h = brep.half_edges[he_01].edge;
    assert!(
        brep.edges[edge_h].curve.is_some(),
        "keep_edge should have curve from discard_edge"
    );
}

// ── stitch.rs coverage: reconcile_edge_curves missing half_edge guard ─────────

/// Exercise the None guard at line 268-270: edge with a curve but with a
/// default (non-existent) half_edge handle. reconcile_edge_curves skips it.
#[test]
fn stitch_reconcile_handles_missing_half_edge() {
    use popkern_core::geometry::{
        Point3,
        curve::{CurveGeometry, CurveOnEdge},
    };
    use popkern_core::topology::{Edge, HalfEdgeHandle, Shell, Solid, Vertex};

    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);

    let solid_h = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell_h = brep.shells.insert(Shell {
        face: Default::default(),
        solid: solid_h,
    });
    brep.solids[solid_h].shell = shell_h;

    let _v0 = brep.vertices.insert(Vertex {
        position: Point3::new(0.0, 0.0, 0.0),
        outgoing: HalfEdgeHandle::default(),
    });

    // Edge with a curve but a default half_edge handle (non-existent).
    let _orphan_edge = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(0.0, 0.0, 0.0),
                direction: Vec3::new(1.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        }),
    });

    // Should not panic — the missing half_edge guard fires and continues.
    stitch::stitch_boundary_edges(&mut brep, &tol);
}
