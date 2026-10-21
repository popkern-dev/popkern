//! Far-face carry (local boolean scope stage 2, PK-814,
//! docs/architecture/local-boolean-scope.md).
//!
//! A scoped boolean builds its graph from a scoped operand's ACTIVE
//! faces only; the far faces — whose boxes miss the other operand's, so
//! they carry no record and lie outside it — are kept by a union (either
//! operand's) or a subtraction (`a`'s, stage 4: PK-816) and come through
//! unchanged. They are appended here,
//! entity for entity: surface and orientation, loops, pcurves, edge
//! curves, names, roles and tolerances. A frontier half-edge (on an
//! edge shared with an active face) takes the free twin of the edge the
//! assembly built from that unsplit original edge; far-only vertices and
//! edges are minted once. A frontier edge with no built edge or no free
//! twin (the active side did not come through as an unsplit original)
//! fails the carry, and the caller reruns the boolean unscoped.

use std::collections::HashMap;

use popkern_core::constructors::BRepBuilder;
use popkern_core::topology::{
    BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, ShellHandle, VertexHandle,
};
use slotmap::SecondaryMap;

use crate::merge::graph::{EdgeSource, NodeSource, Operand};
use crate::scope::ActiveScope;
use crate::serve_region::Carried;

/// One operand's original entities already built by the assembly. The
/// edge and vertex maps are dense on that operand's handles: the carry
/// looks them up per half-edge (~1500 at the slot disc's 64th cut).
struct CarryMaps {
    /// The operand whose entities these are.
    operand: Operand,
    /// Unsplit original edges and the result edge built from each.
    edges: SecondaryMap<EdgeHandle, EdgeHandle>,
    /// Original vertices and their result vertices.
    vertices: SecondaryMap<VertexHandle, VertexHandle>,
    /// Result faces built from each face's kept fragments, with their
    /// flip.
    faces: HashMap<FaceHandle, Vec<(FaceHandle, bool)>>,
}

impl CarryMaps {
    /// Empty maps for `operand`.
    fn of(operand: Operand) -> Self {
        Self {
            operand,
            edges: SecondaryMap::new(),
            vertices: SecondaryMap::new(),
            faces: HashMap::new(),
        }
    }

    /// Record a kept fragment's result face.
    fn kept_faces(&mut self, face: crate::merge::graph::FaceRef, result: FaceHandle, flip: bool) {
        if face.operand == self.operand {
            self.faces
                .entry(face.face)
                .or_default()
                .push((result, flip));
        }
    }
}

/// Both operands' carry bookkeeping for one assembly. `b`'s structure
/// counts as carried only when its far faces are (a scoped union,
/// PK-816); otherwise its faces stay touched as before.
pub(super) struct Carry<'s> {
    a: CarryMaps,
    b: CarryMaps,
    far: (Option<&'s ActiveScope>, Option<&'s ActiveScope>),
}

impl<'s> Carry<'s> {
    /// Bookkeeping for the far faces `far` of `a` and `b` that the
    /// operation keeps.
    pub(super) fn new(far: (Option<&'s ActiveScope>, Option<&'s ActiveScope>)) -> Self {
        Self {
            a: CarryMaps::of(Operand::A),
            b: CarryMaps::of(Operand::B),
            far,
        }
    }

    /// Record the result vertex built from a graph node.
    pub(super) fn vertex(&mut self, source: &NodeSource, vh: VertexHandle) {
        match *source {
            NodeSource::OriginalA(v) => {
                self.a.vertices.insert(v, vh);
            }
            NodeSource::OriginalB(v) if self.far.1.is_some() => {
                self.b.vertices.insert(v, vh);
            }
            _ => {}
        }
    }

    /// Record the result edge built from a graph edge.
    pub(super) fn edge(&mut self, source: &EdgeSource, edge: EdgeHandle, carried: &mut Carried) {
        if source.is_unsplit_original() {
            carried.edges.insert(edge, ());
        }
        let (maps, e) = match *source {
            EdgeSource::OriginalA(e) => (&mut self.a, e),
            EdgeSource::OriginalB(e) if self.far.1.is_some() => (&mut self.b, e),
            _ => return,
        };
        maps.edges.insert(e, edge);
        carried.of(maps.operand).0.insert(edge, ());
    }

    /// Record a kept fragment's result face.
    pub(super) fn face(
        &mut self,
        face: crate::merge::graph::FaceRef,
        result: FaceHandle,
        flip: bool,
        carried: &mut Carried,
    ) {
        self.a.kept_faces(face, result, flip);
        self.b.kept_faces(face, result, flip);
        if face.operand == Operand::A || (self.far.1.is_some() && !flip) {
            carried.of(face.operand).1.insert(result, ());
        }
    }

    /// Append the kept far faces of both operands (see [`carry`]).
    pub(super) fn far_faces(
        &mut self,
        builder: &mut BRepBuilder,
        shell: ShellHandle,
        (a, b): (&BRep, &BRep),
        carried: &mut Carried,
    ) -> Result<(), crate::error::BoolDiagnostics> {
        if let Some(far) = self.far.0 {
            carry(builder, shell, a, far, &mut self.a, carried)?;
        }
        if let Some(far) = self.far.1 {
            carry(builder, shell, b, far, &mut self.b, carried)?;
        }
        Ok(())
    }
}

/// [`append_far_faces`] for the assembler: a failure is a merge
/// diagnostic, and every carried edge (an unsplit original) joins
/// `carried`.
fn carry(
    builder: &mut BRepBuilder,
    shell: ShellHandle,
    a: &BRep,
    scope: &ActiveScope,
    maps: &mut CarryMaps,
    carried: &mut Carried,
) -> Result<(), crate::error::BoolDiagnostics> {
    append_far_holes(builder, a, &scope.far_holes, maps)
        .and_then(|()| append_far_faces(builder, shell, a, &scope.far, maps, carried))
        .map_err(|reason| {
            let mut d = crate::error::BoolDiagnostics::default();
            d.warn(
                crate::error::Stage::Merge,
                format!("far-face carry: {reason}"),
            );
            d
        })?;
    carried.edges.extend(maps.edges.values().map(|&e| (e, ())));
    carried
        .of(maps.operand)
        .0
        .extend(maps.edges.values().map(|&e| (e, ())));
    Ok(())
}

/// Append `far` (faces of the operand `a`, in its order) to the assembly.
fn append_far_faces(
    builder: &mut BRepBuilder,
    shell: ShellHandle,
    a: &BRep,
    far: &[FaceHandle],
    maps: &mut CarryMaps,
    carried: &mut Carried,
) -> Result<(), String> {
    for &fh in far {
        let face = a.faces.get(fh).ok_or("far face missing")?;
        let rf = builder.add_face(shell);
        carried.of(maps.operand).1.insert(rf, ());
        {
            let r = builder.brep_mut();
            r.faces[rf].surface = face.surface.clone();
            r.faces[rf].name = face.name;
            r.faces[rf].tolerance = face.tolerance;
            r.faces[rf].role = face.role;
        }
        let loops: Vec<LoopHandle> = std::iter::once(face.outer_loop)
            .chain(face.inner_loops.iter().copied())
            .collect();
        for (k, lh) in loops.into_iter().enumerate() {
            let hes = loop_half_edges(a, lh)?;
            let mut out = Vec::with_capacity(hes.len());
            for h in hes {
                out.push(carry_half_edge(builder, a, h, maps)?);
            }
            if k == 0 {
                builder.add_loop(&out, rf).map_err(|e| e.to_string())?;
            } else {
                builder
                    .add_inner_loop(&out, rf)
                    .map_err(|e| e.to_string())?;
            }
        }
    }
    Ok(())
}

/// Append each active face's far holes (inner loops the graph never
/// saw) to the one result face its kept fragment became: the hole lies
/// away from the tool, so in a union or a subtraction it stays in the
/// kept material around it. A face that came through as several
/// fragments, flipped, or not at all declines the carry.
fn append_far_holes(
    builder: &mut BRepBuilder,
    a: &BRep,
    holes: &std::collections::BTreeMap<FaceHandle, Vec<LoopHandle>>,
    maps: &mut CarryMaps,
) -> Result<(), String> {
    for (face, loops) in holes {
        let target = match maps.faces.get(face).map(Vec::as_slice) {
            Some(&[(rf, false)]) => rf,
            other => return Err(format!("far holes of {face:?}: kept as {other:?}")),
        };
        for &lh in loops {
            let hes = loop_half_edges(a, lh)?;
            let mut out = Vec::with_capacity(hes.len());
            for h in hes {
                out.push(carry_half_edge(builder, a, h, maps)?);
            }
            builder
                .add_inner_loop(&out, target)
                .map_err(|e| e.to_string())?;
        }
    }
    Ok(())
}

fn loop_half_edges(a: &BRep, lh: LoopHandle) -> Result<Vec<HalfEdgeHandle>, String> {
    let start = a.loops.get(lh).ok_or("far loop missing")?.half_edge;
    let mut out = Vec::new();
    let mut he = start;
    for _ in 0..=a.half_edges.len() {
        out.push(he);
        he = a.half_edges.get(he).ok_or("far half-edge missing")?.next;
        if he == start {
            return Ok(out);
        }
    }
    Err("far loop does not close".into())
}

fn carry_vertex(
    builder: &mut BRepBuilder,
    a: &BRep,
    v: VertexHandle,
    maps: &mut CarryMaps,
) -> Result<VertexHandle, String> {
    if let Some(&rv) = maps.vertices.get(v) {
        return Ok(rv);
    }
    let src = a.vertices.get(v).ok_or("far vertex missing")?;
    let rv = builder.add_vertex(src.position);
    let r = builder.brep_mut();
    r.vertices[rv].name = src.name;
    r.vertices[rv].tolerance = src.tolerance;
    maps.vertices.insert(v, rv);
    Ok(rv)
}

/// The result half-edge for `a`'s half-edge `h` (see the module docs).
fn carry_half_edge(
    builder: &mut BRepBuilder,
    a: &BRep,
    h: HalfEdgeHandle,
    maps: &mut CarryMaps,
) -> Result<HalfEdgeHandle, String> {
    let src = a.half_edges.get(h).ok_or("far half-edge missing")?;
    let e = src.edge;
    let origin = carry_vertex(builder, a, src.origin, maps)?;
    let re = match maps.edges.get(e) {
        Some(&re) => re,
        None => {
            // A far-only edge: built once, its first half-edge from the
            // edge's canonical half-edge.
            let edge = a.edges.get(e).ok_or("far edge missing")?;
            let canon = a
                .half_edges
                .get(edge.half_edge)
                .ok_or("far half-edge missing")?;
            let twin = a
                .half_edges
                .get(canon.twin)
                .ok_or("far half-edge missing")?;
            let (v0, v1) = (canon.origin, twin.origin);
            let (rv0, rv1) = (
                carry_vertex(builder, a, v0, maps)?,
                carry_vertex(builder, a, v1, maps)?,
            );
            let (re, _, _) = builder.add_edge(rv0, rv1);
            let r = builder.brep_mut();
            r.edges[re].curve = edge.curve.clone();
            r.edges[re].name = edge.name;
            r.edges[re].tolerance = edge.tolerance;
            r.edges[re].canonical = edge.canonical;
            r.edges[re].role = edge.role;
            maps.edges.insert(e, re);
            re
        }
    };
    let r = builder.brep_mut();
    let first = r.edges[re].half_edge;
    let second = r.half_edges[first].twin;
    let unwired = |x: HalfEdgeHandle| r.half_edges[x].loop_ == LoopHandle::default();
    // A far-only edge pairs by canonical identity (a closed edge has one
    // vertex, so origin cannot tell its sides apart); a frontier edge
    // takes its free side, which must start at this half-edge's vertex.
    let pick =
        if a.edges[e].half_edge == h && unwired(first) && r.half_edges[first].origin == origin {
            first
        } else if unwired(second) && r.half_edges[second].origin == origin {
            second
        } else if unwired(first) && r.half_edges[first].origin == origin {
            first
        } else {
            return Err(format!(
                "frontier edge {e:?} has no free side at its far half-edge's origin"
            ));
        };
    r.half_edges[pick].pcurve = src.pcurve.clone();
    Ok(pick)
}
