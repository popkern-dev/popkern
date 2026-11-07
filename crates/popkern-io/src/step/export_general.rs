//! Exporting a general body (PK-1108, PK-1114).
//!
//! A general body is written as the representation items its parts
//! need, in one representation, the entities they share written once:
//!
//! - the faces bounding solid regions as OCC writes touching solids:
//!   one `CLOSED_SHELL` over every solid's faces, the shared edges and
//!   vertices used by each, a face between two solids listed once, in
//!   one `MANIFOLD_SOLID_BREP`; with cavities a `BREP_WITH_VOIDS`, each
//!   cavity shell (cavities touching each other one shell) an
//!   `ORIENTED_CLOSED_SHELL` written as `export_voids` writes voids;
//! - faces with one region on both sides (sheets in the void, laminae
//!   in a solid or a cavity) in a `SHELL_BASED_SURFACE_MODEL`, an
//!   `OPEN_SHELL` per shell holding them, reusing the solid faces' edges
//!   where they meet, and a closed sheet around a void (a sheet body of
//!   that kind too) as a `CLOSED_SHELL` there;
//! - wire edges as `TRIMMED_CURVE`s (trimmed at their end points) and
//!   acorns as `CARTESIAN_POINT`s in a `GEOMETRIC_CURVE_SET`, as OCC
//!   writes free edges and vertices.
//!
//! OCC reads the items as a compound (its solids are the body's solid
//! regions); popkern reads them back under `StepNonManifold::General`
//! as one general body again (`import_general::merge_touching_items`:
//! the items glued where they coincide, free sheets, wires and acorns
//! placed in the cells holding them).

use std::collections::BTreeMap;

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, FaceHandle, RegionKind, ShellHandle};

use crate::error::IoError;
use crate::step::export::{ExportContext, bool_, enum_, list, ref_, str_param};
use crate::step::parser::types::{EntityId, Parameter};

/// Emit the items of `brep` (a general, wire or acorn body) over its
/// exported faces `faces` (handle and `ADVANCED_FACE` id).
pub(super) fn emit_general(
    ctx: &mut ExportContext,
    brep: &BRep,
    faces: &[(FaceHandle, EntityId)],
) -> Result<Vec<EntityId>, IoError> {
    let region = |s| brep.shell_region(s);
    let solid = |s| region(s).and_then(|r| brep.region_kind(r)) == Some(RegionKind::Solid);
    let mut outer = Vec::new();
    let mut voids: BTreeMap<ShellHandle, Vec<EntityId>> = BTreeMap::new();
    let mut sheets: BTreeMap<ShellHandle, Vec<EntityId>> = BTreeMap::new();
    for &(fh, id) in faces {
        let Some(face) = brep.faces.get(fh) else {
            continue;
        };
        let two_sided = region(face.back) == region(face.front);
        if two_sided || !(solid(face.back) || solid(face.front)) {
            sheets.entry(face.back).or_default().push(id);
        } else if super::export_voids::is_cavity(brep, face.back) {
            voids.entry(face.back).or_default().push(id);
        } else {
            outer.push(id);
        }
    }
    let mut items = Vec::new();
    let closed = |ctx: &mut ExportContext, ids: Vec<EntityId>| {
        let refs: Vec<_> = ids.into_iter().map(ref_).collect();
        ctx.emit("CLOSED_SHELL", vec![str_param(""), list(refs)])
    };
    if !outer.is_empty() {
        let shell = closed(ctx, outer);
        if voids.is_empty() {
            items.push(ctx.emit("MANIFOLD_SOLID_BREP", vec![str_param(""), ref_(shell)]));
        } else {
            let mut void_refs = Vec::with_capacity(voids.len());
            for ids in voids.into_values() {
                let void = closed(ctx, ids);
                void_refs.push(ref_(ctx.emit(
                    "ORIENTED_CLOSED_SHELL",
                    vec![str_param(""), Parameter::Derived, ref_(void), bool_(false)],
                )));
            }
            items.push(ctx.emit(
                "BREP_WITH_VOIDS",
                vec![str_param(""), ref_(shell), list(void_refs)],
            ));
        }
    } else if !voids.is_empty() {
        return Err(IoError::Export(
            "general body has cavity faces but no outer boundary".into(),
        ));
    }
    if !sheets.is_empty() {
        let shells: Vec<_> = sheets
            .into_iter()
            .map(|(shell, ids)| {
                let refs: Vec<_> = ids.into_iter().map(ref_).collect();
                let kind = if is_closed_sheet(brep, shell) {
                    "CLOSED_SHELL"
                } else {
                    "OPEN_SHELL"
                };
                ref_(ctx.emit(kind, vec![str_param(""), list(refs)]))
            })
            .collect();
        items.push(ctx.emit(
            "SHELL_BASED_SURFACE_MODEL",
            vec![str_param(""), list(shells)],
        ));
    }
    let mut elements = Vec::new();
    for edge in brep.wire_edges() {
        elements.push(ref_(emit_wire_curve(ctx, brep, edge)?));
    }
    for shell in brep.shells.values() {
        if let Some(v) = shell.acorn.and_then(|v| brep.vertices.get(v)) {
            elements.push(ref_(ctx.export_point(&v.position)));
        }
    }
    if !elements.is_empty() {
        items.push(ctx.emit("GEOMETRIC_CURVE_SET", vec![str_param(""), list(elements)]));
    }
    Ok(items)
}

/// Whether the sheet faces behind `shell` close up: a closed sheet
/// bounding a void (its faces between the void inside and the region
/// outside), no edge of its faces free. Written as a `CLOSED_SHELL` of
/// the surface model, it reads back as a sheet around a void.
fn is_closed_sheet(brep: &BRep, shell: ShellHandle) -> bool {
    let region = |s| brep.shell_region(s);
    brep.faces.values().filter(|f| f.back == shell).all(|f| {
        region(f.back) != region(f.front)
            && std::iter::once(f.outer_loop)
                .chain(f.inner_loops.iter().copied())
                .flat_map(|l| loop_edges(brep, l))
                .all(|e| !brep.is_free_edge(e))
    })
}

/// The edges of loop `l`, walked from its first fin (bounded by the
/// fin count, so a broken cycle ends).
fn loop_edges(
    brep: &BRep,
    l: popkern_core::topology::LoopHandle,
) -> Vec<popkern_core::topology::EdgeHandle> {
    let mut out = Vec::new();
    let Some(start) = brep.loops.get(l).map(|lp| lp.half_edge) else {
        return out;
    };
    let mut he = start;
    for _ in 0..brep.half_edges.len() {
        let Some(data) = brep.half_edges.get(he) else {
            break;
        };
        out.push(data.edge);
        he = data.next;
        if he == start {
            break;
        }
    }
    out
}

/// A wire edge as a `TRIMMED_CURVE` from its curve's start to its end
/// (a line between its vertices when it has no curve), trimmed at the
/// end points: popkern's curve parameter need not be the written basis
/// curve's (a circle's reference direction is the exporter's), the
/// points are the same. A wire edge's ends are distinct vertices, so
/// the points name the run unambiguously along the curve's sense.
fn emit_wire_curve(
    ctx: &mut ExportContext,
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
) -> Result<EntityId, IoError> {
    let curve = match brep.edges.get(edge).and_then(|e| e.curve.clone()) {
        Some(curve) => curve,
        None => {
            let (a, b) = brep
                .edge_vertices(edge)
                .ok_or_else(|| IoError::Export("wire edge without vertices".into()))?;
            let (p, q) = (brep.vertices[a].position, brep.vertices[b].position);
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: p,
                    direction: q - p,
                },
                t_start: 0.0,
                t_end: 1.0,
            }
        }
    };
    let basis = ctx.export_curve(&curve);
    let start = ctx.export_point(&curve.geometry.point_at(curve.t_start));
    let end = ctx.export_point(&curve.geometry.point_at(curve.t_end));
    Ok(ctx.emit(
        "TRIMMED_CURVE",
        vec![
            str_param(""),
            ref_(basis),
            list(vec![ref_(start)]),
            list(vec![ref_(end)]),
            bool_(curve.t_end >= curve.t_start),
            enum_("CARTESIAN"),
        ],
    ))
}
