//! The quadtree lattice mesher: leaf cells of a chart triangulated with
//! conforming closure across refinement transitions (split from
//! `tessellate.rs`). It reads the surface only through `SurfaceEval`, so
//! it meshes a face chart or any regular chart of the atlas alike.

use popkern_core::geometry::{SurfaceEval, Vec3};
use popkern_math::hash::FastHashMap as HashMap;

use crate::adaptive::UvCell;
use crate::mesh::MeshVertex;

/// The leaf `cells` of a chart quadtree as a conforming triangle mesh on
/// `surface`: a cell with no finer neighbour splits along its shorter
/// diagonal in space, a transition cell fans round its centre over its
/// neighbours' corners (hanging nodes), so refinement leaves no crack.
/// Vertex normals face along `orientation`; triangles wind
/// counter-clockwise in the chart. `floor` is the length under which a
/// normal is not normalised.
pub(crate) fn lattice_mesh(
    surface: &dyn SurfaceEval,
    cells: &[UvCell],
    orientation: bool,
    floor: f64,
) -> (Vec<MeshVertex>, Vec<[usize; 3]>) {
    // Deduplicate vertices by UV coordinate using bit-exact keys.
    let mut vertex_map: HashMap<(u64, u64), usize> = HashMap::default();
    let mut vertices: Vec<MeshVertex> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // Conforming closure for the adaptive quadtree (PK-333): a coarse
    // cell adjacent to finer neighbours carries their corner points as
    // HANGING NODES on its edges — triangulating each quad
    // independently left an open crack along every refinement
    // transition (measured: an 80-edge open ring at the cone's
    // near-apex refinement boundary, invisible to the crack zipper,
    // which is scoped to BRep-edge welds). Iso-line registries of all
    // cell corner coordinates let each cell include its neighbours'
    // points; children share bit-exact parent arithmetic, so lookup is
    // exact. Cells with no hanging nodes keep the original
    // shorter-diagonal split bit-for-bit.
    let sortable = |f: f64| -> u64 {
        let b = f.to_bits();
        if b >> 63 == 1 { !b } else { b ^ (1u64 << 63) }
    };
    let mut corners_on_u: HashMap<u64, std::collections::BTreeMap<u64, f64>> = HashMap::default();
    let mut corners_on_v: HashMap<u64, std::collections::BTreeMap<u64, f64>> = HashMap::default();
    for cell in cells {
        for (u, v) in [
            (cell.u_min, cell.v_min),
            (cell.u_max, cell.v_min),
            (cell.u_min, cell.v_max),
            (cell.u_max, cell.v_max),
        ] {
            corners_on_u
                .entry(u.to_bits())
                .or_default()
                .insert(sortable(v), v);
            corners_on_v
                .entry(v.to_bits())
                .or_default()
                .insert(sortable(u), u);
        }
    }
    let hanging = |registry: &HashMap<u64, std::collections::BTreeMap<u64, f64>>,
                   line: f64,
                   lo: f64,
                   hi: f64|
     -> Vec<f64> {
        use std::ops::Bound::Excluded;
        // A zero-extent cell side (a seam-touching sliver trim) has
        // lo == hi; an equal excluded range panics in BTreeMap.
        if sortable(lo) >= sortable(hi) {
            return Vec::new();
        }
        registry
            .get(&line.to_bits())
            .map(|m| {
                m.range((Excluded(sortable(lo)), Excluded(sortable(hi))))
                    .map(|(_, &x)| x)
                    .collect()
            })
            .unwrap_or_default()
    };

    fn intern_vertex(
        u: f64,
        v: f64,
        surface: &dyn SurfaceEval,
        orientation: bool,
        floor: f64,
        vertex_map: &mut HashMap<(u64, u64), usize>,
        vertices: &mut Vec<MeshVertex>,
    ) -> usize {
        let key = (u.to_bits(), v.to_bits());
        *vertex_map.entry(key).or_insert_with(|| {
            let pos = surface.point_at(u, v);
            let mut normal = surface.normal_at(u, v);
            let len = normal.length();
            if len > floor {
                normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
            }
            if !orientation {
                normal = -normal;
            }
            let idx = vertices.len();
            vertices.push(MeshVertex {
                position: pos,
                normal,
                uv: (u, v),
            });
            idx
        })
    }

    for cell in cells {
        let bottom = hanging(&corners_on_v, cell.v_min, cell.u_min, cell.u_max);
        let right = hanging(&corners_on_u, cell.u_max, cell.v_min, cell.v_max);
        let top = hanging(&corners_on_v, cell.v_max, cell.u_min, cell.u_max);
        let left = hanging(&corners_on_u, cell.u_min, cell.v_min, cell.v_max);

        if bottom.is_empty() && right.is_empty() && top.is_empty() && left.is_empty() {
            // Four corners of the quad cell.
            let corners = [
                (cell.u_min, cell.v_min),
                (cell.u_max, cell.v_min),
                (cell.u_min, cell.v_max),
                (cell.u_max, cell.v_max),
            ];
            let mut indices = [0_usize; 4];
            for (i, &(u, v)) in corners.iter().enumerate() {
                indices[i] = intern_vertex(
                    u,
                    v,
                    surface,
                    orientation,
                    floor,
                    &mut vertex_map,
                    &mut vertices,
                );
            }

            // Split quad into 2 triangles along the shorter diagonal.
            // Indices: 0=bottom-left, 1=bottom-right, 2=top-left, 3=top-right
            let diag_02_13 = vertices[indices[0]]
                .position
                .distance_to(&vertices[indices[3]].position);
            let diag_01_23 = vertices[indices[1]]
                .position
                .distance_to(&vertices[indices[2]].position);

            // A chart-degenerate border (a NURBS revolution cap's pole
            // row: every u = 1 sample is the same 3D point) gives the
            // pole-side cell two coincident corners — one of its two
            // triangles is a zero-length-edge sliver. Emitting it
            // leaves a permanent (P, P) half-edge no neighbour can
            // cancel (measured, PK-455: eight per dimple pole, flagged
            // as an open edge by every watertight census). The other
            // triangle of the cell is the honest pole fan — keep it.
            // A cell whose BOTH triangles are degenerate (a 1-cell
            // grid over a bipolar chart) keeps them: an empty cell is
            // a hole, and the caller asked for a mesh this coarse.
            let degenerate = |a: usize, b: usize, c: usize| -> bool {
                let pa = vertices[a].position;
                let pb = vertices[b].position;
                let pc = vertices[c].position;
                let (ab, bc, ca) = (
                    pa.distance_to(&pb),
                    pb.distance_to(&pc),
                    pc.distance_to(&pa),
                );
                let max = ab.max(bc).max(ca);
                ab.min(bc).min(ca) < 1.0e-9 * max
            };
            let pair = if diag_02_13 <= diag_01_23 {
                // Split along 0-3 diagonal.
                [
                    [indices[0], indices[1], indices[3]],
                    [indices[0], indices[3], indices[2]],
                ]
            } else {
                // Split along 1-2 diagonal.
                [
                    [indices[0], indices[1], indices[2]],
                    [indices[1], indices[3], indices[2]],
                ]
            };
            let kept: Vec<[usize; 3]> = pair
                .into_iter()
                .filter(|t| !degenerate(t[0], t[1], t[2]))
                .collect();
            if kept.is_empty() {
                triangles.extend(pair);
            } else {
                triangles.extend(kept);
            }
            continue;
        }

        // Transition cell: fan-triangulate around the cell centre over
        // the counter-clockwise boundary ring including hanging nodes.
        let mut ring: Vec<(f64, f64)> = Vec::new();
        ring.push((cell.u_min, cell.v_min));
        ring.extend(bottom.iter().map(|&u| (u, cell.v_min)));
        ring.push((cell.u_max, cell.v_min));
        ring.extend(right.iter().map(|&v| (cell.u_max, v)));
        ring.push((cell.u_max, cell.v_max));
        ring.extend(top.iter().rev().map(|&u| (u, cell.v_max)));
        ring.push((cell.u_min, cell.v_max));
        ring.extend(left.iter().rev().map(|&v| (cell.u_min, v)));

        let u_mid = 0.5 * (cell.u_min + cell.u_max);
        let v_mid = 0.5 * (cell.v_min + cell.v_max);
        let center = intern_vertex(
            u_mid,
            v_mid,
            surface,
            orientation,
            floor,
            &mut vertex_map,
            &mut vertices,
        );
        let ring_idx: Vec<usize> = ring
            .iter()
            .map(|&(u, v)| {
                intern_vertex(
                    u,
                    v,
                    surface,
                    orientation,
                    floor,
                    &mut vertex_map,
                    &mut vertices,
                )
            })
            .collect();
        // Same zero-length-edge guard as the quad branch (distinct
        // indices can still intern the same 3D point on a
        // chart-degenerate border row), with the same escape hatch: a
        // cell must never end up empty.
        let mut kept: Vec<[usize; 3]> = Vec::new();
        let mut dropped: Vec<[usize; 3]> = Vec::new();
        for i in 0..ring_idx.len() {
            let a = ring_idx[i];
            let b = ring_idx[(i + 1) % ring_idx.len()];
            if a == b {
                continue;
            }
            let (pa, pb, pc) = (
                vertices[center].position,
                vertices[a].position,
                vertices[b].position,
            );
            let (ab, bc, ca) = (
                pa.distance_to(&pb),
                pb.distance_to(&pc),
                pc.distance_to(&pa),
            );
            if ab.min(bc).min(ca) >= 1.0e-9 * ab.max(bc).max(ca) {
                kept.push([center, a, b]);
            } else {
                dropped.push([center, a, b]);
            }
        }
        if kept.is_empty() {
            triangles.extend(dropped);
        } else {
            triangles.extend(kept);
        }
    }

    (vertices, triangles)
}
