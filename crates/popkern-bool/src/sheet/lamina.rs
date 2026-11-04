//! One sheet form for the cell complex (PK-978, design §2): a face
//! holding both half-edges of each boundary edge, the outer loop and
//! the same edges reversed as an inner loop (`BodyKind::Sheet`).
//!
//! The core's sheet constructors build the other form, a front face
//! and a surfaceless backside face sharing the boundary edges (a closed
//! shell that reads as `BodyKind::Solid`). The backside face's loop
//! already traverses the boundary reversed: moved onto the front face
//! as an inner loop, it is the single-face form. Pcurves follow on
//! planar fronts; a curved front keeps its own and the moved half-edges
//! take theirs from the front's twins.

use popkern_core::topology::{BRep, FaceHandle};

/// Whether `a` and `b` lie on one plane with opposite oriented normals.
fn opposite_planes(
    a: &popkern_core::geometry::surface::SurfaceOnFace,
    b: &popkern_core::geometry::surface::SurfaceOnFace,
) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry::Plane;
    let (
        Plane {
            origin: oa,
            normal: na,
            ..
        },
        Plane {
            origin: ob,
            normal: nb,
            ..
        },
    ) = (&a.geometry, &b.geometry)
    else {
        return false;
    };
    let sign =
        |s: &popkern_core::geometry::surface::SurfaceOnFace| if s.orientation { 1.0 } else { -1.0 };
    let (na, nb) = (*na * sign(a), *nb * sign(b));
    let (Some(na), Some(nb)) = (na.normalized(), nb.normalized()) else {
        return false;
    };
    na.dot(&nb) < -1.0 + 1e-9
        && (*ob - *oa).dot(&na).abs() < 1e-9 * (1.0 + oa.x.abs().max(oa.y.abs()).max(oa.z.abs()))
}

/// `body` in the single-face sheet form, when it is the two-face form:
/// every shell a pair of faces, exactly one with a surface, sharing all
/// their edges. `None` for any other body (a solid, or already a
/// sheet).
pub fn single_face_sheet(body: &BRep) -> Option<BRep> {
    let mut out = body.clone();
    let pairs: Vec<(FaceHandle, FaceHandle)> = out
        .shells
        .keys()
        .map(|shell| {
            let faces: Vec<FaceHandle> = out
                .faces
                .iter()
                .filter(|(_, f)| f.shell == shell)
                .map(|(h, _)| h)
                .collect();
            let [a, b] = faces[..] else {
                return None;
            };
            match (&out.faces[a].surface, &out.faces[b].surface) {
                (Some(_), None) => Some((a, b)),
                (None, Some(_)) => Some((b, a)),
                // Both carry the plane, facing opposite ways: the front
                // is the face the constructor made first.
                (Some(sa), Some(sb)) if opposite_planes(sa, sb) => Some((a.min(b), a.max(b))),
                _ => None,
            }
        })
        .collect::<Option<_>>()?;
    if pairs.is_empty() {
        return None;
    }
    for (front, back) in pairs {
        let back_face = out.faces.get(back)?;
        if !back_face.inner_loops.is_empty() {
            return None;
        }
        let ring = back_face.outer_loop;
        // Every edge of the backside must be shared with the front.
        let start = out.loops.get(ring)?.half_edge;
        let mut cursor = start;
        let mut moved = Vec::new();
        for _ in 0..=out.half_edges.len() {
            let he = out.half_edges.get(cursor)?;
            let twin_face = out.loops.get(out.half_edges.get(he.twin)?.loop_)?.face;
            if twin_face != front {
                return None;
            }
            moved.push(cursor);
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        // The ring becomes the front's inner loop; the backside goes.
        out.loops.get_mut(ring)?.face = front;
        out.faces.get_mut(front)?.inner_loops.push(ring);
        out.faces.remove(back);
        // Each moved half-edge runs along its twin reversed: its pcurve
        // on the front is the twin's, reversed.
        for h in moved {
            let twin = out.half_edges[h].twin;
            let reversed = out.half_edges[twin].pcurve.as_ref().map(|p| p.reversed());
            out.half_edges[h].pcurve = reversed;
        }
    }
    (out.body_kind() == popkern_core::topology::BodyKind::Sheet).then_some(out)
}

/// Give every boundary half-edge without a twin (a sheet whose faces
/// were retired to an open boundary) a twin of its own: reversed, on
/// the same face, in a loop no face lists — the mirror form the view
/// already presents, so the graph registers the free edge once.
fn twin_open_boundary(view: &mut BRep) {
    use popkern_core::topology::{HalfEdge, Loop};
    let open: Vec<_> = view
        .half_edges
        .iter()
        .filter(|(_, h)| !view.half_edges.contains_key(h.twin))
        .map(|(k, _)| k)
        .collect();
    for h in open {
        let (face, end, edge, pcurve) = {
            let he = &view.half_edges[h];
            let Some(face) = view.loops.get(he.loop_).map(|l| l.face) else {
                continue;
            };
            let Some(end) = view.half_edges.get(he.next).map(|n| n.origin) else {
                continue;
            };
            (face, end, he.edge, he.pcurve.as_ref().map(|p| p.reversed()))
        };
        let ring = view.loops.insert(Loop { half_edge: h, face });
        let twin = view.half_edges.insert(HalfEdge {
            origin: end,
            twin: h,
            next: h,
            prev: h,
            edge,
            loop_: ring,
            pcurve,
        });
        view.half_edges[twin].next = twin;
        view.half_edges[twin].prev = twin;
        view.loops[ring].half_edge = twin;
        view.half_edges[h].twin = twin;
    }
}

#[cfg(test)]
#[path = "lamina_tests.rs"]
mod tests;

/// `sheet` as the boolean pipeline reads it: each face lists no
/// free-boundary mirror among its inner loops (`BRep::is_mirror_loop`).
/// The mirror's half-edges stay, twinned to the outer loop's, but no
/// face lists them, so no hole or containment test reads the mirror
/// as a hole and the graph registers each free edge once — a lamina's
/// free edge in the complex.
pub(crate) fn pipeline_view(sheet: &BRep) -> BRep {
    let mut view = sheet.clone();
    twin_open_boundary(&mut view);
    let faces: Vec<FaceHandle> = view.faces.keys().collect();
    for face in faces {
        let mirrors: Vec<_> = view.faces[face]
            .inner_loops
            .iter()
            .copied()
            .filter(|&l| sheet.is_mirror_loop(face, l))
            .collect();
        view.faces[face]
            .inner_loops
            .retain(|l| !mirrors.contains(l));
    }
    view
}
