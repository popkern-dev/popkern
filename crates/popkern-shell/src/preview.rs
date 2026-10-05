//! Interactive shell/offset preview (PK-328 row 483, design §4).
//!
//! A preview NEVER fits offset surfaces: it samples the exact offset
//! map `S(u,v) + d·m·n̂(u,v)` on a UV grid over each face's trim
//! polygons and emits a light triangle mesh — cheap enough to drive a
//! slider, honest enough to show the true offset skin. Gap edges get
//! approximate pipe strips (normal-slerp arcs around the source edge);
//! corners are left open. The result is a display artifact only:
//! never a `BRep`, never serialized as a model.
//!
//! popkern-shell cannot depend on popkern-mesh (dependency direction),
//! so [`PreviewMesh`] is deliberately minimal; a UI uploads it
//! directly.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

use crate::body_distance::FaceTrimEntry;
use crate::classify::classify_faces;
use crate::error::ShellError;
use crate::types::ShellDirection;

/// UV grid resolution per face axis.
const GRID: usize = 16;
/// Stations along a gap edge strip.
const STRIP_STATIONS: usize = 16;
/// Arc samples across a gap strip.
const STRIP_ARC: usize = 6;

/// A light preview triangle mesh (display-only; see module docs).
#[derive(Debug, Clone, Default)]
pub struct PreviewMesh {
    /// Vertex positions.
    pub positions: Vec<Point3>,
    /// Per-vertex unit normals (offset-map normals; strip arc normals).
    pub normals: Vec<Vec3>,
    /// Triangle indices into `positions`.
    pub triangles: Vec<[usize; 3]>,
    /// Source face per triangle (`None` for join-strip triangles).
    pub face_map: Vec<Option<FaceHandle>>,
}

/// Sample the offset skin of `brep` at `thickness` along `direction`
/// into a preview mesh: exact offset-map samples per face, approximate
/// pipe strips along separating (gap) edges, no trim, no validation.
///
/// `Midplane` previews the two half-distance skins together.
///
/// # Errors
///
/// [`ShellError::NoFaces`] /
/// [`ShellError::MaterialSenseUnresolved`] from classification.
// `tolerance` is carried for API symmetry with the real operations (a
// finer preview grid may consume it later); today only the recursive
// midplane split passes it through.
#[allow(clippy::only_used_in_recursion)]
pub fn shell_preview(
    brep: &BRep,
    thickness: f64,
    direction: ShellDirection,
    tolerance: Tolerance,
) -> Result<PreviewMesh, ShellError> {
    if matches!(direction, ShellDirection::Midplane) {
        let mut mesh = shell_preview(brep, thickness / 2.0, ShellDirection::Inward, tolerance)?;
        let outward = shell_preview(brep, thickness / 2.0, ShellDirection::Outward, tolerance)?;
        let base = mesh.positions.len();
        mesh.positions.extend(outward.positions);
        mesh.normals.extend(outward.normals);
        mesh.triangles.extend(
            outward
                .triangles
                .iter()
                .map(|t| [t[0] + base, t[1] + base, t[2] + base]),
        );
        mesh.face_map.extend(outward.face_map);
        return Ok(mesh);
    }

    let classification = classify_faces(brep, thickness, &[], direction, &[])?;
    let motions: std::collections::HashMap<FaceHandle, f64> = classification
        .offset_faces
        .iter()
        .map(|&(fh, d, inward)| (fh, if inward { -d } else { d }))
        .collect();

    let mut mesh = PreviewMesh::default();
    let mut entries: std::collections::HashMap<FaceHandle, FaceTrimEntry> =
        std::collections::HashMap::new();
    let mut face_handles: Vec<FaceHandle> = motions.keys().copied().collect();
    face_handles.sort();

    for &fh in &face_handles {
        let Some(sof) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(entry) = FaceTrimEntry::build(brep, fh) else {
            continue;
        };
        let Some((u0, v0, u1, v1)) = entry.uv_bbox() else {
            continue;
        };
        let signed = motions[&fh];
        // Grid vertices: Some(index) where in-trim with a healthy
        // normal.
        let mut index: Vec<Option<usize>> = vec![None; (GRID + 1) * (GRID + 1)];
        for i in 0..=GRID {
            for j in 0..=GRID {
                let u = u0 + (u1 - u0) * i as f64 / GRID as f64;
                let v = v0 + (v1 - v0) * j as f64 / GRID as f64;
                if !entry.contains(u, v) {
                    continue;
                }
                let Some(n) = entry.surface_normal((u, v)) else {
                    continue;
                };
                let p = sof.geometry.point_at(u, v) + n * signed;
                index[i * (GRID + 1) + j] = Some(mesh.positions.len());
                mesh.positions.push(p);
                mesh.normals.push(n);
            }
        }
        for i in 0..GRID {
            for j in 0..GRID {
                let q = [
                    index[i * (GRID + 1) + j],
                    index[(i + 1) * (GRID + 1) + j],
                    index[(i + 1) * (GRID + 1) + j + 1],
                    index[i * (GRID + 1) + j + 1],
                ];
                if let [Some(a), Some(b), Some(c), Some(d)] = q {
                    mesh.triangles.push([a, b, c]);
                    mesh.face_map.push(Some(fh));
                    mesh.triangles.push([a, c, d]);
                    mesh.face_map.push(Some(fh));
                }
            }
        }
        entries.insert(fh, entry);
    }

    // Approximate joins: pipe strips along edges whose adjacent offset
    // images SEPARATE (normal-slerp arc around the source edge).
    let mut edge_handles: Vec<_> = brep.edges.keys().collect();
    edge_handles.sort();
    for eh in edge_handles {
        let edge = &brep.edges[eh];
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(l1), Some(l2)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_)) else {
            continue;
        };
        let (fa, fb) = (l1.face, l2.face);
        if fa == fb {
            continue; // seam: one face, no gap
        }
        let (Some(ea), Some(eb)) = (entries.get(&fa), entries.get(&fb)) else {
            continue;
        };
        let (Some(&ma), Some(&mb)) = (motions.get(&fa), motions.get(&fb)) else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        use popkern_core::geometry::curve::CurveEval;
        let (t0, t1) = (curve.t_start, curve.t_end);
        // GAP test once per edge (midpoint): the offsets separate on
        // CONVEX edges moving out of the material and CONCAVE edges
        // moving in. Convexity: B's material-outward normal sees A's
        // interior tangent on its negative side (material angle < π).
        let mid = curve.geometry.point_at((t0 + t1) * 0.5);
        let outward_sign = if matches!(direction, ShellDirection::Outward) {
            1.0
        } else {
            -1.0
        };
        let gap_edge = (|| -> Option<bool> {
            let near_a = ea.nearest(&mid)?;
            let nb = eb.nearest(&mid).and_then(|n| eb.surface_normal(n.uv))?;
            let into_a = interior_tangent(&brep.faces.get(fa)?.surface.as_ref()?.geometry, ea, near_a.uv)?;
            let ob = nb * (mb.signum() * outward_sign);
            let convex = ob.dot(&into_a) < -1e-6;
            if popkern_core::trace::trace_enabled("preview") {
                eprintln!(
                    "[preview-gap] {eh:?} {fa:?}x{fb:?} ob=({:.2},{:.2},{:.2}) into_a=({:.2},{:.2},{:.2}) convex={convex}",
                    ob.x, ob.y, ob.z, into_a.x, into_a.y, into_a.z
                );
            }
            Some(convex == matches!(direction, ShellDirection::Outward))
        })()
        .unwrap_or(false);
        if !gap_edge {
            continue;
        }
        let mut prev_ring: Option<Vec<usize>> = None;
        for s in 0..=STRIP_STATIONS {
            let t = t0 + (t1 - t0) * s as f64 / STRIP_STATIONS as f64;
            let p = curve.geometry.point_at(t);
            let (Some(na), Some(nb)) = (
                ea.nearest(&p).and_then(|n| ea.surface_normal(n.uv)),
                eb.nearest(&p).and_then(|n| eb.surface_normal(n.uv)),
            ) else {
                prev_ring = None;
                continue;
            };
            let (da, db) = (na * ma, nb * mb);
            let mag = ma.abs().max(mb.abs()).max(1e-12);
            // Slerp the two displaced directions around the station.
            let (ua, ub) = (da * (1.0 / da.length()), db * (1.0 / db.length()));
            let dot = ua.dot(&ub).clamp(-1.0, 1.0);
            let angle = dot.pacos();
            let ring: Vec<usize> = (0..=STRIP_ARC)
                .map(|k| {
                    let f = k as f64 / STRIP_ARC as f64;
                    let dir = if angle < 1e-6 {
                        ua
                    } else {
                        let (wa, wb) = (
                            ((1.0 - f) * angle).psin() / angle.psin(),
                            (f * angle).psin() / angle.psin(),
                        );
                        let v = ua * wa + ub * wb;
                        v * (1.0 / v.length().max(1e-12))
                    };
                    let idx = mesh.positions.len();
                    mesh.positions.push(p + dir * mag);
                    mesh.normals.push(dir);
                    idx
                })
                .collect();
            if let Some(prev) = prev_ring {
                for k in 0..STRIP_ARC {
                    mesh.triangles.push([prev[k], ring[k], ring[k + 1]]);
                    mesh.face_map.push(None);
                    mesh.triangles.push([prev[k], ring[k + 1], prev[k + 1]]);
                    mesh.face_map.push(None);
                }
            }
            prev_ring = Some(ring);
        }
    }
    Ok(mesh)
}

/// A unit tangent from `uv` toward the face interior: step the UV in
/// the four axis directions and take the first in-trim step's 3D
/// image direction.
fn interior_tangent(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    entry: &FaceTrimEntry,
    uv: (f64, f64),
) -> Option<Vec3> {
    let (u0, v0, u1, v1) = entry.uv_bbox()?;
    let (du, dv) = (
        ((u1 - u0).abs()).max(1e-9) * 0.02,
        ((v1 - v0).abs()).max(1e-9) * 0.02,
    );
    let p = surface.point_at(uv.0, uv.1);
    let dir_of = |su: f64, sv: f64| -> Option<Vec3> {
        let q = surface.point_at(uv.0 + su, uv.1 + sv);
        let d = q - p;
        let len = d.length();
        (len > 1e-12).then(|| d * (1.0 / len))
    };
    // A genuine INTO-face direction crosses the trim boundary: the
    // step is in-trim while its opposite is out (a step parallel to
    // the boundary keeps both sides in-trim and is NOT into the face
    // — measured: a cube edge read its own edge direction and the
    // convexity dot came back 0).
    for (su, sv) in [(du, 0.0), (-du, 0.0), (0.0, dv), (0.0, -dv)] {
        if entry.contains(uv.0 + su, uv.1 + sv)
            && !entry.contains(uv.0 - su, uv.1 - sv)
            && let Some(d) = dir_of(su, sv)
        {
            return Some(d);
        }
    }
    for (su, sv) in [(du, 0.0), (-du, 0.0), (0.0, dv), (0.0, -dv)] {
        if entry.contains(uv.0 + su, uv.1 + sv)
            && let Some(d) = dir_of(su, sv)
        {
            return Some(d);
        }
    }
    None
}

#[cfg(test)]
#[path = "preview_tests.rs"]
mod tests;
