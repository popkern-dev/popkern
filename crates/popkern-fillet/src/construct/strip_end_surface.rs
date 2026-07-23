//! Strip-end SURFACE support: the deviation-test abstraction over a
//! blend strip's surface (cylindrical fillet strips and planar
//! chamfer bevels share the end-corner defect shape), and the
//! missing-pcurve backfill that follows end-corner contraction.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::strip_end_arcs::radial_dist;

/// The strip's surface for end-corner deviation tests: how far a
/// point sits off the blend surface. Cylinder strips are fillets;
/// PLANE strips are flat chamfer bevels — their mitred end corners
/// are the same defect shape (corner retained off-surface between
/// two on-surface feet) with a planar distance function.
pub(super) enum StripSurface {
    Cylinder {
        origin: Point3,
        axis: Vec3,
        radius: f64,
    },
    Plane {
        origin: Point3,
        normal: Vec3,
    },
}

impl StripSurface {
    pub(super) fn of(brep: &BRep, face: FaceHandle) -> Option<StripSurface> {
        match brep.faces.get(face)?.surface.as_ref()?.geometry {
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => Some(StripSurface::Cylinder {
                origin,
                axis,
                radius,
            }),
            SurfaceGeometry::Plane { origin, normal, .. } => {
                Some(StripSurface::Plane { origin, normal })
            }
            _ => None,
        }
    }

    /// Unsigned deviation of `p` from the surface.
    pub(super) fn deviation(&self, p: &Point3) -> f64 {
        match self {
            StripSurface::Cylinder {
                origin,
                axis,
                radius,
            } => (radial_dist(origin, axis, p) - radius).abs(),
            StripSurface::Plane { origin, normal } => {
                let n = normal.length();
                if n < 1e-300 {
                    f64::INFINITY
                } else {
                    ((*p - *origin).dot(normal) / n).abs()
                }
            }
        }
    }

    /// Signed outward excess (used by the corner test: the corner
    /// must sit OFF the surface, beyond the gate).
    pub(super) fn excess(&self, p: &Point3) -> f64 {
        match self {
            StripSurface::Cylinder {
                origin,
                axis,
                radius,
            } => radial_dist(origin, axis, p) - radius,
            StripSurface::Plane { .. } => self.deviation(p),
        }
    }
}

/// Attach exact straight pcurves to half-edges that have NONE, on a
/// planar face: line edges on a plane project to line segments in
/// the chart. Non-planar faces and half-edges that already carry a
/// pcurve are untouched.
pub(super) fn fill_missing_planar_pcurves(brep: &mut BRep, face: FaceHandle) {
    use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let Some(SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    }) = face_data.surface.as_ref().map(|s| s.geometry.clone())
    else {
        return;
    };
    let mut loops = vec![face_data.outer_loop];
    loops.extend(face_data.inner_loops.iter().copied());
    for lh in loops {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(d) = brep.half_edges.get(he) else {
                break;
            };
            let next = d.next;
            if d.pcurve.is_none() {
                let p0 = brep.vertices[d.origin].position;
                let p1 = brep.vertices[brep.half_edges[next].origin].position;
                let uv = |p: popkern_core::geometry::Point3| {
                    let dvec = p - origin;
                    (dvec.dot(&u_axis), dvec.dot(&v_axis))
                };
                let pc = Pcurve {
                    geometry: PcurveGeometry::Linear {
                        start: uv(p0),
                        end: uv(p1),
                    },
                    t_start: 0.0,
                    t_end: 1.0,
                };
                if let Some(hm) = brep.half_edges.get_mut(he) {
                    hm.pcurve = Some(pc);
                }
            }
            he = next;
            if he == start {
                break;
            }
        }
    }
}

/// Phase 3: remove zero-width spur slits (a half-edge immediately
/// followed by its own twin, dangling vertex valence 1) from the end
/// face's loops. Earlier rebuilds left the amputated L-legs doubled
/// in the base loop; the slit has zero area, so removal preserves
/// the face point set by construction.
pub(crate) fn remove_spur_slits(brep: &mut BRep, plane_face: FaceHandle) -> usize {
    let mut removed = 0;
    for _ in 0..16 {
        let Some(face) = brep.faces.get(plane_face) else {
            return removed;
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        let mut target: Option<HalfEdgeHandle> = None;
        'outer: for lh in loops {
            let Some(l) = brep.loops.get(lh) else {
                continue;
            };
            let start_he = l.half_edge;
            let mut he = start_he;
            for _ in 0..brep.half_edges.len() + 1 {
                let d = &brep.half_edges[he];
                if d.next == d.twin {
                    target = Some(he);
                    break 'outer;
                }
                he = d.next;
                if he == start_he {
                    break;
                }
            }
        }
        let Some(he) = target else {
            return removed;
        };
        let mut candidate = brep.clone();
        match candidate.kev(he) {
            Ok(_) => {
                if candidate.validate().is_ok() {
                    *brep = candidate;
                    removed += 1;
                    if std::env::var_os("PK_CARVE_TRACE").is_some() {
                        eprintln!("  [strip-end] {plane_face:?}: spur slit removed");
                    }
                } else {
                    if std::env::var_os("PK_CARVE_TRACE").is_some() {
                        eprintln!("  [strip-end] {plane_face:?}: spur removal ROLLED BACK");
                    }
                    return removed;
                }
            }
            Err(e) => {
                if std::env::var_os("PK_CARVE_TRACE").is_some() {
                    eprintln!("  [strip-end] {plane_face:?}: kev declined: {e:?}");
                }
                return removed;
            }
        }
    }
    removed
}
