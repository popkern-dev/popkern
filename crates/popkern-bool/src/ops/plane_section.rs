//! Exact plane sections (PK-856): the regions where a plane cuts a
//! solid, as closed loops of the body's own exact curves.
//!
//! The section is the boolean's own answer: the body intersected with
//! the plane's half-space (`half_space`, a prism whose cap is the
//! plane), whose faces ON the plane are the section regions — trimmed
//! by exact intersection curves, holes and all, validated like any
//! boolean result. Each region's area and centroid integrate Green's
//! theorem along its loops' exact curves (Gauss-Legendre on each
//! half-edge), so they are exact to quadrature, not to a polyline; the
//! loops are also sampled within a chord height for drawing and export.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use crate::error::BoolError;

/// A section plane: through `origin`, with unit `normal`, and the 2-D
/// frame `(u_axis, normal x u_axis)` its loops are written in.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SectionPlane {
    /// A point on the plane (the 2-D frame's origin).
    pub origin: Point3,
    /// The plane's normal; the section regions face along it.
    pub normal: Vec3,
    /// The 2-D frame's x direction (in the plane).
    pub u_axis: Vec3,
}

impl SectionPlane {
    /// The plane through `origin` with `normal`, its 2-D x axis the
    /// in-plane direction least parallel to a world axis
    /// (`orthonormal_basis`). `None` for a zero normal.
    pub fn new(origin: Point3, normal: Vec3) -> Option<SectionPlane> {
        let normal = normal.normalized()?;
        let (u_axis, _) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(normal);
        Some(SectionPlane {
            origin,
            normal,
            u_axis,
        })
    }

    /// The 2-D frame's y direction.
    pub fn v_axis(&self) -> Vec3 {
        self.normal.cross(&self.u_axis)
    }

    /// `p` in the plane's 2-D frame.
    pub fn to_2d(&self, p: &Point3) -> [f64; 2] {
        let r = *p - self.origin;
        [r.dot(&self.u_axis), r.dot(&self.v_axis())]
    }
}

/// One closed loop of a section region, in the plane's 2-D frame.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionLoop {
    /// Points along the loop within the chord height, first not
    /// repeated: counter-clockwise for an outer boundary, clockwise for
    /// a hole (seen from the normal's tip).
    pub points: Vec<[f64; 2]>,
}

/// One connected region where the plane cuts the body.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionRegion {
    /// The outer boundary.
    pub outer: SectionLoop,
    /// The holes.
    pub holes: Vec<SectionLoop>,
    /// The region's area (exact curves, holes subtracted).
    pub area: f64,
    /// The region's centroid in the 2-D frame.
    pub centroid: [f64; 2],
}

/// A body's plane section.
#[derive(Debug, Clone, PartialEq)]
pub struct PlaneSection {
    /// The plane and its 2-D frame.
    pub plane: SectionPlane,
    /// The regions, largest first.
    pub regions: Vec<SectionRegion>,
}

impl PlaneSection {
    /// Total area of the regions.
    pub fn area(&self) -> f64 {
        self.regions.iter().map(|r| r.area).sum()
    }
}

/// The section of solid `body` by `plane`, its loops sampled within
/// `chord` (a length). Empty when the plane misses the body.
///
/// # Errors
///
/// [`BoolError`] when the boolean declines (the half-space cut of a
/// body it cannot serve), or a region face has no planar surface.
pub fn plane_section(
    body: &BRep,
    plane: &SectionPlane,
    chord: f64,
    tolerance: &Tolerance,
) -> Result<PlaneSection, BoolError> {
    let empty = PlaneSection {
        plane: *plane,
        regions: Vec::new(),
    };
    // The body below the plane; the half-space covers all of it when
    // the body lies below, and none when it lies above.
    let Some(cutter) = super::half_space::half_space(body, plane.origin, plane.normal, *tolerance)?
    else {
        return Ok(empty);
    };
    let below = crate::intersection_ref(body, &cutter, *tolerance, None)?;
    let on_plane = |f: FaceHandle| -> bool {
        let Some(SurfaceGeometry::Plane { origin, normal, .. }) = below
            .faces
            .get(f)
            .and_then(|x| x.surface.as_ref())
            .map(|s| &s.geometry)
        else {
            return false;
        };
        let parallel = normal.cross(&plane.normal).length() <= tolerance.angular();
        parallel && (*origin - plane.origin).dot(&plane.normal).abs() <= 10.0 * tolerance.linear()
    };
    let chord = if chord.is_finite() && chord > 0.0 {
        chord
    } else {
        1e-3
    };
    let mut regions = Vec::new();
    for (fh, face) in &below.faces {
        if !on_plane(fh) {
            continue;
        }
        let outer = section_loop(&below, face.outer_loop, plane, chord)?;
        let holes = face
            .inner_loops
            .iter()
            .map(|&l| section_loop(&below, l, plane, chord))
            .collect::<Result<Vec<_>, _>>()?;
        // Green's theorem over every loop, each walked as the face
        // walks it: the cap faces along the plane's normal, so the
        // holes subtract by their orientation.
        let (mut a2, mut mx, mut my) = (0.0, 0.0, 0.0);
        for lh in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
            let (da, dx, dy) = loop_moments(&below, lh, plane)?;
            a2 += da;
            mx += dx;
            my += dy;
        }
        let area = 0.5 * a2;
        let centroid = if area.abs() > 0.0 {
            [mx / (2.0 * area), my / (2.0 * area)]
        } else {
            [0.0, 0.0]
        };
        regions.push(SectionRegion {
            outer: outer.0,
            holes: holes.into_iter().map(|h| h.0).collect(),
            area: area.abs(),
            centroid,
        });
    }
    regions.sort_by(|a, b| b.area.total_cmp(&a.area));
    Ok(PlaneSection {
        plane: *plane,
        regions,
    })
}

/// The half-edges of loop `lh`, in order.
fn loop_half_edges(brep: &BRep, lh: LoopHandle) -> Result<Vec<HalfEdgeHandle>, BoolError> {
    let bad = || BoolError::DegenerateInput("section loop does not close".into());
    let start = brep.loops.get(lh).ok_or_else(bad)?.half_edge;
    let mut out = Vec::new();
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(he);
        he = brep.half_edges.get(he).ok_or_else(bad)?.next;
        if he == start {
            return Ok(out);
        }
    }
    Err(bad())
}

/// Half-edge `he`'s curve and the parameter interval it walks, in its
/// own direction (`pcurve_recompute`'s convention: an open edge starts
/// at its origin; a closed edge's own half-edge runs forward, its twin
/// backward). `None` for a curveless edge.
fn walked_curve(
    brep: &BRep,
    he: HalfEdgeHandle,
) -> Option<(&popkern_core::geometry::curve::CurveGeometry, f64, f64)> {
    let h = brep.half_edges.get(he)?;
    let edge = brep.edges.get(h.edge)?;
    let curve = edge.curve.as_ref()?;
    let (t0, t1) = (curve.t_start, curve.t_end);
    let (p0, p1) = (curve.geometry.point_at(t0), curve.geometry.point_at(t1));
    let origin = brep.vertices.get(h.origin)?.position;
    let forward = if p0.distance_to(&p1) < 1e-9 {
        edge.half_edge == he
    } else {
        origin.distance_to(&p0) <= origin.distance_to(&p1)
    };
    Some(if forward {
        (&curve.geometry, t0, t1)
    } else {
        (&curve.geometry, t1, t0)
    })
}

/// Loop `lh` sampled within `chord`.
fn section_loop(
    brep: &BRep,
    lh: LoopHandle,
    plane: &SectionPlane,
    chord: f64,
) -> Result<(SectionLoop,), BoolError> {
    let mut points = Vec::new();
    for he in loop_half_edges(brep, lh)? {
        let Some((curve, ta, tb)) = walked_curve(brep, he) else {
            // A curveless edge is a straight chord between its ends.
            let h = &brep.half_edges[he];
            points.push(plane.to_2d(&brep.vertices[h.origin].position));
            continue;
        };
        let mut ts = vec![ta, tb];
        // Bisect every span whose midpoint sags more than the chord
        // height, breadth first under a cap.
        for _ in 0..16 {
            let mut next = vec![ts[0]];
            let mut split = false;
            for w in ts.windows(2) {
                let (a, b) = (w[0], w[1]);
                let m = 0.5 * (a + b);
                let (pa, pb, pm) = (curve.point_at(a), curve.point_at(b), curve.point_at(m));
                let mid = Point3::new(
                    0.5 * (pa.x + pb.x),
                    0.5 * (pa.y + pb.y),
                    0.5 * (pa.z + pb.z),
                );
                // A closed curve's first span is the whole loop: always
                // split it.
                if pm.distance_to(&mid) > chord || pa.distance_to(&pb) < 1e-12 {
                    next.push(m);
                    split = true;
                }
                next.push(b);
            }
            ts = next;
            if !split || ts.len() > 1 << 16 {
                break;
            }
        }
        // Every sample but the last: the next half-edge starts there.
        for &t in &ts[..ts.len() - 1] {
            points.push(plane.to_2d(&curve.point_at(t)));
        }
    }
    Ok((SectionLoop { points },))
}

/// Twice the signed area, and the first moments `(∮ x^2 dy, -∮ y^2 dx)`
/// of loop `lh`, by Gauss-Legendre on each half-edge's exact curve.
fn loop_moments(
    brep: &BRep,
    lh: LoopHandle,
    plane: &SectionPlane,
) -> Result<(f64, f64, f64), BoolError> {
    // 8-point Gauss-Legendre nodes and weights on [-1, 1].
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let (u, v) = (plane.u_axis, plane.v_axis());
    let (mut a2, mut mx, mut my) = (0.0, 0.0, 0.0);
    for he in loop_half_edges(brep, lh)? {
        let Some((curve, ta, tb)) = walked_curve(brep, he) else {
            continue;
        };
        // Sub-intervals: exact for lines, converged for conics and the
        // smooth NURBS of intersection curves.
        let pieces = 32;
        for k in 0..pieces {
            let a = ta + (tb - ta) * k as f64 / pieces as f64;
            let b = ta + (tb - ta) * (k + 1) as f64 / pieces as f64;
            let (c, h) = (0.5 * (a + b), 0.5 * (b - a));
            for i in 0..4 {
                for s in [-1.0, 1.0] {
                    let t = c + s * h * X[i];
                    let p = plane.to_2d(&curve.point_at(t));
                    let d = curve.derivative_at(t, 1);
                    let (dx, dy) = (d.dot(&u) * h, d.dot(&v) * h);
                    let w = W[i];
                    a2 += w * (p[0] * dy - p[1] * dx);
                    mx += w * p[0] * p[0] * dy;
                    my -= w * p[1] * p[1] * dx;
                }
            }
        }
    }
    Ok((a2, mx, my))
}

#[cfg(test)]
#[path = "plane_section_tests.rs"]
mod tests;
