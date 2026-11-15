//! The shell of a convex body as the body less its cavity (PK-1212).
//!
//! Shelling moves each kept face inward by the thickness and runs the
//! moved faces on, extended, until they meet one another and the
//! surfaces of the openings. Where that cavity never reaches an opening
//! (a ball cut below its offset sphere's top, a frustum whose offset
//! cone closes to its apex under the top) the wall's rim solve has no
//! rim to find, and the skin assembly declined with open loops; where
//! the cavity is empty (a frustum too thick to hollow, a box whose walls
//! meet) or a consumed round's neighbours meet past it, the offset
//! declined as exceeding its maximum.
//!
//! For a convex body the cavity is exactly the intersection of the inner
//! sides of the kept faces' offset surfaces: each face's supporting
//! region (a half-space, a ball, a solid cylinder, a cone's nappe) holds
//! the whole body, so the body is their intersection, and moving each
//! inward by the thickness and intersecting gives the points of the body
//! farther than the thickness from every kept face's extended surface. A
//! curved face whose offset collapses (radius at or under the tolerance)
//! is consumed and drops out: its neighbours' extended offsets meet past
//! it. An opening adds no constraint: beyond it the cavity lies outside
//! the body. The shell is then the body less the cavity, one boolean:
//! open where the cavity crosses an opening, a void where it does not,
//! the body itself where the cavity is empty.
//!
//! Applicability is measured, not assumed: every face's region must hold
//! the body's vertices, samples along every edge, and samples inside
//! every face's trim (a bump on a face holds its rim on the plane and
//! its crown beyond it). Only uniform inward shells take this path, and
//! only after the skin path declined; the caller records it as
//! [`crate::ShellWarning::ConvexCavity`].

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, PortableMath, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

/// A face's supporting region, the side holding the body.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Region {
    /// `(p − point) · out ≤ 0`.
    Half { point: Point3, out: Vec3 },
    /// Within `radius` of `centre`.
    Ball { centre: Point3, radius: f64 },
    /// Within `radius` of the line through `origin` along `axis`.
    Rod {
        origin: Point3,
        axis: Vec3,
        radius: f64,
    },
    /// The nappe from `apex` along `dir` with half-angle `angle`.
    Nappe { apex: Point3, dir: Vec3, angle: f64 },
}

impl Region {
    /// How far `p` lies outside the region (negative inside), in length.
    fn excess(&self, p: &Point3) -> f64 {
        match *self {
            Region::Half { point, out } => (*p - point).dot(&out),
            Region::Ball { centre, radius } => p.distance_to(&centre) - radius,
            Region::Rod {
                origin,
                axis,
                radius,
            } => {
                let d = *p - origin;
                (d - axis * d.dot(&axis)).length() - radius
            }
            Region::Nappe { apex, dir, angle } => {
                // Signed distance to the nappe's surface, behind the apex
                // the distance to the apex itself.
                let d = *p - apex;
                let h = d.dot(&dir);
                let r = (d - dir * h).length();
                let (s, c) = (angle.psin(), angle.pcos());
                let along = h * c + r * s;
                if along < 0.0 {
                    d.length()
                } else {
                    r * c - h * s
                }
            }
        }
    }

    /// The region moved inward by `t`; `None` when it collapses (a
    /// consumed face).
    fn offset(&self, t: f64, lin: f64) -> Option<Region> {
        match *self {
            Region::Half { point, out } => Some(Region::Half {
                point: point - out * t,
                out,
            }),
            Region::Ball { centre, radius } => (radius - t > lin).then_some(Region::Ball {
                centre,
                radius: radius - t,
            }),
            Region::Rod {
                origin,
                axis,
                radius,
            } => (radius - t > lin).then_some(Region::Rod {
                origin,
                axis,
                radius: radius - t,
            }),
            Region::Nappe { apex, dir, angle } => Some(Region::Nappe {
                apex: apex + dir * (t / angle.psin()),
                dir,
                angle,
            }),
        }
    }

    /// The region as a solid reaching `reach` past `centre` every way.
    pub(crate) fn solid(&self, centre: Point3, reach: f64, tol: Tolerance) -> Option<BRep> {
        let placed = |origin: Point3, z: Vec3| -> Option<Transform3> {
            let seed = if z.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let x = (seed - z * seed.dot(&z)).normalized()?;
            Transform3::from_axis_placement(origin, z, x).ok()
        };
        match *self {
            Region::Half { point, out } => {
                // A slab from the plane `reach` deep into the inner side,
                // `reach` wide about the body's centre.
                let foot = centre - out * (centre - point).dot(&out);
                let l = reach;
                let square = [
                    Point3::new(-l, -l, 0.0),
                    Point3::new(l, -l, 0.0),
                    Point3::new(l, l, 0.0),
                    Point3::new(-l, l, 0.0),
                ];
                BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 2.0 * l, tol)
                    .ok()?
                    .transform(&placed(foot, -out)?)
                    .ok()
            }
            Region::Ball { centre: c, radius } => BRep::sphere(radius, tol)
                .ok()?
                .transform(&Transform3::translation(c.x, c.y, c.z))
                .ok(),
            Region::Rod {
                origin,
                axis,
                radius,
            } => {
                let mid = origin + axis * (centre - origin).dot(&axis);
                BRep::cylinder(radius, 2.0 * reach, tol)
                    .ok()?
                    .transform(&placed(mid - axis * reach, axis)?)
                    .ok()
            }
            Region::Nappe { apex, dir, angle } => {
                // The core cone stands on its base at z = 0 with its apex
                // at z = height: its z runs from the base to the apex,
                // against the nappe's opening.
                let height = apex.distance_to(&centre) + 2.0 * reach;
                let radius = height * angle.psin() / angle.pcos();
                BRep::cone(radius, height, tol)
                    .ok()?
                    .transform(&placed(apex + dir * height, -dir)?)
                    .ok()
            }
        }
    }
}

/// The face's supporting region holding `samples`, within `reach`.
pub(crate) fn region_of(
    surface: &SurfaceGeometry,
    samples: &[Point3],
    reach: f64,
) -> Option<Region> {
    let holds = |r: &Region| samples.iter().all(|p| r.excess(p) <= reach);
    let region = match *surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let n = normal.normalized()?;
            let up = Region::Half {
                point: origin,
                out: n,
            };
            let down = Region::Half {
                point: origin,
                out: -n,
            };
            return [up, down].into_iter().find(holds);
        }
        SurfaceGeometry::Sphere { center, radius, .. } => Region::Ball {
            centre: center,
            radius,
        },
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Region::Rod {
            origin,
            axis: axis.normalized()?,
            radius,
        },
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let a = axis.normalized()?;
            let angle = half_angle.abs();
            if angle <= 0.0 || angle >= std::f64::consts::FRAC_PI_2 {
                return None;
            }
            let up = Region::Nappe {
                apex,
                dir: a,
                angle,
            };
            let down = Region::Nappe {
                apex,
                dir: -a,
                angle,
            };
            return [up, down].into_iter().find(holds);
        }
        _ => return None,
    };
    holds(&region).then_some(region)
}

/// Points of `body`: its vertices, samples along each edge, and samples
/// inside each face's trim.
pub(crate) fn body_samples(body: &BRep) -> Vec<Point3> {
    let mut out: Vec<Point3> = body.vertices.values().map(|v| v.position).collect();
    for edge in body.edges.values() {
        if let Some(c) = edge.curve.as_ref() {
            for k in 0..=EDGE_SAMPLES {
                let s = k as f64 / EDGE_SAMPLES as f64;
                out.push(c.geometry.point_at(c.t_start + s * (c.t_end - c.t_start)));
            }
        }
    }
    for (face, data) in &body.faces {
        let Some(sof) = data.surface.as_ref() else {
            continue;
        };
        out.extend(face_samples(body, face, &sof.geometry));
    }
    out
}

/// Samples along each edge.
const EDGE_SAMPLES: usize = 16;
/// Grid lines across a face's parameter box.
const FACE_GRID: usize = 9;

/// The box of `face`'s outer loop in its chart, from the loop's own
/// pcurves (one consistent branch on a periodic chart; the vertex-based
/// bounded range picked a quarter-turn round's window off its period).
fn loop_uv_box(body: &BRep, face: FaceHandle) -> Option<((f64, f64), (f64, f64))> {
    let start = body.loops.get(body.faces.get(face)?.outer_loop)?.half_edge;
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    let mut h = start;
    for _ in 0..body.half_edges.len() {
        let he = body.half_edges.get(h)?;
        let pc = he.pcurve.as_ref()?;
        for k in 0..=EDGE_SAMPLES {
            let q =
                pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * k as f64 / EDGE_SAMPLES as f64);
            lo = (lo.0.min(q.x), lo.1.min(q.y));
            hi = (hi.0.max(q.x), hi.1.max(q.y));
        }
        h = he.next;
        if h == start {
            return (lo.0 < hi.0 && lo.1 < hi.1).then_some(((lo.0, hi.0), (lo.1, hi.1)));
        }
    }
    None
}

/// Points of `surface` inside `face`'s trim on a grid over its box.
fn face_samples(body: &BRep, face: FaceHandle, surface: &SurfaceGeometry) -> Vec<Point3> {
    let Some(((u0, u1), (v0, v1))) = loop_uv_box(body, face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for i in 1..FACE_GRID {
        for j in 1..FACE_GRID {
            let u = u0 + (u1 - u0) * i as f64 / FACE_GRID as f64;
            let v = v0 + (v1 - v0) * j as f64 / FACE_GRID as f64;
            if body.face_uv_trim_containment(face, surface, u, v)
                == Some(FaceUvContainment::Interior)
            {
                out.push(surface.point_at(u, v));
            }
        }
    }
    out
}

/// Whether `face` is a convex round or corner ball whose offset by
/// `consumed_at` collapses: its radius at or under it, and the tangent
/// plane at every sample of it a supporting plane of the body with the
/// axis or centre on the body's side. Its tangent half-spaces moved
/// inward are then implied by the two extreme ones, its neighbours'
/// (a rounded plate shelled thicker than its corner radius meets in
/// sharp inner corners), so it adds no constraint to the cavity. A
/// round whose offset survives is not one region and declines.
pub(crate) fn consumed_convex_blend(
    body: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    samples: &[Point3],
    consumed_at: f64,
    band: f64,
) -> bool {
    let outward = |q: &Point3| -> Option<Vec3> {
        match *surface {
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } if radius <= consumed_at => {
                let a = axis.normalized()?;
                let d = *q - origin;
                (d - a * d.dot(&a)).normalized()
            }
            SurfaceGeometry::Sphere { center, radius, .. } if radius <= consumed_at => {
                (*q - center).normalized()
            }
            _ => None,
        }
    };
    let points = face_samples(body, face, surface);
    !points.is_empty()
        && points
            .iter()
            .all(|q| outward(q).is_some_and(|n| samples.iter().all(|p| (*p - *q).dot(&n) <= band)))
}

/// The shell of `body` with `removed` open, `thickness` inward, as the
/// body less its cavity; `None` when the body is not convex in this
/// sense or a boolean declines (the caller keeps its own error).
pub(crate) fn convex_cavity_shell(
    body: &BRep,
    thickness: f64,
    removed: &[FaceHandle],
    tolerance: Tolerance,
) -> Option<BRep> {
    use popkern_bool::error::BoolOperation;
    let lin = tolerance.linear();
    if thickness.is_nan() || thickness <= lin {
        return None;
    }
    let samples = body_samples(body);
    let (lo, hi) = samples.iter().fold(
        (
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            Point3::new(f64::MIN, f64::MIN, f64::MIN),
        ),
        |(lo, hi), p| {
            (
                Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
                Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
            )
        },
    );
    let size = lo.distance_to(&hi);
    if !(size.is_finite() && size > lin) {
        return None;
    }
    let centre = Point3::new(
        0.5 * (lo.x + hi.x),
        0.5 * (lo.y + hi.y),
        0.5 * (lo.z + hi.z),
    );
    let reach = 2.0 * (size + thickness);
    // Containment within ten tolerances: samples on a face's own surface
    // and on edges where two regions meet read zero up to rounding.
    let band = 10.0 * lin;
    let mut cavity: Option<BRep> = None;
    for (face, data) in &body.faces {
        let surface = &data.surface.as_ref()?.geometry;
        let Some(region) = region_of(surface, &samples, band) else {
            if consumed_convex_blend(body, face, surface, &samples, thickness + lin, band) {
                continue;
            }
            return None;
        };
        if removed.contains(&face) {
            continue;
        }
        let Some(moved) = region.offset(thickness, lin) else {
            continue;
        };
        let solid = moved.solid(centre, reach, tolerance)?;
        cavity = Some(match cavity {
            None => solid,
            Some(c) => {
                match popkern_bool::boolean(BoolOperation::Intersection, c, solid, tolerance, None)
                {
                    Ok(b) if b.faces.is_empty() => return empty_cavity(body, removed),
                    Ok(b) => b,
                    // Convex regions meet in a convex set: touching
                    // without overlap is no volume, an empty cavity.
                    Err(popkern_bool::error::BoolError::NonManifoldResult { .. }) => {
                        return empty_cavity(body, removed);
                    }
                    Err(_) => return None,
                }
            }
        });
    }
    // Every face removed or consumed: nothing bounds the cavity.
    let cavity = cavity?;
    let shelled = match popkern_bool::boolean(
        BoolOperation::Subtraction,
        body.clone(),
        cavity,
        tolerance,
        None,
    ) {
        Ok(b) => b,
        Err(_) => return None,
    };
    shelled
        .validate_collect_structural()
        .is_empty()
        .then_some(shelled)
}

/// The shell whose cavity is empty: with an opening, the part left
/// whole (the walls meet before any hollow forms, as the frustum too
/// thick to hollow is specified); with none, no void exists to make and
/// the skin path's decline stands (`OffsetTooLarge` with its bound, the
/// slider contract of `max_safe_offset`).
fn empty_cavity(body: &BRep, removed: &[FaceHandle]) -> Option<BRep> {
    (!removed.is_empty()).then(|| body.clone())
}

/// After the skin path declined `error`: a uniform inward shell of a
/// convex body served as the body less its cavity, recorded as
/// [`crate::ShellWarning::ConvexCavity`]; otherwise `error` stands.
pub(crate) fn after_decline(
    body: &BRep,
    thickness: f64,
    direction: crate::ShellDirection,
    removed: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    tolerance: Tolerance,
    error: crate::ShellError,
) -> Result<crate::ShellResult, crate::ShellError> {
    if direction != crate::ShellDirection::Inward || !face_thicknesses.is_empty() {
        return Err(error);
    }
    match convex_cavity_shell(body, thickness, removed, tolerance) {
        Some(brep) => Ok(crate::ShellResult {
            brep,
            warnings: vec![crate::ShellWarning::ConvexCavity {
                declined: error.to_string(),
            }],
        }),
        None => Err(error),
    }
}
