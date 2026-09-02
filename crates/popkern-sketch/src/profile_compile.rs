//! Compilation of an [`AnalyticProfile<f64>`] to exact geometry.
//!
//! Two consumers, two outputs, one description:
//!
//! - [`AnalyticProfile::local_wire`] / [`AnalyticProfile::placed_wire`]:
//!   one closed rational NURBS wire (lines exact with unit weights,
//!   arcs exact rational quadratics, the circle in its standard
//!   9-control-point form) for the sweep kernel, whose canonical
//!   recognition maps the arc spans to true `Cylinder` walls.
//! - [`AnalyticProfile::to_brep`]: a planar sheet face whose edges
//!   carry exact `Line` / `Arc` curves, for the history replay path
//!   (the extrude operator promotes an `Arc` edge extruded along its
//!   normal to a `Cylinder`, PK-377).
//!
//! The face path reuses the constraint dialect's tracer and builder
//! with the outline loaded into a [`Sketch`] as a plain geometry
//! container: zero constraints, never solved, never persisted — the
//! document keeps the analytic description (design D2 forbids
//! synthesising a constraint sketch), and the solver is never on the
//! path (design invariant 6).
//!
//! **Tolerance:** the face path takes the caller's explicit
//! [`Tolerance`] for endpoint merging; wire compilation is tolerance
//! free (exact arithmetic on the description).

use popkern_constraint::Sketch;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;

use crate::error::SketchError;
use crate::exact_arc::{
    QuadraticSpan, arc_spans, assemble_quadratic, circle_wire, line_span, polygon_wire,
};
use crate::plane::SketchPlane;
use crate::profile::{AnalyticProfile, PathSegment, ProfileError, ProfileShape};
use crate::{SketchBrepResult, SketchOnPlane};

/// A circular arc in sketch coordinates: centre, radius, start angle,
/// and SIGNED sweep (positive = counter-clockwise).
#[derive(Debug, Clone, Copy)]
struct Arc2 {
    cx: f64,
    cy: f64,
    r: f64,
    a0: f64,
    sweep: f64,
}

/// The circle through `a`, `via`, `b`, traversed from `a` to `b` on
/// the side that passes through `via`.
fn three_point_arc(a: (f64, f64), via: (f64, f64), b: (f64, f64)) -> Result<Arc2, ProfileError> {
    let (ax, ay) = a;
    let (vx, vy) = via;
    let (bx, by) = b;
    let d = 2.0 * (ax * (vy - by) + vx * (by - ay) + bx * (ay - vy));
    if d.abs() < 1.0e-12 {
        return Err(ProfileError::CollinearArc);
    }
    let aa = ax * ax + ay * ay;
    let vv = vx * vx + vy * vy;
    let bb = bx * bx + by * by;
    let cx = (aa * (vy - by) + vv * (by - ay) + bb * (ay - vy)) / d;
    let cy = (aa * (bx - vx) + vv * (ax - bx) + bb * (vx - ax)) / d;
    let r = ((ax - cx).powi(2) + (ay - cy).powi(2)).sqrt();
    let a0 = (ay - cy).atan2(ax - cx);
    let a_v = (vy - cy).atan2(vx - cx);
    let a1 = (by - cy).atan2(bx - cx);
    let norm = |t: f64| t.rem_euclid(std::f64::consts::TAU);
    let ccw_v = norm(a_v - a0);
    let ccw_1 = norm(a1 - a0);
    let sweep = if ccw_v <= ccw_1 {
        ccw_1
    } else {
        -(std::f64::consts::TAU - ccw_1)
    };
    Ok(Arc2 {
        cx,
        cy,
        r,
        a0,
        sweep,
    })
}

/// The outline as a run of 2-D primitives: `Ok(Left)` lines by
/// endpoints, `Ok(Right)` arcs — the single source both compilers
/// walk. Coordinates are offset by `(dx, dy)`.
enum Piece {
    Line((f64, f64), (f64, f64)),
    Arc(Arc2),
}

fn positive(what: &'static str, value: f64) -> Result<f64, ProfileError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(ProfileError::NonPositiveDimension { what, value })
    }
}

/// Decompose the profile into pieces, offset by `(dx, dy)`. A circle
/// yields a single full-turn arc.
fn pieces(profile: &AnalyticProfile<f64>, dx: f64, dy: f64) -> Result<Vec<Piece>, ProfileError> {
    let sh = |p: &[f64; 2]| (p[0] + dx, p[1] + dy);
    let polygon = |pts: Vec<(f64, f64)>| -> Vec<Piece> {
        (0..pts.len())
            .map(|i| Piece::Line(pts[i], pts[(i + 1) % pts.len()]))
            .collect()
    };
    Ok(match &profile.shape {
        ProfileShape::Rect { width, height } => {
            let hw = positive("rect width", *width)? / 2.0;
            let hh = positive("rect height", *height)? / 2.0;
            polygon(vec![
                (dx - hw, dy - hh),
                (dx + hw, dy - hh),
                (dx + hw, dy + hh),
                (dx - hw, dy + hh),
            ])
        }
        ProfileShape::Circle { radius } => vec![Piece::Arc(Arc2 {
            cx: dx,
            cy: dy,
            r: positive("circle radius", *radius)?,
            a0: 0.0,
            sweep: std::f64::consts::TAU,
        })],
        ProfileShape::Polygon { points } => {
            if points.len() < 3 {
                return Err(ProfileError::PolygonTooShort { got: points.len() });
            }
            polygon(points.iter().map(sh).collect())
        }
        ProfileShape::Path { start, segments } => {
            if segments.len() < 2 {
                return Err(ProfileError::PathTooShort);
            }
            let start = sh(start);
            let mut cur = start;
            let mut out = Vec::with_capacity(segments.len() + 1);
            for seg in segments {
                match seg {
                    PathSegment::Line { to } => {
                        let to = sh(to);
                        out.push(Piece::Line(cur, to));
                        cur = to;
                    }
                    PathSegment::Arc { via, to } => {
                        let to = sh(to);
                        out.push(Piece::Arc(three_point_arc(cur, sh(via), to)?));
                        cur = to;
                    }
                }
            }
            let gap = ((cur.0 - start.0).powi(2) + (cur.1 - start.1).powi(2)).sqrt();
            if gap > 1.0e-12 {
                out.push(Piece::Line(cur, start));
            }
            out
        }
    })
}

impl AnalyticProfile<f64> {
    /// The exact closed wire in the profile's LOCAL frame — centred at
    /// the plane origin, placement ignored (callers that carry the
    /// placement elsewhere, e.g. through a sweep path).
    ///
    /// # Errors
    ///
    /// [`ProfileError`] for a description that encloses no area or
    /// carries a non-positive dimension.
    pub fn local_wire(&self, plane: &SketchPlane) -> Result<NurbsCurve, ProfileError> {
        wire(&pieces(self, 0.0, 0.0)?, &self.shape, plane)
    }

    /// The exact closed wire in the sketch plane's own 2-D frame
    /// (`x` along `u_axis`, `y` along `v_axis`, `z = 0`), unplaced and
    /// unlifted — the form the sweep kernel's wire contract takes
    /// (profile `x` → frame normal, `y` → binormal). A caller that
    /// sweeps on any plane but XY sweeps this wire along local `z` and
    /// places the body with the plane's rigid frame afterwards; lifting
    /// the wire first collapses every span along the plane's `v_axis`
    /// (PK-416: a rectangle on YZ extruded to two zero-normal faces).
    ///
    /// # Errors
    ///
    /// As [`Self::local_wire`].
    pub fn flat_wire(&self) -> Result<NurbsCurve, ProfileError> {
        wire_2d(&pieces(self, 0.0, 0.0)?, &self.shape)
    }

    /// The exact closed wire with the placement applied.
    ///
    /// # Errors
    ///
    /// As [`Self::local_wire`].
    pub fn placed_wire(&self, plane: &SketchPlane) -> Result<NurbsCurve, ProfileError> {
        wire(
            &pieces(self, self.center[0], self.center[1])?,
            &self.shape,
            plane,
        )
    }

    /// The outline loaded into a [`Sketch`] as a plain geometry
    /// container (placement applied): lines between shared points,
    /// arcs by centre/radius/angles, a circle as one full-turn arc.
    /// Zero constraints; never solve it — it exists only to drive the
    /// tracer/builder pipeline.
    ///
    /// # Errors
    ///
    /// As [`Self::local_wire`].
    pub fn to_sketch_elements(&self) -> Result<Sketch, ProfileError> {
        let mut sk = Sketch::new();
        let mut last_point: Option<((f64, f64), popkern_constraint::PointHandle)> = None;
        let mut first_point: Option<((f64, f64), popkern_constraint::PointHandle)> = None;
        let same =
            |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12;
        for piece in pieces(self, self.center[0], self.center[1])? {
            match piece {
                Piece::Line(a, b) => {
                    let pa = match last_point {
                        Some((p, h)) if same(p, a) => h,
                        _ => sk.add_point(a.0, a.1),
                    };
                    if first_point.is_none() {
                        first_point = Some((a, pa));
                    }
                    let pb = match first_point {
                        Some((p, h)) if same(p, b) => h,
                        _ => sk.add_point(b.0, b.1),
                    };
                    sk.add_line(pa, pb);
                    last_point = Some((b, pb));
                }
                Piece::Arc(arc) => {
                    // Elements are always spelled counter-clockwise; the
                    // tracer resolves traversal direction from geometry.
                    let (a0, a1) = if arc.sweep >= 0.0 {
                        (arc.a0, arc.a0 + arc.sweep)
                    } else {
                        (arc.a0 + arc.sweep, arc.a0)
                    };
                    sk.add_arc(arc.cx, arc.cy, arc.r, a0, a1);
                    last_point = None;
                }
            }
        }
        Ok(sk)
    }

    /// Compile the profile to a planar sheet body on `plane` — the same
    /// lamina the constraint dialect produces (front face
    /// `Generated/0`, exact `Line` / `Arc` edge curves, exact rational
    /// pcurves). With `tracking`, element and entity names are minted
    /// under the context's operation exactly as for a constraint
    /// sketch.
    ///
    /// # Errors
    ///
    /// [`SketchError::Profile`] for an uncompilable description;
    /// otherwise as [`SketchOnPlane::to_brep`].
    pub fn to_brep(
        &self,
        plane: SketchPlane,
        tolerance: Tolerance,
        tracking: Option<&mut TrackingContext>,
    ) -> Result<SketchBrepResult, SketchError> {
        let mut on_plane = SketchOnPlane {
            sketch: self.to_sketch_elements()?,
            plane,
            element_names: Default::default(),
        };
        match tracking {
            Some(ctx) => {
                on_plane.assign_element_names(ctx.operation);
                on_plane.to_brep(tolerance, Some(ctx))
            }
            None => on_plane.to_brep(tolerance, None),
        }
    }
}

/// Assemble the pieces into one closed wire on `plane`.
fn wire(
    pieces: &[Piece],
    shape: &ProfileShape<f64>,
    plane: &SketchPlane,
) -> Result<NurbsCurve, ProfileError> {
    let mut curve = wire_2d(pieces, shape)?;
    // Lifted by the plane's affine map — exact for a rational curve.
    for p in &mut curve.control_points {
        *p = plane.to_3d(p.x, p.y);
    }
    Ok(curve)
}

/// Assemble the pieces into one closed wire in the plane's 2-D frame
/// (`z = 0`). Rect and polygon: degree-1 wires; circle: the standard
/// form; path: one quadratic span per piece.
fn wire_2d(pieces: &[Piece], shape: &ProfileShape<f64>) -> Result<NurbsCurve, ProfileError> {
    let curve = match shape {
        ProfileShape::Rect { .. } | ProfileShape::Polygon { .. } => {
            let pts: Vec<(f64, f64)> = pieces
                .iter()
                .filter_map(|p| match p {
                    Piece::Line(a, _) => Some(*a),
                    Piece::Arc(_) => None,
                })
                .collect();
            polygon_wire(&pts, 0.0)
        }
        ProfileShape::Circle { .. } => match pieces.first() {
            Some(Piece::Arc(a)) => circle_wire(a.cx, a.cy, a.r, 0.0),
            _ => return Err(ProfileError::PathTooShort),
        },
        ProfileShape::Path { .. } => {
            let mut spans: Vec<QuadraticSpan> = Vec::new();
            for piece in pieces {
                match piece {
                    Piece::Line(a, b) => spans.push(line_span(
                        Point3::new(a.0, a.1, 0.0),
                        Point3::new(b.0, b.1, 0.0),
                        0.0,
                        0.0,
                    )),
                    Piece::Arc(arc) => {
                        let (a0, a1) = (
                            arc.a0.min(arc.a0 + arc.sweep),
                            arc.a0.max(arc.a0 + arc.sweep),
                        );
                        let mut run = arc_spans(arc.cx, arc.cy, arc.r, a0, a1, 0.0);
                        if arc.sweep < 0.0 {
                            // Traverse clockwise: reverse the span run.
                            run.reverse();
                            for s in &mut run {
                                std::mem::swap(&mut s.p0, &mut s.p2);
                            }
                        }
                        spans.extend(run);
                    }
                }
            }
            // Uniform knot labels, one interval per span on [0, 1].
            let n = spans.len() as f64;
            for (i, s) in spans.iter_mut().enumerate() {
                s.t0 = i as f64 / n;
                s.t1 = (i as f64 + 1.0) / n;
            }
            assemble_quadratic(&spans)
        }
    };
    Ok(curve)
}

#[cfg(test)]
#[path = "profile_compile_tests.rs"]
mod tests;
