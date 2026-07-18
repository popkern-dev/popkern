//! PK-142 stage 1 — open sheet body for a marched face–face blend
//! strip. One face, four boundary edges (two contact curves, two end
//! arcs), boundary-registered, with the skinned strip surface and
//! chart-border pcurves — one geometric authority from construction
//! (the PK-148 lesson).

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{CurveEval, NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::contact::ContactSolution;
use crate::error::BlendError;

pub(crate) const ARC_SAMPLES: usize = 9;

/// A boundary polyline paired with its chart border endpoints.
type BoundarySide<'a> = (&'a [Point3], ((f64, f64), (f64, f64)));

/// Strip skin surface plus its four boundary polylines in chart order
/// (`side_a` = v=0 border, `side_b` = v=1, arcs = u=0 / u=1 borders).
pub(crate) struct StripSkin {
    /// Skinned strip surface (degree-1 polyline skin — stage-1
    /// geometric authority; sections and borders interpolate).
    pub surface: NurbsSurface,
    /// Contact polyline on face A (v = 0 border).
    pub side_a: Vec<Point3>,
    /// Contact polyline on face B (v = 1 border).
    pub side_b: Vec<Point3>,
    /// First section arc (u = 0 border), from side A to side B.
    pub start_arc: Vec<Point3>,
    /// Last section arc (u = 1 border), from side A to side B.
    pub end_arc: Vec<Point3>,
}

/// Skin the section arcs of a marched chain into the strip surface
/// and extract its boundary polylines. Shared by the sheet emitter
/// and the TrimBoth attach surgery so both use ONE geometric
/// authority for the strip. `bulge` resolves the antipodal
/// (half-circle) section ambiguity of line-contact (parallel-support)
/// blends — the arc bows toward it.
pub(crate) fn strip_skin(chain: &[ContactSolution], radius: f64, bulge: Option<Vec3>) -> StripSkin {
    let n_sections = chain.len();
    let mut control_points = Vec::with_capacity(n_sections * ARC_SAMPLES);
    for sol in chain {
        let arc = section_arc_points(sol, radius, bulge);
        control_points.extend_from_slice(&arc);
    }
    let degree_u = 1.min(n_sections - 1);
    let degree_v = 1.min(ARC_SAMPLES - 1);
    let surface = NurbsSurface {
        degree_u,
        degree_v,
        knots_u: clamped_uniform_knots(n_sections, degree_u),
        knots_v: clamped_uniform_knots(ARC_SAMPLES, degree_v),
        control_points: control_points.clone(),
        num_u: n_sections,
        num_v: ARC_SAMPLES,
        weights: vec![1.0; n_sections * ARC_SAMPLES],
    };
    let side_a: Vec<Point3> = (0..n_sections)
        .map(|i| control_points[i * ARC_SAMPLES])
        .collect();
    let side_b: Vec<Point3> = (0..n_sections)
        .map(|i| control_points[i * ARC_SAMPLES + ARC_SAMPLES - 1])
        .collect();
    let start_arc: Vec<Point3> = control_points[..ARC_SAMPLES].to_vec();
    let end_arc: Vec<Point3> = control_points[(n_sections - 1) * ARC_SAMPLES..].to_vec();
    StripSkin {
        surface,
        side_a,
        side_b,
        start_arc,
        end_arc,
    }
}

/// Build the strip sheet body from marched contact solutions.
pub(crate) fn strip_sheet_body(
    chain: &[ContactSolution],
    radius: f64,
    bulge: Option<Vec3>,
    tolerance: Tolerance,
) -> Result<BRep, BlendError> {
    let wrap = |reason: String| BlendError::InvalidParameter {
        parameter: "face_blend".into(),
        value: "strip sheet".into(),
        reason,
    };

    let skin = strip_skin(chain, radius, bulge);
    let StripSkin {
        surface,
        side_a,
        side_b,
        start_arc,
        end_arc,
    } = skin;

    let corner_a0 = side_a[0];
    let corner_a1 = *side_a.last().ok_or_else(|| wrap("empty side".into()))?;
    let corner_b0 = side_b[0];
    let corner_b1 = *side_b.last().ok_or_else(|| wrap("empty side".into()))?;

    // Topology: quad sheet a0 → a1 → b1 → b0 (start arc runs a0→b0,
    // end arc a1→b1; loop order a0→a1 (contact A), a1→b1 (end arc),
    // b1→b0 (contact B reversed), b0→a0 (start arc reversed)).
    let mut sheet = BRep::new(tolerance);
    let mvfs = sheet.mvfs(corner_a0).map_err(|e| wrap(format!("{e}")))?;
    let mev1 = sheet
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, corner_a1)
        .map_err(|e| wrap(format!("{e}")))?;
    let mev2 = sheet
        .mev(mev1.half_edge_twin, corner_b1)
        .map_err(|e| wrap(format!("{e}")))?;
    let mev3 = sheet
        .mev(mev2.half_edge_twin, corner_b0)
        .map_err(|e| wrap(format!("{e}")))?;
    let mef = sheet
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .map_err(|e| wrap(format!("{e}")))?;

    // popkern sheets are TWO-SIDED (sheet_polygon pattern): mef
    // leaves a front face (the original) and a backside face (the
    // mef result). The strip surface goes on both, opposite
    // orientations; every edge is interior (both half-edges live).
    let back = mef.face;
    let front = sheet
        .faces
        .keys()
        .find(|f| *f != back)
        .ok_or_else(|| wrap("front face missing after mef".into()))?;
    sheet.faces[front].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(surface.clone()),
        orientation: true,
    });
    sheet.faces[back].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(surface),
        orientation: false,
    });

    // Curves + chart-border pcurves assigned by ENDPOINT MATCHING —
    // the front loop's start half-edge and direction after mef are
    // an Euler-op detail, not a contract (assuming an order left a
    // slot-width OpenLoop gap). Each boundary polyline knows its
    // chart border; orientation follows each edge's PRIMARY
    // half-edge, pcurves follow each half-edge's own direction.
    let sides: [BoundarySide; 4] = [
        (&side_a, ((0.0, 0.0), (1.0, 0.0))),
        (&end_arc, ((1.0, 0.0), (1.0, 1.0))),
        (&side_b, ((0.0, 1.0), (1.0, 1.0))),
        (&start_arc, ((0.0, 0.0), (0.0, 1.0))),
    ];
    let tol_match = tolerance.linear().max(1.0e-9) * 100.0;
    let front_loop = sheet.faces[front].outer_loop;
    let start_he = sheet.loops[front_loop].half_edge;
    let mut he = start_he;
    for _ in 0..4 {
        let he_data = sheet.half_edges[he].clone();
        let origin = sheet.vertices[he_data.origin].position;
        let dest = sheet.vertices[sheet.half_edges[he_data.twin].origin].position;
        let (points, border) = sides
            .iter()
            .find_map(|(pts, border)| {
                let first = pts[0];
                let last = *pts.last()?;
                if first.distance_to(&origin) < tol_match && last.distance_to(&dest) < tol_match {
                    Some((pts.to_vec(), (border.0, border.1)))
                } else if first.distance_to(&dest) < tol_match
                    && last.distance_to(&origin) < tol_match
                {
                    Some((reversed(pts), (border.1, border.0)))
                } else {
                    None
                }
            })
            .ok_or_else(|| wrap(format!("no boundary polyline matches edge at {origin:?}")))?;
        // points now run origin→dest of THIS half-edge; orient the
        // stored curve along the edge's primary half-edge.
        let primary_forward = sheet.edges[he_data.edge].half_edge == he;
        let curve_points = if primary_forward {
            points.clone()
        } else {
            reversed(&points)
        };
        sheet.edges[he_data.edge].curve = Some(polyline_curve(&curve_points).map_err(wrap)?);
        sheet.half_edges[he].pcurve = Some(Pcurve::linear(border.0, border.1, 0.0, 1.0));
        sheet.half_edges[he_data.twin].pcurve = Some(Pcurve::linear(border.1, border.0, 0.0, 1.0));
        he = he_data.next;
    }
    let _ = radius;
    Ok(sheet)
}

fn reversed(points: &[Point3]) -> Vec<Point3> {
    points.iter().rev().copied().collect()
}

fn polyline_curve(points: &[Point3]) -> Result<CurveOnEdge, String> {
    let nurbs: NurbsCurve = crate::face_blend::polyline_nurbs(points)?;
    let (t0, t1) = nurbs.parameter_range();
    Ok(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: t0,
        t_end: t1,
    })
}

pub(crate) fn section_arc_points(
    sol: &ContactSolution,
    radius: f64,
    bulge: Option<Vec3>,
) -> Vec<Point3> {
    // Arc from contact0 to contact1 about the ball centre, minor
    // sweep, explicit basis (no frame convention assumed).
    let u = sol.contact0 - sol.center;
    let to_end = sol.contact1 - sol.center;
    let normal = u.cross(&to_end);
    let mut points = Vec::with_capacity(ARC_SAMPLES);
    match normal.normalized() {
        None => {
            // Antipodal contacts: the section is a half-circle and
            // its bow side is not determined by the two contacts
            // alone — use the caller's bulge hint. Without one, fall
            // back to a straight chord (degenerate/legacy behavior).
            let half_circle = bulge
                .and_then(|d| u.normalized().map(|uh| d - uh * d.dot(&uh)))
                .and_then(|d| d.normalized());
            if let Some(v_dir) = half_circle {
                let v = v_dir * u.length();
                for k in 0..ARC_SAMPLES {
                    let theta = std::f64::consts::PI * k as f64 / (ARC_SAMPLES - 1) as f64;
                    points.push(sol.center + u * theta.cos() + v * theta.sin());
                }
                let _ = radius;
                return points;
            }
            for k in 0..ARC_SAMPLES {
                let t = k as f64 / (ARC_SAMPLES - 1) as f64;
                points.push(Point3::new(
                    sol.contact0.x + (sol.contact1.x - sol.contact0.x) * t,
                    sol.contact0.y + (sol.contact1.y - sol.contact0.y) * t,
                    sol.contact0.z + (sol.contact1.z - sol.contact0.z) * t,
                ));
            }
        }
        Some(n) => {
            let v = n.cross(&u);
            let sweep = to_end.dot(&v).atan2(to_end.dot(&u));
            let sweep = if sweep >= 0.0 {
                sweep
            } else {
                sweep + std::f64::consts::TAU
            };
            for k in 0..ARC_SAMPLES {
                let t = sweep * k as f64 / (ARC_SAMPLES - 1) as f64;
                points.push(sol.center + u * t.cos() + v * t.sin());
            }
        }
    }
    let _ = radius;
    points
}

pub(crate) fn clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let interior = n - degree - 1;
    let mut knots = vec![0.0; degree + 1];
    for i in 1..=interior {
        knots.push(i as f64 / (interior + 1) as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}
