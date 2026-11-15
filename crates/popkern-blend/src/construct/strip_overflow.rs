//! Fillet overflow onto a neighbouring strip (PK-1189, design
//! `docs/architecture/fillet-overflow-onto-neighbour-strip.md`).
//!
//! A concave strip whose foot on a shared support face F runs past a
//! convex strip's rail into its band rolls on: past the rail F has ended
//! and the convex strip continues it tangentially, so the concave strip's
//! ball keeps touching its other support W and rolls onto the convex
//! strip S. This module plans that overflow before the apply phase — the
//! foot's crossings of the rail, the ball solved against W and S in each
//! section plane of the concave strip, and the overflow surface skinned
//! through those exact sections — so every geometric decline happens
//! before topology changes. `strip_overflow_surgery` splices it into the
//! applied body.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::strip_overflow_sections::{
    FootPath, arc_section, g2_ramp, g2_section, rolling_section, skin_sections, solve_section,
};
use crate::analysis::Convexity;
use crate::construct::EdgeGeometry;

/// Sections first solved across the overflow (the skin interpolates
/// them cubically and doubles them until it meets the exact sections
/// between, up to `MAX_STATIONS`).
const STATIONS: usize = 17;

/// The most sections an overflow skin takes.
const MAX_STATIONS: usize = 513;

/// Samples of the foot when locating its crossings of the rail.
const FOOT_SAMPLES: usize = 512;

/// One concave strip rolling over a convex one, planned on the input
/// body and spliced into the applied one.
#[derive(Debug, Clone)]
pub(crate) struct OverflowPlan {
    /// Spec index of the concave (rolling) strip.
    pub rolling: usize,
    /// Spec index of the convex strip it rolls onto.
    pub onto: usize,
    /// The shared support face F (its handle survives the apply phase).
    pub face_f: FaceHandle,
    /// The rolling strip's other support W.
    pub face_w: FaceHandle,
    /// F's surface.
    pub support_surface: SurfaceGeometry,
    /// The rolling strip's surface (a torus or cylinder).
    pub rolling_surface: SurfaceGeometry,
    /// The convex strip's surface.
    pub onto_surface: SurfaceGeometry,
    /// The overflow surface: `u` across the overflow (P1 side at 0), `v`
    /// across the section (on S at 0, on W at 1).
    pub surface: NurbsSurface,
    /// Where the foot crosses the rail entering the band (on F and S).
    pub p1: Point3,
    /// Where the foot crosses the rail leaving the band.
    pub p2: Point3,
    /// The overflow's end on W at P1's section.
    pub q1: Point3,
    /// The overflow's end on W at P2's section.
    pub q2: Point3,
    /// G2 sections (degree 5, matching each support's curvature) instead
    /// of the rolling ball's circular arcs.
    pub g2: bool,
}

impl OverflowPlan {
    /// The overflow's contact curve on S, P1 to P2.
    pub fn curve_on_onto(&self) -> NurbsCurve {
        self.surface.iso_curve_at_v(0.0)
    }
    /// The overflow's contact curve on W, Q1 to Q2.
    pub fn curve_on_wall(&self) -> NurbsCurve {
        self.surface.iso_curve_at_v(1.0)
    }
    /// The section at P1, from P1 to Q1 (the rolling strip's own iso-curve).
    pub fn section_start(&self) -> NurbsCurve {
        self.surface.iso_curve_at_u(0.0)
    }
    /// The section at P2, from P2 to Q2.
    pub fn section_end(&self) -> NurbsCurve {
        self.surface.iso_curve_at_u(1.0)
    }
}

/// Plan the overflow of an interfering pair whose concave strip rolls
/// over its convex one; `Err` says why the pair does not.
pub(crate) fn plan_overflow(
    edge_geometries: &[EdgeGeometry],
    brep: &BRep,
    pair: (usize, usize),
    continuity: crate::types::ContinuityLevel,
    tolerance: &Tolerance,
) -> Result<OverflowPlan, String> {
    let g2 = continuity == crate::types::ContinuityLevel::G2;
    let (i, j) = pair;
    let (rolling, onto) = match (
        edge_geometries[i].analysis.convexity,
        edge_geometries[j].analysis.convexity,
    ) {
        (Convexity::Concave, Convexity::Convex) => (i, j),
        (Convexity::Convex, Convexity::Concave) => (j, i),
        _ => return Err("not a concave strip against a convex one".into()),
    };
    let t = &edge_geometries[rolling];
    let s = &edge_geometries[onto];
    // The shared support face and each strip's trim on it (with its
    // exact circle when the contact is circular).
    let sides_t = [
        (
            t.analysis.face_1,
            &t.trim_3d_1,
            &t.analysis.surface_1,
            t.analysis.face_2,
            &t.analysis.surface_2,
        ),
        (
            t.analysis.face_2,
            &t.trim_3d_2,
            &t.analysis.surface_2,
            t.analysis.face_1,
            &t.analysis.surface_1,
        ),
    ];
    let sides_s = [
        (s.analysis.face_1, &s.trim_3d_1, s.trim_circle_1),
        (s.analysis.face_2, &s.trim_3d_2, s.trim_circle_2),
    ];
    let (face_f, foot_trim, f_surface, face_w, wall, rail, rail_circle) = sides_t
        .iter()
        .find_map(|&(ft, foot, fsurf, fw, wall)| {
            sides_s
                .iter()
                .find(|(fs, _, _)| *fs == ft)
                .map(|&(_, rail, circle)| (ft, foot, fsurf, fw, wall, rail, circle))
        })
        .ok_or("the strips share no support face")?;
    // The rolling strip: a torus or cylinder (G1), or the swept NURBS of
    // G2 sections, read through its edge (a circle or line).
    let r = t.effective_radius;
    let path = match &t.surface {
        SurfaceGeometry::Torus { .. } | SurfaceGeometry::Cylinder { .. } if !g2 => {
            FootPath::of(&t.surface, foot_trim, tolerance)?
        }
        SurfaceGeometry::Nurbs(_) if g2 => {
            let edge = brep
                .edges
                .get(t.analysis.edge)
                .and_then(|e| e.curve.as_ref())
                .ok_or("the rolling strip's edge has no curve")?;
            FootPath::of_edge(edge, foot_trim, (wall, f_surface), r, tolerance)?
        }
        _ => {
            return Err(
                "the rolling strip is neither a torus or cylinder (G1) nor G2 sections".into(),
            );
        }
    };
    let foot_at = |u: f64| -> Result<Point3, String> {
        let spine = path.spine(u);
        let uv = (f_surface as &dyn SurfaceEval)
            .closest_point(&spine)
            .ok_or("the spine does not project onto F")?;
        Ok(f_surface.point_at(uv.0, uv.1))
    };
    let rail_point = |p: &Point3| -> Point3 {
        match rail_circle {
            Some((c, n, radius)) => {
                let d = *p - c;
                match (d - n * d.dot(&n)).normalized() {
                    Some(radial) => c + radial * radius,
                    None => crate::center_curve_station::closest_point_on_nurbs(rail, p),
                }
            }
            None => crate::center_curve_station::closest_point_on_nurbs(rail, p),
        }
    };
    let onto_edge = brep
        .edges
        .get(s.analysis.edge)
        .and_then(|e| e.curve.as_ref())
        .ok_or("the convex strip's edge has no curve")?;
    let edge_dist = |p: &Point3| -> (f64, f64) {
        let tt = popkern_core::geometry::curve::closest_point_in_window(
            &onto_edge.geometry,
            p,
            onto_edge.t_start,
            onto_edge.t_end,
        );
        (tt, onto_edge.geometry.point_at(tt).distance_to(p))
    };
    // How far into S's band a point of F lies: the band's width at the
    // rail point nearest it less the point's distance to S's edge.
    let into_band = |p: &Point3| -> f64 { edge_dist(&rail_point(p)).1 - edge_dist(p).1 };
    let g_at = |u: f64| -> Result<f64, String> { Ok(into_band(&foot_at(u)?)) };
    let param = |k: usize| k as f64 / FOOT_SAMPLES as f64;
    let g: Vec<f64> = (0..=FOOT_SAMPLES)
        .map(|k| g_at(param(k)))
        .collect::<Result<_, _>>()?;
    // Sign changes, refined by bisection: entering (g rising through 0)
    // and leaving.
    let refine = |mut a: f64, mut b: f64| -> Result<f64, String> {
        let ga = g_at(a)?;
        for _ in 0..60 {
            let m = 0.5 * (a + b);
            if (g_at(m)? > 0.0) == (ga > 0.0) {
                a = m;
            } else {
                b = m;
            }
        }
        Ok(0.5 * (a + b))
    };
    let mut enter = Vec::new();
    let mut leave = Vec::new();
    for k in 0..FOOT_SAMPLES {
        if g[k] <= 0.0 && g[k + 1] > 0.0 {
            enter.push(refine(param(k), param(k + 1))?);
        } else if g[k] > 0.0 && g[k + 1] <= 0.0 {
            leave.push(refine(param(k), param(k + 1))?);
        }
    }
    let (t_enter, t_leave) = match (enter.as_slice(), leave.as_slice()) {
        ([a], [b]) => (*a, *b),
        _ => {
            return Err(format!(
                "the foot crosses the rail {} times entering and {} leaving (one each is rolled over)",
                enter.len(),
                leave.len()
            ));
        }
    };
    // The overflow runs forward from entering to leaving, through the
    // foot's closure on a closed ring.
    let span = if t_leave > t_enter {
        t_leave - t_enter
    } else if path.closed() {
        t_leave - t_enter + 1.0
    } else {
        return Err("the overflow runs off the foot's end".into());
    };
    let at = |f: f64| {
        let tt = t_enter + span * f;
        if tt > 1.0 { tt - 1.0 } else { tt }
    };
    // The exact section at a fraction of the overflow, its S contact
    // checked on the convex strip.
    let section = |f: f64| -> Result<crate::contact::ContactSolution, String> {
        let p = foot_at(at(f))?;
        let (spine, plane_normal) = rolling_section(&path, &p)?;
        let solved = solve_section(wall, &s.surface, r, &spine, &plane_normal, tolerance)?;
        // The S contact must lie on the convex strip: no farther from
        // its rail on F than the strip's own width.
        let rail_pt = rail_point(&solved.contact1);
        if solved.contact1.distance_to(&rail_pt) > edge_dist(&rail_pt).1 * 1.5 {
            return Err(format!(
                "section at {f:.3}: the ball leaves the convex strip (contact {:.3e} from its rail)",
                solved.contact1.distance_to(&rail_pt)
            ));
        }
        if f == 0.0 || f == 1.0 {
            // At the crossings the ball is the rolling strip's own.
            let off = solved.center.distance_to(&spine);
            if off > (tolerance.linear() * 10.0).max(1e-9) {
                return Err(format!(
                    "end section: the solved ball is {off:.3e} off the strip's spine"
                ));
            }
        }
        Ok(solved)
    };
    // A section's curve from its S contact (v = 0) to its W contact
    // (v = 1) at a fraction `f` of the overflow: the ball's rational arc
    // (G1), or the degree-5 G2 section. At the overflow's ends the G2
    // section is read against F, so it is the rolling strip's own there:
    // F and S share the tangent plane but not the curvature, so the
    // support's curvature ramps from F's to S's over a fifth of the
    // overflow at each end (a section's interior control points are
    // linear in its end curvature, so the blend of the two sections is
    // the section of the blended curvature).
    let curve_of = |solved: &crate::contact::ContactSolution, f: f64| {
        if g2 {
            g2_section(solved, wall, (f_surface, &s.surface), g2_ramp(f))
        } else {
            arc_section(solved)
        }
    };
    // Sections doubled until the skin meets the exact section midway
    // between each pair within a quarter of the tolerance (17 sections
    // left 6.7e-6 between them on a 0.37-long overflow).
    let gate = tolerance.linear() * 0.25;
    let mut count = STATIONS;
    let (surface, ends) = loop {
        let solved: Vec<crate::contact::ContactSolution> = (0..count)
            .map(|k| section(k as f64 / (count - 1) as f64))
            .collect::<Result<_, _>>()?;
        let curves: Vec<NurbsCurve> = solved
            .iter()
            .enumerate()
            .map(|(k, x)| curve_of(x, k as f64 / (count - 1) as f64))
            .collect::<Result<_, _>>()?;
        let surface = skin_sections(&curves)?;
        let skin = SurfaceGeometry::Nurbs(surface.clone());
        let mut worst = 0.0_f64;
        for k in 0..count - 1 {
            let f = (k as f64 + 0.5) / (count - 1) as f64;
            let exact = curve_of(&section(f)?, f)?;
            let (e0, e1) = exact.parameter_range();
            for v in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let e = exact.point_at(e0 + (e1 - e0) * v);
                worst = worst.max(e.distance_to(&skin.point_at(f, v)));
            }
        }
        if worst <= gate {
            let ends: Vec<(Point3, Point3)> =
                solved.iter().map(|x| (x.contact1, x.contact0)).collect();
            break (surface, ends);
        }
        if count >= MAX_STATIONS {
            return Err(format!(
                "the overflow skin stays {worst:.3e} off its sections at {count} sections"
            ));
        }
        count = 2 * count - 1;
    };
    let (p1, q1) = ends[0];
    let (p2, q2) = ends[ends.len() - 1];
    Ok(OverflowPlan {
        rolling,
        onto,
        face_f,
        face_w,
        support_surface: f_surface.clone(),
        rolling_surface: t.surface.clone(),
        onto_surface: s.surface.clone(),
        surface,
        p1,
        p2,
        q1,
        q2,
        g2,
    })
}
