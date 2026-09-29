//! The variational-fit sides of a setback corner (PK-574): one
//! [`VariationalCapSide`] per setback arc, from the extractor's exact
//! section curve (or the saddle bridge) and the surface the cap meets
//! there, in the loop order the core fit expects (counter-clockwise
//! about the cap's outward normal).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::variational_cap::VariationalCapSide;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::kp20_paper_border_extract::Kp20SetbackArc;
use crate::fillet_surface::FilletSurfaceResult;

/// How to orient each side's neighbour: `outward(side_index, point,
/// natural_normal)` returns whether the natural normal is outward.
pub(crate) type OutwardTest<'a> = dyn Fn(usize, Point3, Vec3) -> bool + 'a;

/// Build the fit sides from the extractor's arcs (angular-walk order).
/// Each side runs `corner_prev` → `corner_next`; the list is reversed
/// (and every curve with it) when the arcs wind clockwise about the
/// outward normal.
pub(crate) fn fit_sides(
    arcs: &[Kp20SetbackArc],
    fillet_surfaces: &[FilletSurfaceResult],
    outward: &OutwardTest<'_>,
) -> Result<Vec<VariationalCapSide>, String> {
    let mut sides = Vec::with_capacity(arcs.len());
    let mut outward_sum = Vec3::new(0.0, 0.0, 0.0);
    for (k, arc) in arcs.iter().enumerate() {
        let (curve, neighbour) = match &arc.bridge {
            Some(bridge) => (bridge.curve.clone(), bridge.plane.clone()),
            None => (
                arc.section_curve(),
                fillet_surfaces
                    .get(arc.fillet_index)
                    .ok_or_else(|| {
                        format!("arc {k}: fillet index {} out of range", arc.fillet_index)
                    })?
                    .surface
                    .clone(),
            ),
        };
        let (t0, t1) = curve.parameter_range();
        let mid = curve.point_at(0.5 * (t0 + t1));
        let (u, v) = project_point_to_uv(&neighbour, &mid)
            .ok_or_else(|| format!("arc {k}: side midpoint does not project onto its neighbour"))?;
        let natural = neighbour.normal_at(u, v);
        let is_outward = outward(k, mid, natural);
        outward_sum = outward_sum + if is_outward { natural } else { natural * -1.0 };
        sides.push(VariationalCapSide {
            curve,
            neighbour,
            neighbour_reversed: !is_outward,
        });
    }
    // Newell normal of the corner loop.
    let corners: Vec<Point3> = sides
        .iter()
        .map(|s| s.curve.point_at(s.curve.parameter_range().0))
        .collect();
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    for (i, a) in corners.iter().enumerate() {
        let b = corners[(i + 1) % corners.len()];
        newell = newell
            + Vec3::new(
                (a.y - b.y) * (a.z + b.z),
                (a.z - b.z) * (a.x + b.x),
                (a.x - b.x) * (a.y + b.y),
            );
    }
    if newell.dot(&outward_sum) < 0.0 {
        sides.reverse();
        for side in &mut sides {
            side.curve = reversed_curve(&side.curve);
        }
    }
    Ok(sides)
}

/// The same curve traversed backwards over the same parameter range.
fn reversed_curve(curve: &NurbsCurve) -> NurbsCurve {
    let (a, b) = curve.parameter_range();
    let knots = curve.knots.iter().rev().map(|k| a + b - k).collect();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points: curve.control_points.iter().rev().copied().collect(),
        weights: curve.weights.iter().rev().copied().collect(),
    }
}
