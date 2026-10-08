//! The served result's last geometric gate (PK-706): every edge carries
//! a curve. The curve-surface validator skips edges with none, so a
//! curveless edge passed every check while the mesher chorded it and
//! opened the display mesh (measured: a post-processing rim weld left
//! one on a curved face of a 300-long piano hinge leaf; its root was
//! PK-707's ray window, but a served body must never carry one).
//!
//! An edge between two planar faces lies on their common line, which IS
//! the chord between its ends: it is given that exact line (a repair,
//! traced). An edge touching any curved face has no geometry to recover
//! here, and the result is refused.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;

/// Give curveless edges between two planes their exact line; refuse
/// `result` when a curveless edge touches a curved face.
///
/// # Errors
///
/// [`BoolError::ValidationFailed`] naming how many edges stay curveless.
pub(crate) fn repair_or_decline_curveless_edges(
    result: &mut BRep,
    tolerance: &Tolerance,
) -> Result<(), BoolError> {
    let planar = |b: &BRep, he| {
        let face = b.loops.get(b.half_edges[he].loop_).map(|l| l.face);
        face.and_then(|f| b.faces.get(f))
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
    };
    let bare: Vec<_> = result
        .edges
        .iter()
        .filter(|(_, e)| e.curve.is_none())
        .map(|(h, _)| h)
        .collect();
    let mut left = 0;
    let mut repaired = 0;
    for eh in bare {
        let he = result.edges[eh].half_edge;
        let twin = result.half_edges[he].twin;
        let line = (planar(result, he) && planar(result, twin))
            .then(|| {
                let a = result.vertices[result.half_edges[he].origin].position;
                let b = result.vertices[result.half_edges[twin].origin].position;
                CurveOnEdge::line_segment(a, b, tolerance.linear())
            })
            .flatten();
        match line {
            Some(curve) => {
                result.edges[eh].curve = Some(curve);
                repaired += 1;
            }
            None => left += 1,
        }
    }
    if repaired > 0 && popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "[boolean] post-processing: {repaired} curveless planar edge(s) given their line"
        );
    }
    if left > 0 {
        return Err(BoolError::ValidationFailed(format!(
            "{left} edge(s) of the result carry no curve geometry on a curved face — \
             declining rather than serving a body the mesher would chord open"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::repair_or_decline_curveless_edges;
    use popkern_core::geometry::curve::CurveGeometry;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn a_planar_edge_gets_its_line_and_a_curved_one_is_refused() {
        let tol = Tolerance::default();
        let mut cube = BRep::unit_cube(tol).expect("cube");
        let edge = cube.edges.keys().next().expect("an edge");
        cube.edges[edge].curve = None;
        repair_or_decline_curveless_edges(&mut cube, &tol).expect("repaired");
        assert!(matches!(
            cube.edges[edge].curve.as_ref().map(|c| &c.geometry),
            Some(CurveGeometry::Line { .. })
        ));
        let mut rod = BRep::cylinder(1.0, 2.0, tol).expect("rod");
        // Any edge on the cylinder's curved wall.
        let rim = rod
            .edges
            .iter()
            .find(|(_, e)| {
                let face = rod.loops[rod.half_edges[e.half_edge].loop_].face;
                rod.faces[face].surface.as_ref().is_some_and(|s| {
                    !matches!(
                        s.geometry,
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
                    )
                })
            })
            .map(|(h, _)| h)
            .expect("a wall edge");
        rod.edges[rim].curve = None;
        let err = repair_or_decline_curveless_edges(&mut rod, &tol).expect_err("refused");
        assert!(err.to_string().contains("1 edge(s)"), "{err}");
    }
}
