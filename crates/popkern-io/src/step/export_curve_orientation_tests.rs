//! PK-554: an exact edge curve stored opposite its edge exports as
//! itself, not as a polyline sampled from a surface trim.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::*;
use crate::step::export::ExportContext;

/// Every edge curve reversed by swapping its parameter bounds: the same
/// point set, traversed opposite the edge.
fn reversed_curves(mut brep: BRep) -> BRep {
    for edge in brep.edges.values_mut() {
        if let Some(curve) = edge.curve.as_mut() {
            std::mem::swap(&mut curve.t_start, &mut curve.t_end);
        }
    }
    brep
}

fn entity_counts(brep: &BRep) -> (usize, usize, usize) {
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(brep)
        .unwrap_or_else(|e| panic!("export: {e:?}"));
    let count = |name: &str| ctx.entities.iter().filter(|e| e.type_name == name).count();
    (
        count("LINE"),
        count("CIRCLE"),
        count("B_SPLINE_CURVE_WITH_KNOTS"),
    )
}

#[test]
fn a_reversed_edge_curve_agrees_with_its_trims() {
    // Co-parameterized edges (the cube's lines): before PK-554 a line
    // stored opposite its edge read its full length off its trims (40
    // on pocketed_tray) and was replaced by a sampled polyline.
    let reversed =
        reversed_curves(BRep::unit_cube(Tolerance::default()).unwrap_or_else(|e| panic!("{e}")));
    let mut measured = 0;
    for (edge, data) in reversed.edges.iter() {
        let Some(curve) = data.curve.as_ref() else {
            continue;
        };
        let gap = edge_curve_surface_trim_disagreement(&reversed, edge, curve)
            .unwrap_or_else(|| panic!("{edge:?}: trims measurable"));
        assert!(
            gap < 1e-9,
            "{edge:?}: reversed curve reads {gap:e} off its trims"
        );
        measured += 1;
    }
    assert_eq!(measured, 12);
}

#[test]
fn reversed_exact_curves_export_as_lines_and_circles() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    let cylinder = BRep::cylinder(1.5, 2.0, Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    for body in [cube, cylinder] {
        let native = entity_counts(&body);
        let reversed = entity_counts(&reversed_curves(body));
        assert_eq!(reversed, native, "(LINE, CIRCLE, B_SPLINE) counts");
    }
}

/// Degree-1 B-spline curves with more than two control points in
/// exported Part 21 text: sampled polylines standing in for curves.
fn polyline_curves(brep: &BRep) -> usize {
    let mut ctx = ExportContext::new();
    ctx.export_brep_items(brep)
        .unwrap_or_else(|e| panic!("export: {e:?}"));
    ctx.entities
        .iter()
        .filter(|e| e.type_name == "B_SPLINE_CURVE_WITH_KNOTS")
        .filter(|e| format!("{:?}", e.params).starts_with("[String(\"\"), Int(1), List(["))
        .filter(|e| {
            let text = format!("{:?}", e.params);
            text.matches("EntityRef").count() > 2
        })
        .count()
}

#[test]
fn a_rim_parameterized_differently_from_its_trim_keeps_its_exact_curve() {
    // The cylinder's rims are rational NURBS circles; their wall trims
    // are linear in angle, so they disagree at equal fractions by
    // 2.4e-2 while agreeing as point sets. Before PK-554 each rim was
    // rewritten as a 17-sample polyline (2.9e-2 sagitta).
    let cylinder = BRep::cylinder(1.5, 2.0, Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    for (edge, data) in cylinder.edges.iter() {
        let Some(curve) = data.curve.as_ref() else {
            continue;
        };
        let gap = edge_curve_surface_trim_geometric_gap(&cylinder, edge, curve)
            .unwrap_or_else(|| panic!("{edge:?}: trims measurable"));
        assert!(gap < 1e-9, "{edge:?}: {gap:e} off its trims as a point set");
    }
    assert_eq!(polyline_curves(&cylinder), 0);
    assert_eq!(polyline_curves(&reversed_curves(cylinder)), 0);
}
