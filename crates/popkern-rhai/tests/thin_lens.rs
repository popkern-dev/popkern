//! Two balls of radius 5 overlapping 1e-6 deep. At the default linear
//! tolerance (1e-6) their common part is thinner than the tolerance and
//! has no inside; a part declaring the tolerance its features need
//! (1e-8) gets the lens: two spherical caps meeting in a circle, each
//! cap's area 2 pi r h with h = 5e-7.

use popkern_core::tolerance::Tolerance;

fn lens(tolerance: &str) -> Option<popkern_core::topology::BRep> {
    let source = format!(
        "let m = model({tolerance});\nlet lens = m.sphere(5.0);\n\
         lens = lens.intersect(m.tool(|m| m.sphere(5.0).translate(9.999999, 0.0, 0.0)));\nm"
    );
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host.run(&source, "lens").ok()?;
    let parts = run.model?.to_parts_all().ok()?;
    parts.first().map(|p| p.brep().clone())
}

#[test]
fn a_lens_thinner_than_the_tolerance_has_no_inside() {
    let empty = lens("").is_none_or(|b| b.faces.is_empty());
    assert!(empty, "a 1e-6 lens at a 1e-6 tolerance has no inside");
}

#[test]
fn a_lens_at_a_tolerance_its_features_need_is_two_caps() {
    let body = lens("#{ tolerance: 1e-8 }").expect("the lens builds");
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    assert_eq!(body.faces.len(), 2);
    let area: f64 = body
        .faces
        .keys()
        .map(|f| popkern_measure::face_area(&body, f, &Tolerance::default()).unwrap())
        .sum();
    let want = 2.0 * (2.0 * std::f64::consts::PI * 5.0 * 5e-7);
    assert!((area - want).abs() < 1e-6 * want, "{area} want {want}");
}
