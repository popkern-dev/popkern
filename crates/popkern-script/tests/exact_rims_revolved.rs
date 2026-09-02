//! PK-418: cuts through revolved bodies leave exact circular rims on
//! sphere and torus faces too — a cap sliced off a sphere selects its
//! rim as `circular()`, a torus cut through its spine plane keeps both
//! equator rims and the pieces of both seams as arcs — and the bodies
//! stay watertight at the coarse end of the tolerance envelope, where
//! the mesher used to put the sphere face's hole rings on the wrong
//! period sheet.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_history::feature::AxisDesc;
use popkern_script::prelude::*;
use popkern_script::{Model, Plane, RevolveEnd};

fn z_axis() -> AxisDesc {
    AxisDesc::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    }
}

fn arc_radii(brep: &popkern_core::topology::BRep) -> Vec<f64> {
    let mut radii: Vec<f64> = brep
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { radius, .. })
            | Some(CurveGeometry::Circle { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    radii.sort_by(|a, b| a.total_cmp(b));
    radii
}

fn nurbs_edges(brep: &popkern_core::topology::BRep) -> usize {
    brep.edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Nurbs(_))
            )
        })
        .count()
}

fn has_surface(brep: &popkern_core::topology::BRep, kind: &str) -> bool {
    brep.faces.values().any(|f| {
        matches!(
            (kind, f.surface.as_ref().map(|s| &s.geometry)),
            ("sphere", Some(SurfaceGeometry::Sphere { .. }))
                | ("torus", Some(SurfaceGeometry::Torus { .. }))
        )
    })
}

/// A half-disc against the axis revolves to a sphere of radius 0.6.
fn sphere(m: &mut Model) -> Body {
    m.profile(Plane::XZ)
        .path(0.0, -0.6)
        .arc_to((0.6, 0.0), (0.0, 0.6))
        .line_to(0.0, -0.6)
        .close()
        .sketch("half_disc")
        .expect("sketch")
        .revolve(m, z_axis(), RevolveEnd::FullCircle, "ball")
        .expect("revolve")
}

#[test]
fn a_cap_cut_off_a_sphere_leaves_one_exact_circular_rim_that_selects() {
    let mut m = Model::new();
    let ball = sphere(&mut m);
    let brep = m.brep_of(ball).expect("brep");
    assert!(
        has_surface(&brep, "sphere"),
        "the revolve is an exact sphere"
    );
    // Slice the cap above z = 0.3 off with a box.
    let cutter = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.3))
                .rect(2.0, 2.0)
                .sketch("cap_box")?
                .extrude(m, End::Blind(1.0.into()))
        })
        .expect("cutter");
    let mut ball = ball.subtract(&mut m, cutter).expect("subtract");
    let brep = m.brep_of(ball).expect("brep");
    let radii = arc_radii(&brep);
    let rim = (0.6_f64.powi(2) - 0.3_f64.powi(2)).sqrt();
    assert!(
        radii.iter().any(|r| (r - rim).abs() < 1e-9),
        "the cap rim is an exact arc of radius {rim:.4}: {radii:?}"
    );
    assert_eq!(
        nurbs_edges(&brep),
        0,
        "no NURBS edge remains on the capped ball"
    );
    assert!(brep.validate_all().is_empty());
    // The rim selects as circular at its height and takes a fillet.
    m.label_edge("cap_rim", &edges().circular().at_height(0.3))
        .expect("the cap rim is circular");
    ball = ball
        .fillet(&mut m, edges().labeled("cap_rim"), 0.05)
        .expect("fillet on the cap rim");
    let _ = ball;
    assert_eq!(m.report().failures().count(), 0);
    m.assert_watertight().expect("watertight");
}

#[test]
fn a_torus_cut_above_its_spine_keeps_both_latitude_rims_and_its_seams_exact() {
    let mut m = Model::new();
    let ring = m
        .profile(Plane::XZ)
        .circle(0.3)
        .at(1.0, 0.0)
        .sketch("tube_section")
        .expect("sketch")
        .revolve(&mut m, z_axis(), RevolveEnd::FullCircle, "ring")
        .expect("revolve");
    let brep = m.brep_of(ring).expect("brep");
    assert!(has_surface(&brep, "torus"), "the revolve is an exact torus");
    // The cut sits above the torus's own equator seam (a cut THROUGH
    // that seam is the coincident-edge class, PK-420).
    let cutter = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.1))
                .rect(3.0, 3.0)
                .sketch("top_box")?
                .extrude(m, End::Blind(1.0.into()))
        })
        .expect("cutter");
    let cut = ring.subtract(&mut m, cutter).expect("subtract");
    let brep = m.brep_of(cut).expect("brep");
    let radii = arc_radii(&brep);
    let half_chord = (0.3_f64.powi(2) - 0.1_f64.powi(2)).sqrt();
    for want in [1.0 - half_chord, 1.0 + half_chord] {
        assert!(
            radii.iter().any(|r| (r - want).abs() < 1e-6),
            "latitude rim of radius {want:.4}: {radii:?}"
        );
    }
    assert_eq!(
        nurbs_edges(&brep),
        0,
        "no NURBS edge remains on the cut torus"
    );
    assert!(brep.validate_all().is_empty());
    // Both rims select as circular at the cut height; the cap is an
    // annulus between them.
    let annulus = m.label_face("cap", &faces().kind(SurfaceKind::Plane).at_height(0.1));
    assert!(annulus.is_ok(), "{annulus:?}");
    m.assert_watertight().expect("watertight");
}

/// The coarse end of the envelope: the sphere face's hole rings ride
/// the outer ring's period sheet, so the mesh stays closed and the
/// volume stays right when the model's tolerance is 1e-3.
#[test]
fn a_bored_ball_at_a_coarse_tolerance_meshes_watertight() {
    let mut m =
        Model::with_tolerance(popkern_core::tolerance::Tolerance::default().with_linear(1.0e-3));
    let ball = sphere(&mut m);
    let bore = m
        .tool(|m| {
            m.profile(Plane::xy_at(-1.0))
                .circle(0.2)
                .sketch("bore")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("bore");
    let bored = ball.subtract(&mut m, bore).expect("subtract");
    let brep = m.brep_of(bored).expect("brep");
    assert!(brep.validate_all().is_empty());
    // Both bore rims (r = 0.2) and the sphere's seam piece are exact
    // arcs; the bore's own straight seam stays NURBS — the recogniser
    // refuses a collinear run at this gate instead of answering a
    // circle of radius 8078 with a zero-degree sweep.
    let radii = arc_radii(&brep);
    assert_eq!(
        radii.iter().filter(|r| (*r - 0.2).abs() < 1e-3).count(),
        2,
        "{radii:?}"
    );
    assert!(radii.iter().any(|r| (r - 0.6).abs() < 1e-3), "{radii:?}");
    assert!(
        radii.iter().all(|r| *r < 1.0),
        "no runaway-radius arc from a straight seam: {radii:?}"
    );
    m.assert_watertight().expect("watertight at 1e-3");
}
