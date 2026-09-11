//! Frame-convention pins for the thread engine (P1 probe promoted):
//! with the axis-referenced auxiliary frame, profile +x maps to the
//! inward radial direction (exactly, at frame samples) and profile
//! +y to the binormal — along +axis for a right-handed helix, −axis
//! for a left-handed one, tilted off the axis by the lead angle (the
//! cross-section is perpendicular to the helix). The engine's
//! `X_PER_RADIAL` / `y_per_axial` constants encode this.
//!
//! Tolerances here are the real second-order lead-angle effects
//! (azimuthal binormal displacement folding into radius, binormal z
//! shortening), NOT construction slack: radius 1e-4, z 1e-3 at
//! radius 10 / lead 1.

use popkern_core::tolerance::Tolerance;
use popkern_sweep::{
    FrameMethod, HelixPitch, HelixSpec, SweepOptions, SweepProfile, helical_sweep,
};

use crate::engine::{axis_line, helical_tool, plan_turns, polygon_wire, sampled_centroid};
use crate::spec::ThreadSpec;
use crate::thread_form;

/// Sweep an asymmetric right triangle 1.25 turns and read the corner
/// vertices at the start station: their radii and heights measure
/// the frame mapping directly. The start station sits at azimuth
/// ~π/2 (the axis basis start); 1.25 turns puts the end station a
/// quarter turn away, so the azimuth filter isolates the start.
fn measure_start_corners(left_handed: bool) -> Vec<(f64, f64)> {
    let points = [[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]];
    let wire = polygon_wire(&points, "probe").expect("wire");
    let turns = 1.25;
    let helix = HelixSpec {
        axis_origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        axis_dir: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        start_radius: 10.0,
        end_radius: 10.0,
        pitch: HelixPitch::Constant(1.0),
        turns,
        start_angle: 0.0,
        left_handed,
    };
    let options = SweepOptions {
        frame: FrameMethod::Auxiliary(axis_line(0.0, turns)),
        ..SweepOptions::default()
    };
    let result = helical_sweep(
        &SweepProfile::Wire(wire),
        &helix,
        &options,
        Tolerance::default(),
        None,
    )
    .expect("helical sweep");
    let mut corners: Vec<(f64, f64)> = result
        .brep
        .vertices
        .iter()
        .map(|(_, v)| v.position)
        .filter(|p| (p.y.atan2(p.x) - std::f64::consts::FRAC_PI_2).abs() < 0.3)
        .map(|p| ((p.x * p.x + p.y * p.y).sqrt(), p.z))
        .collect();
    // Sort by (z, radius): radii of distinct corners can collide to
    // within the lead-angle effects, z values separate cleanly.
    corners.sort_by(|a, b| (a.1, a.0).partial_cmp(&(b.1, b.0)).expect("finite"));
    eprintln!("[frame-probe] lh={left_handed} corners={corners:?}");
    corners
}

fn assert_convention(left_handed: bool) {
    let corners = measure_start_corners(left_handed);
    assert_eq!(
        corners.len(),
        3,
        "three profile corners at the start station"
    );
    let points = [[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]];
    let wire = polygon_wire(&points, "probe").expect("wire");
    let (cx, cy) = sampled_centroid(&wire);
    // x → inward radial: radius = 10 − (x − cx).
    // y → binormal: z = ±(y − cy) (+ right-handed, − left-handed).
    let axial = if left_handed { -1.0 } else { 1.0 };
    let mut expected: Vec<(f64, f64)> = points
        .iter()
        .map(|p| (10.0 - (p[0] - cx), axial * (p[1] - cy)))
        .collect();
    expected.sort_by(|a, b| (a.1, a.0).partial_cmp(&(b.1, b.0)).expect("finite"));
    for (got, want) in corners.iter().zip(expected.iter()) {
        assert!(
            (got.0 - want.0).abs() < 1e-4 && (got.1 - want.1).abs() < 1e-3,
            "start-station corner {got:?} != expected {want:?} \
             (lh={left_handed}; frame convention drifted — update engine constants)"
        );
    }
}

#[test]
fn frame_convention_right_handed() {
    assert_convention(false);
}

#[test]
fn frame_convention_left_handed() {
    assert_convention(true);
}

/// The finished tool lands the authored ABSOLUTE radial coordinates:
/// vertex radii span [minor radius, major radius + embed] to within
/// the lead-angle second-order effects.
#[test]
fn helical_tool_lands_absolute_radii() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let profile = thread_form::external_groove(&spec).expect("groove profile");
    let (z0, turns) = plan_turns(0.0, 3.0, spec.lead());
    let tool = helical_tool(&spec, &profile, z0, turns, 0, "pin", Tolerance::default())
        .expect("tool body");
    let radii: Vec<f64> = tool
        .vertices
        .iter()
        .map(|(_, v)| {
            let p = v.position;
            (p.x * p.x + p.y * p.y).sqrt()
        })
        .collect();
    let r_min = radii.iter().cloned().fold(f64::MAX, f64::min);
    let r_max = radii.iter().cloned().fold(f64::MIN, f64::max);
    let want_min = spec.minor_diameter() / 2.0;
    let want_max = spec.major_diameter() / 2.0 + thread_form::embed_margin(&spec);
    assert!(
        (r_min - want_min).abs() < 1e-4,
        "tool inner radius {r_min} != minor radius {want_min}"
    );
    assert!(
        (r_max - want_max).abs() < 1e-4,
        "tool outer radius {r_max} != embed radius {want_max}"
    );
}

#[test]
fn plan_turns_quarter_integral_with_runout() {
    let (z0, turns) = plan_turns(0.0, 3.0, 1.0);
    assert!((z0 + 0.5).abs() < 1e-12);
    assert!((turns - 4.0).abs() < 1e-12, "3 span + 1 runout = 4 turns");
    let (_, turns) = plan_turns(0.0, 3.1, 1.0);
    assert!((turns - 4.25).abs() < 1e-12, "rounds UP to quarter turns");
    assert!((turns * 4.0).fract().abs() < 1e-12);
}

/// PK-475 pin: both helical thread tools are folded-OUTWARD solids
/// (flag-folded normals out of the material). The sweep normalizes
/// closed-profile winding, so the tool is outward regardless of how
/// the `thread_form` profile is authored. A folded-inward tool
/// enters the boolean with every fragment walking shared edges
/// reversed; the flush-crest shank wall's band cells then assemble
/// wrong-twin seam ladders whose loops gape one pitch (the
/// external_thread_on_existing_shank decline class). Winding-honesty
/// audits cannot see this — a coherently inside-out solid is honest
/// per face — so the material sense is pinned here directly.
#[test]
fn thread_tools_fold_outward() {
    use popkern_core::validation::MaterialSense;
    let spec = ThreadSpec::parse("M8x1").expect("M8x1");
    let (z0, turns) = plan_turns(6.0, 10.0, spec.lead());
    for (name, profile) in [
        (
            "external groove",
            thread_form::external_groove(&spec).expect("groove profile"),
        ),
        (
            "internal ridge",
            thread_form::internal_ridge(&spec, 0.0).expect("ridge profile"),
        ),
    ] {
        let tool = helical_tool(&spec, &profile, z0, turns, 0, name, Tolerance::default())
            .expect("tool body");
        let seed = tool.faces.iter().next().expect("tool has faces").0;
        let sense = tool
            .component_material_sense(seed)
            .expect("material sense measurable");
        assert_eq!(
            sense,
            MaterialSense::FoldedOutward,
            "{name} tool is folded-inward (inside-out solid)"
        );
    }
}
