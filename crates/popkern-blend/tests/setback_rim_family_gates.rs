//! PK-575: rim (cap↔strip) G2 across the setback apex family, measured
//! on the whole seam with exact jets (`seam_profile`, 201 samples,
//! 1e-3 end inset), plus each fixture's evidence verdict and its
//! setback-corner compatibility record.
//!
//! Before PK-575 (honest instrument): worst rim g2 1.6e-2 flagship,
//! 0.30 pentagon5_blunt, 7.7e-2 hexagon6_blunt, 0.70 hexagon6_mixed,
//! 6.0 pentagon5_variable; worst rim g1 up to 8e-5. The binding term
//! was the corner clusters' shared translation. It could not give two
//! incident sides different slopes, so it was frozen wherever they
//! disagreed. It is replaced by affine cluster fields where the fan
//! seam at the junction is a parametric subdivision seam
//! (`construct/setback_cluster_field.rs`), and by dual-jet corners plus
//! a shared junction translation at spoke junctions
//! (`promote_spoke_junctions`). Bars are 1.25× the 2026-09-27
//! measurement.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::dispatch_overrides::with_cap_fairness_ladder_disabled;
use popkern_blend::seam_profile::seam_profile;
use popkern_blend::{
    BlendResult, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};
use popkern_core::geometry::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

struct Fixture {
    name: &'static str,
    n: usize,
    height: f64,
    radii: Vec<f64>,
    /// Worst rim transverse-curvature bar away from the rim ends (form
    /// norm, 1/length).
    rim_g2_bar: f64,
    /// Worst rim transverse-curvature bar in the end layers.
    rim_end_bar: f64,
    /// Worst rim g1 bar (normal-angle sine).
    rim_g1_bar: f64,
    /// Evidence verdict the fixture must reach.
    achieved: ContinuityLevel,
}

fn family() -> Vec<Fixture> {
    vec![
        Fixture {
            name: "pentagon_flagship",
            n: 5,
            height: 1.2,
            radii: vec![0.16; 5],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[0],
            rim_end_bar: 1.25 * RIM_END_MEASURED[0],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[0]),
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "pentagon5_blunt",
            n: 5,
            height: 3.0,
            radii: vec![0.05; 5],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[1],
            rim_end_bar: 1.25 * RIM_END_MEASURED[1],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[1]),
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "hexagon6_blunt",
            n: 6,
            height: 3.0,
            radii: vec![0.05; 6],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[2],
            rim_end_bar: 1.25 * RIM_END_MEASURED[2],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[2]),
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "hexagon6_mixed",
            n: 6,
            height: 0.8,
            radii: vec![0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[3],
            rim_end_bar: 1.25 * RIM_END_MEASURED[3],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[3]),
            achieved: ContinuityLevel::G2,
        },
        // Rim 1.3e-3 (PK-581), g1 3.1e-8; the verdict is G1 because
        // three fan-internal cap-cap seams carry a certified G2 failure
        // (PK-580).
        Fixture {
            name: "pentagon5_variable",
            n: 5,
            height: 0.8,
            radii: vec![0.02, 0.04, 0.06, 0.04, 0.02],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[4],
            rim_end_bar: 1.25 * RIM_END_MEASURED[4],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[4]),
            // G2 since PK-582 (PK-580): its failing "fan-internal"
            // cap-cap seams were the spokes carrying the rim's junction
            // translation (certified curvature-jump LOWER bounds
            // 0.050-0.075 on edges 8/10/50/52). The spoke re-reconcile
            // certifies every seam under the 5e-2 gate.
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "hexagon6_variable",
            n: 6,
            height: 0.8,
            radii: vec![0.02, 0.04, 0.06, 0.06, 0.04, 0.02],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[5],
            rim_end_bar: 1.25 * RIM_END_MEASURED[5],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[5]),
            // G1 without PK-582's spoke re-reconcile, like
            // pentagon5_variable (PK-580 row 905).
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "octagon8_mixed",
            n: 8,
            height: 0.8,
            radii: vec![0.04, 0.06, 0.05, 0.07, 0.05, 0.06, 0.04, 0.06],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[6],
            rim_end_bar: 1.25 * RIM_END_MEASURED[6],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[6]),
            achieved: ContinuityLevel::G2,
        },
    ]
}

/// Worst rim form-norm curvature jump away from the rim ends, in
/// family order (measured 2026-10-11, PK-776). History: 1.545e-5 /
/// 1.700e-4 / 2.426e-4 / 1.212e-4 / 1.329e-3 / 2.401e-3 / 4.326e-4 at
/// PK-773. The curvature stations' exact row second derivative halved
/// the blunt fixtures; before it, a central difference clamped into
/// `[h, 1 − h]` read the wrong station at the row ends.
///
/// The rim is gated in two parts, interior and END LAYERS (the last ~2%
/// of each rim, a spoke-junction dual end; PK-581). The end layers used
/// to be roundoff-ill-conditioned: a one-ulp change of the pyramid
/// height moved them by up to 94%. That noise came from three
/// finite-difference support normals in the rim passes and from the
/// seam measurement's exact rational form of the dual-jet corner
/// pieces, ill-conditioned next to activated corners. Both are exact
/// now (PK-776), and `rim_measures_are_well_conditioned` pins the
/// move under 2%, so each constant is a single measurement.
const RIM_G2_MEASURED: [f64; 7] = [
    6.899e-6, 7.616e-5, 9.818e-5, 8.853e-5, 1.237e-3, 2.167e-3, 2.079e-4,
];

/// Worst rim form-norm curvature jump in the end layers (see above).
/// Was 4.267e-5 / 6.576e-4 / 1.079e-3 / 1.780e-4 / 1.287e-3 / 2.379e-3 /
/// 8.021e-4, each the larger of two perturbed measurements.
const RIM_END_MEASURED: [f64; 7] = [
    3.930e-5, 5.306e-4, 8.080e-4, 1.162e-4, 1.168e-3, 2.071e-3, 1.788e-4,
];

/// Worst rim g1, in family order, after the PK-584 band refit (was
/// 6.278e-9 / 1.464e-8 / 1.780e-8 / 8.946e-7 / 1.309e-5: the variable
/// fixtures' rim tangent no longer tilts off the strip).
const RIM_G1_MEASURED: [f64; 7] = [
    7.324e-9, 2.358e-8, 2.552e-8, 1.675e-8, 3.121e-8, 3.120e-8, 3.868e-8,
];

/// Rim G1 floor (normal-angle sine): every fixture is under it with at
/// least 3x margin since PK-584. It was 1e-5, with pentagon5_variable
/// (1.3e-5) excepted.
const RIM_G1_BAR: f64 = 1.0e-7;

fn build(fx: &Fixture) -> (BlendResult, BTreeSet<FaceHandle>) {
    build_at(fx, fx.height)
}

/// [`build`] with the pyramid height overridden.
fn build_at(fx: &Fixture, height: f64) -> (BlendResult, BTreeSet<FaceHandle>) {
    let tol = Tolerance::default();
    let pyr = BRep::n_gonal_pyramid(fx.n, 1.0, height, tol).expect("pyramid");
    let apex = pyr
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    let mut edges: Vec<EdgeHandle> = pyr
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(fx.radii.iter())
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let pre: BTreeSet<FaceHandle> = pyr.faces.keys().collect();
    // The family pins the KP20 fan's rims, so it runs with the
    // interior-fairness ladder off (PK-612): hexagon6_variable's fan
    // bends to 4.73 / r_min, past the bound, and in production the
    // fairer tangent fit (2.89) serves it at G1. Every other member's
    // fan is under the bound (octagon8_mixed 3.92, hexagon6_mixed 2.00).
    let result =
        with_cap_fairness_ladder_disabled(|| fillet_edges(&pyr, &specs, &options, tol, None))
            .expect("fillet");
    (result, pre)
}

/// Samples at each rim end read as the END LAYER (4 of 201, ~2% of
/// the rim at either end).
const END_LAYER: usize = 4;

/// Worst rim form jump and g1 over every cap↔non-support seam.
fn rim_worst(brep: &BRep, caps: &BTreeSet<FaceHandle>, pre: &BTreeSet<FaceHandle>) -> (f64, f64) {
    let (interior, end, g1) = rim_split(brep, caps, pre);
    (interior.max(end), g1)
}

/// Worst rim form jump in the interior and in the end layers, and g1.
fn rim_split(
    brep: &BRep,
    caps: &BTreeSet<FaceHandle>,
    pre: &BTreeSet<FaceHandle>,
) -> (f64, f64, f64) {
    let (mut g2, mut end, mut g1) = (0.0f64, 0.0f64, 0.0f64);
    for (eh, e) in brep.edges.iter() {
        let he = &brep.half_edges[e.half_edge];
        let fa = brep.loops[he.loop_].face;
        let fb = brep.loops[brep.half_edges[he.twin].loop_].face;
        let is_strip = |f: FaceHandle| !caps.contains(&f) && !pre.contains(&f);
        let rim = (caps.contains(&fa) && is_strip(fb)) || (caps.contains(&fb) && is_strip(fa));
        if !rim {
            continue;
        }
        let profile = seam_profile(brep, eh, 201, 1.0e-3).expect("rim seam profile");
        let n = profile.len();
        for (i, s) in profile.iter().enumerate() {
            let jump = s.form_jump.or(s.sampled_jump()).unwrap_or(0.0);
            if i < END_LAYER || i + END_LAYER >= n {
                end = end.max(jump);
            } else {
                g2 = g2.max(jump);
            }
            g1 = g1.max(s.g1.unwrap_or(0.0));
        }
    }
    (g2, end, g1)
}

#[test]
fn setback_family_rims_meet_g2_and_every_corner_is_compatible() {
    let gates = ContinuityGates::commercial();
    for fx in family() {
        let (result, pre) = build(&fx);
        let caps: BTreeSet<FaceHandle> = result.corner_patches.iter().map(|p| p.face).collect();
        let (g2, end, g1) = rim_split(&result.brep, &caps, &pre);
        let blend: Vec<FaceHandle> = result
            .brep
            .faces
            .keys()
            .filter(|f| !pre.contains(f))
            .collect();
        let evidence = measure_blend_continuity_evidence(
            &result.brep,
            &blend,
            ContinuityLevel::G2,
            &gates,
            33,
        );
        let corners = &result.repair_report.setback_corners;
        let worst_corner = corners.iter().map(|c| c.form_jump).fold(0.0, f64::max);
        eprintln!(
            "{}: rim g2 {g2:.3e} (bar {:.3e}) end {end:.3e} (bar {:.3e}) g1 {g1:.3e}; achieved {:?}; {} corners, worst |ΔII| {worst_corner:.3e}",
            fx.name,
            fx.rim_g2_bar,
            fx.rim_end_bar,
            evidence.achieved,
            corners.len()
        );
        assert!(g2 <= fx.rim_g2_bar, "{}: rim g2 {g2:.3e}", fx.name);
        assert!(
            end <= fx.rim_end_bar,
            "{}: rim end layer {end:.3e}",
            fx.name
        );
        assert!(g1 <= fx.rim_g1_bar, "{}: rim g1 {g1:.3e}", fx.name);
        assert_eq!(evidence.achieved, fx.achieved, "{}", fx.name);
        // Row 886: one record per setback corner (n at an n-gon apex),
        // each second-order compatible (strips flat at the corner).
        assert!(
            corners.len() >= fx.n,
            "{}: {} corners",
            fx.name,
            corners.len()
        );
        for c in corners {
            assert_eq!(c.compatible_order, 2, "{}: {c:?}", fx.name);
            assert!(c.form_jump < 1e-6, "{}: {c:?}", fx.name);
        }
    }
}

/// PK-776: the rim measures are a well-conditioned function of the
/// input. A one-ulp change of the pyramid height moves the interior and
/// end-layer rim G2 by under 2% on the two fixtures that moved most
/// before (octagon8_mixed's end layers by 94%, pentagon5_blunt's by 24%;
/// both now under 0.7%). Negative controls, each removing one mechanism
/// on its own: a finite-difference support normal in the rim jet (40%),
/// the curvature stations (25%) or the band refit (46%), and the seam
/// measurement's exact rational form for dual-jet patches (38%).
#[test]
fn rim_measures_are_well_conditioned() {
    for name in ["octagon8_mixed", "pentagon5_blunt"] {
        let fx = family()
            .into_iter()
            .find(|f| f.name == name)
            .expect("family member");
        let measure = |height: f64| {
            let (result, pre) = build_at(&fx, height);
            let caps: BTreeSet<FaceHandle> = result.corner_patches.iter().map(|p| p.face).collect();
            let (g2, end, _) = rim_split(&result.brep, &caps, &pre);
            (g2, end)
        };
        let (g2, end) = measure(fx.height);
        let (g2_up, end_up) = measure(fx.height.next_up());
        let moved = |a: f64, b: f64| (a - b).abs() / a.max(b);
        eprintln!(
            "{name}: rim g2 {g2:.4e} -> {g2_up:.4e}, end {end:.4e} -> {end_up:.4e} under one ulp of height"
        );
        assert!(
            moved(g2, g2_up) < 0.02,
            "{name}: rim g2 {g2:.4e} vs {g2_up:.4e}"
        );
        assert!(
            moved(end, end_up) < 0.02,
            "{name}: end layer {end:.4e} vs {end_up:.4e}"
        );
    }
}

/// Row 887's negative control: a cap whose depth-2 rows are pushed off
/// the reconciled solution keeps G0 and G1 but must read over the
/// flagship's rim bar. A gate that cannot fail is not evidence.
#[test]
fn a_seeded_depth2_error_reads_over_the_rim_bar() {
    let fx = family().into_iter().next().expect("flagship");
    let (mut result, pre) = build(&fx);
    let caps: BTreeSet<FaceHandle> = result.corner_patches.iter().map(|p| p.face).collect();
    let (clean, _) = rim_worst(&result.brep, &caps, &pre);
    for &face in &caps {
        let Some(sof) = result.brep.faces[face].surface.as_mut() else {
            continue;
        };
        let SurfaceGeometry::Gregory(patch) = &mut sof.geometry else {
            continue;
        };
        let Some(grid) = patch.binonic_grid_mut() else {
            continue;
        };
        // Depth-2 rows of all four sides, off the corners: G0 and G1
        // along every side are untouched.
        let n = {
            let a = grid[0][0];
            let b = grid[9][0];
            let c = grid[0][9];
            (b - a).cross(&(c - a)).normalized().expect("patch normal")
        };
        for k in 3..=6usize {
            for (i, j) in [(2, k), (7, k), (k, 2), (k, 7)] {
                let p = grid[i][j];
                grid[i][j] = popkern_core::geometry::Point3::new(
                    p.x + 1e-5 * n.x,
                    p.y + 1e-5 * n.y,
                    p.z + 1e-5 * n.z,
                );
            }
        }
    }
    let (seeded, _, g1) = rim_split(&result.brep, &caps, &pre);
    eprintln!("clean {clean:.3e} seeded interior {seeded:.3e} (g1 {g1:.3e})");
    assert!(seeded > 10.0 * fx.rim_g2_bar, "seeded {seeded:.3e}");
    assert!(g1 <= RIM_G1_BAR, "the seed must stay G1: {g1:.3e}");
}
