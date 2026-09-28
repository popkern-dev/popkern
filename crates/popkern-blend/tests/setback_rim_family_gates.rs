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
    /// Worst rim transverse-curvature bar (form norm, 1/length).
    rim_g2_bar: f64,
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
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[0]),
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "pentagon5_blunt",
            n: 5,
            height: 3.0,
            radii: vec![0.05; 5],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[1],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[1]),
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "hexagon6_blunt",
            n: 6,
            height: 3.0,
            radii: vec![0.05; 6],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[2],
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[2]),
            achieved: ContinuityLevel::G2,
        },
        Fixture {
            name: "hexagon6_mixed",
            n: 6,
            height: 0.8,
            radii: vec![0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
            rim_g2_bar: 1.25 * RIM_G2_MEASURED[3],
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
            rim_g1_bar: RIM_G1_BAR.max(1.25 * RIM_G1_MEASURED[4]),
            achieved: ContinuityLevel::G1,
        },
    ]
}

/// Worst rim form-norm curvature jump, in family order, measured
/// 2026-09-27 after the PK-584 band refit and the PK-581 curvature-unit
/// weighting of the rim curvature pass, against the 4.088e-5 / 5.978e-4
/// / 9.026e-4 / 1.825e-4 / 2.281e-3 before them. pentagon5_variable fell
/// 2.28e-3 -> 1.32e-3. Its remainder is an end layer ~0.02 from a
/// spoke-junction dual end, where the strip's transverse curvature ramps
/// to 0 and the junction point is pinned by the spoke boundary cell.
/// pentagon5_blunt stays +13% from the band refit: the setback-corner
/// k = 2 cell is the neighbour rim's fixed cross-curvature end cell.
const RIM_G2_MEASURED: [f64; 5] = [4.054e-5, 6.753e-4, 8.488e-4, 1.817e-4, 1.323e-3];

/// Worst rim g1, in family order, after the PK-584 band refit (was
/// 6.278e-9 / 1.464e-8 / 1.780e-8 / 8.946e-7 / 1.309e-5: the variable
/// fixtures' rim tangent no longer tilts off the strip).
const RIM_G1_MEASURED: [f64; 5] = [7.324e-9, 2.358e-8, 2.552e-8, 1.675e-8, 3.121e-8];

/// Rim G1 floor (normal-angle sine): every fixture is under it with at
/// least 3x margin since PK-584. It was 1e-5, with pentagon5_variable
/// (1.3e-5) excepted.
const RIM_G1_BAR: f64 = 1.0e-7;

fn build(fx: &Fixture) -> (BlendResult, BTreeSet<FaceHandle>) {
    let tol = Tolerance::default();
    let pyr = BRep::n_gonal_pyramid(fx.n, 1.0, fx.height, tol).expect("pyramid");
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
    let result = fillet_edges(&pyr, &specs, &options, tol, None).expect("fillet");
    (result, pre)
}

/// Worst rim form jump and g1 over every cap↔non-support seam.
fn rim_worst(brep: &BRep, caps: &BTreeSet<FaceHandle>, pre: &BTreeSet<FaceHandle>) -> (f64, f64) {
    let (mut g2, mut g1) = (0.0f64, 0.0f64);
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
        for s in &profile {
            g2 = g2.max(s.form_jump.or(s.sampled_jump()).unwrap_or(0.0));
            g1 = g1.max(s.g1.unwrap_or(0.0));
        }
    }
    (g2, g1)
}

#[test]
fn setback_family_rims_meet_g2_and_every_corner_is_compatible() {
    let gates = ContinuityGates::commercial();
    for fx in family() {
        let (result, pre) = build(&fx);
        let caps: BTreeSet<FaceHandle> = result.corner_patches.iter().map(|p| p.face).collect();
        let (g2, g1) = rim_worst(&result.brep, &caps, &pre);
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
            "{}: rim g2 {g2:.3e} (bar {:.3e}) g1 {g1:.3e}; achieved {:?}; {} corners, worst |ΔII| {worst_corner:.3e}",
            fx.name,
            fx.rim_g2_bar,
            evidence.achieved,
            corners.len()
        );
        assert!(g2 <= fx.rim_g2_bar, "{}: rim g2 {g2:.3e}", fx.name);
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
    let (seeded, g1) = rim_worst(&result.brep, &caps, &pre);
    eprintln!("clean {clean:.3e} seeded {seeded:.3e} (g1 {g1:.3e})");
    assert!(seeded > 10.0 * fx.rim_g2_bar, "seeded {seeded:.3e}");
    assert!(g1 <= RIM_G1_BAR, "the seed must stay G1: {g1:.3e}");
}
