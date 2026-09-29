//! PK-599: the KP20 setback cap served from the constrained bi-9
//! resampling keeps every certified claim and rim bound of the bi-6 cap
//! it replaces (A/B through `with_kp20_bi9_resample_disabled`), while
//! the sampling artefact leaves the interior (the fairness report,
//! `popkern-examples/tests/variational_corner_fairness.rs`).

use popkern_blend::dispatch_overrides::with_kp20_bi9_resample_disabled;
use popkern_blend::{
    BlendResult, ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::edge_seam_certificate::{certify_edge_seam, seam_certify_options};
use popkern_core::topology::{BRep, VertexHandle};

fn fillet(body: &BRep, vertex: VertexHandle, radii: &[f64]) -> BlendResult {
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    fillet_edges(body, &specs, &options, Tolerance::default(), None).expect("fillet")
}

/// Worst certified bounds `(g0, g1, g2)` over the cap's rim and spoke seams.
fn cap_bounds(result: &BlendResult) -> ([f64; 3], [f64; 3]) {
    let caps: std::collections::HashSet<_> = result.corner_patches.iter().map(|p| p.face).collect();
    let brep = &result.brep;
    let options = seam_certify_options(brep);
    let (mut rim, mut spoke) = ([0.0f64; 3], [0.0f64; 3]);
    for (edge, data) in brep.edges.iter() {
        let he = &brep.half_edges[data.half_edge];
        let (fa, fb) = (
            brep.loops[he.loop_].face,
            brep.loops[brep.half_edges[he.twin].loop_].face,
        );
        let (ca, cb) = (caps.contains(&fa), caps.contains(&fb));
        if !(ca || cb) || fa == fb {
            continue;
        }
        let c = certify_edge_seam(brep, edge, &options);
        assert!(
            c.bounds_cover_seam(),
            "seam {edge:?} not covered: {:?}",
            c.status
        );
        let slot = if ca && cb { &mut spoke } else { &mut rim };
        for (s, b) in slot.iter_mut().zip([c.g0_upper, c.g1_upper, c.g2_upper]) {
            *s = s.max(b);
        }
    }
    (rim, spoke)
}

fn assert_same_claims(name: &str, build: &dyn Fn() -> BlendResult) {
    let resampled = build();
    let bi6 = with_kp20_bi9_resample_disabled(build);
    let claims = |r: &BlendResult| {
        r.corner_patches
            .iter()
            .map(|p| p.achieved_continuity)
            .collect::<Vec<_>>()
    };
    assert_eq!(claims(&resampled), claims(&bi6), "{name}: claims");
    assert!(
        claims(&resampled).iter().all(|l| *l == ContinuityLevel::G2),
        "{name}: {:?}",
        claims(&resampled)
    );
    let ((rim_a, spoke_a), (rim_b, _)) = (cap_bounds(&resampled), cap_bounds(&bi6));
    // Rim seams are carried by the held rim bands: bounds identical.
    for k in 0..3 {
        assert!(
            (rim_a[k] - rim_b[k]).abs() <= 1e-9 + 1e-6 * rim_b[k],
            "{name}: rim {rim_a:?} vs {rim_b:?}"
        );
    }
    // Spokes stay certified G2 under the commercial gates.
    assert!(
        spoke_a[0] <= 1e-5 && spoke_a[1] <= 5e-2 && spoke_a[2] <= 5e-2,
        "{name}: spoke {spoke_a:?}"
    );
    resampled
        .brep
        .validate()
        .unwrap_or_else(|e| panic!("{name}: {e}"));
}

#[test]
fn bi9_cap_keeps_the_certified_claims_on_the_cube_corner() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9)
        .map(|(h, _)| h)
        .expect("corner");
    assert_same_claims("cube corner 0.2", &|| fillet(&cube, corner, &[0.2]));
}

#[test]
fn bi9_cap_keeps_the_certified_claims_on_the_variable_pentagon() {
    let body = BRep::n_gonal_pyramid(5, 1.0, 0.8, Tolerance::default()).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    assert_same_claims("pentagon5 variable", &|| {
        fillet(&body, apex, &[0.02, 0.04, 0.06, 0.04, 0.02])
    });
}

#[test]
fn bi9_cap_keeps_the_certified_claims_on_the_blunt_pentagon() {
    let body = BRep::n_gonal_pyramid(5, 1.0, 3.0, Tolerance::default()).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    let build = || fillet(&body, apex, &[0.05]);
    let (a, b) = (
        cap_bounds(&build()),
        cap_bounds(&with_kp20_bi9_resample_disabled(build)),
    );
    eprintln!(
        "blunt bi9 rim {:?} spoke {:?}; bi6 rim {:?} spoke {:?}",
        a.0, a.1, b.0, b.1
    );
    assert_same_claims("pentagon5 blunt", &build);
}
