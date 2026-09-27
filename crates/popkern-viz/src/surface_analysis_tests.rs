//! Tests for the analysis renderer: palette semantics, robust
//! ranging, and that zebra/curvature views of analytic fixtures
//! produce the expected image structure.

use super::*;
use popkern_core::tolerance::Tolerance;

#[test]
fn diverging_palette_endpoints_and_center() {
    let half = 1.0;
    let neg = diverging_color(-1.0, half);
    let zero = diverging_color(0.0, half);
    let pos = diverging_color(1.0, half);
    assert!(neg[2] > neg[0], "negative maps blue");
    assert!(pos[0] > pos[2], "positive maps red");
    assert!(
        zero[0] > 200 && zero[1] > 200 && zero[2] > 200,
        "zero maps near-white"
    );
}

#[test]
fn robust_range_ignores_spikes_and_flat_bulk() {
    // A part that is mostly flat with a curved band at |κ| = 1 and
    // one degenerate spike: the range must key on the band.
    let mut values = vec![0.0_f64; 5000];
    values.extend(std::iter::repeat_n(1.0, 800));
    values.push(1.0e9);
    let half = robust_half_range(&values);
    assert!(
        (half - 1.25).abs() < 1e-9,
        "nonzero median must key the range on the curved band, got {half}"
    );
}

#[test]
fn sphere_zebra_has_both_bands_and_curvature_is_uniform_red() {
    let tol = Tolerance::default();
    let brep = popkern_core::topology::BRep::sphere(0.8, tol).expect("sphere");
    let mesh =
        popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default()).expect("m");
    let camera = Camera::default();
    let zebra = render_surface_analysis(
        &brep,
        &mesh,
        &camera,
        160,
        120,
        SurfaceAnalysisMode::Zebra {
            stripe_count: 8,
            normals: ZebraNormals::Mesh,
        },
    )
    .expect("zebra");
    let (mut dark, mut light) = (0usize, 0usize);
    for px in zebra.chunks(4) {
        match px[0] {
            24 => dark += 1,
            236 => light += 1,
            _ => {}
        }
    }
    assert!(
        dark > 100 && light > 100,
        "a curved body must show both stripe bands (dark {dark}, light {light})"
    );

    // Positive-K body: the Gaussian map must be warm (red-dominant)
    // on every rendered surface pixel.
    let curv = render_surface_analysis(
        &brep,
        &mesh,
        &camera,
        160,
        120,
        SurfaceAnalysisMode::Curvature {
            kind: popkern_measure::CurvatureKind::Gaussian,
        },
    )
    .expect("curvature");
    let mut cold_pixels = 0usize;
    let mut surface_pixels = 0usize;
    for px in curv.chunks(4) {
        if px == [38, 38, 38, 255] {
            continue; // background
        }
        surface_pixels += 1;
        if px[2] > px[0].saturating_add(10) {
            cold_pixels += 1;
        }
    }
    assert!(surface_pixels > 500, "sphere must cover pixels");
    assert!(
        cold_pixels == 0,
        "a sphere has no negative Gaussian curvature ({cold_pixels}/{surface_pixels} cold)"
    );
}
