//! PK-503 attribution gates: on the setback-cap fixtures, which zebra
//! stripe defects are the SURFACE's and which the DISPLAY MESH's, as
//! measured by the deferred zebra frame's two views
//! (`zebra_attribution` bin prints the full table).
//!
//! Measured 2026-09-27 at 640 px, 16 stripes (worst-seam medians;
//! normal jump in radians, slope jump per pixel):
//!
//! | fixture | seam | surf nj | mesh nj | surf slope | mesh slope |
//! |---|---|---|---|---|---|
//! | pentagon flagship | strip-support | 1.2e-5 | 5.3e-4 | 1.7e-4 | 2.0e-3 |
//! | pentagon flagship | cap-rim | 1.5e-3 | 6.9e-3 | 9.0e-3 | 1.8e-2 |
//! | pentagon flagship | cap-cap | 3.5e-4 | 1.9e-3 | 1.1e-3 | 2.8e-3 |
//! | pentagon5_blunt G2 | cap-cap | 3.4e-4 | 3.1e-3 | 4.3e-4 | 8.0e-3 |
//! | pentagon5_blunt G2 | cap-rim | 3.2e-2 | 2.4e-2 | 1.6e-1 | 1.4e-1 |
//! | cube corner | cap-rim | 2.7e-4 | 4.0e-3 | 1.3e-3 | 7.2e-3 |
//! | cube corner | cap-cap | 7.8e-4 | 3.9e-3 | 3.7e-3 | 1.6e-2 |
//! | pentagon5_variable G2 | cap-rim | 5.0e-3 | 7.7e-3 | 1.9e-2 | 1.6e-2 |
//! | pentagon5_variable G2 | cap-cap | 2.6e-3 | 7.3e-3 | 7.1e-3 | 1.5e-2 |
//!
//! Attribution: every strip-support kink the display shows is the
//! tessellation's (the surface view is G1/G2-clean there), and so is
//! most of every cap-cap kink; the cap-rim kink on the flagship is
//! half surface (the rim G2 residual PK-575 drives down) and half
//! tessellation; pentagon5_blunt's cap-rim kink is the surface's
//! (PK-575).
//!
//! The surface-view bars below are ratchets at 1.25x the measurement:
//! a construction change that roughens a cap reads over them. PK-575
//! tightens the cap-rim rows as the residual falls.

use popkern_examples::zebra_attribution::{Attribution, attribute, attribution_fixtures};

const SIZE: u32 = 640;
const STRIPES: u32 = 16;

fn attributions() -> Vec<(&'static str, Attribution)> {
    attribution_fixtures()
        .into_iter()
        .map(|fx| (fx.name, attribute(&fx, SIZE, STRIPES)))
        .collect()
}

#[test]
fn cap_fixtures_attribute_stripe_defects_to_surface_or_tessellation() {
    let all = attributions();
    let get = |name: &str| &all.iter().find(|(n, _)| *n == name).expect(name).1;
    for (name, at) in &all {
        let d = &at.disagreement;
        eprintln!(
            "{name}: mode max {:.3}° p99 {:.3}°, ray hits {}/{}",
            d.max_angle_rad.to_degrees(),
            d.p99_angle_rad.to_degrees(),
            d.ray_hit_pixels,
            d.evaluated_pixels
        );
        for (class, s, m) in &at.classes {
            eprintln!(
                "  {:<14} seams {:>2}: nj {:.2e}/{:.2e} slope {:.2e}/{:.2e} (surface/mesh)",
                class.label(),
                s.seams,
                s.normal_p50,
                m.normal_p50,
                s.slope_p50,
                m.slope_p50
            );
        }
        // Every pixel's surface point is located on the surface (the
        // ray hit, or a seam pixel's domain-boundary point).
        assert!(
            d.ray_hit_pixels as f64 >= 0.999 * d.evaluated_pixels as f64,
            "{name}: {} of {} pixels located",
            d.ray_hit_pixels,
            d.evaluated_pixels
        );
    }

    // Strip-support seams are exact G1 (cylinder/G2-strip on a plane):
    // the surface view reads them clean and every display kink there
    // is the tessellation's.
    for name in ["pentagon_flagship", "pentagon5_blunt_g2", "cube_corner"] {
        let (s, m) = get(name).class("strip-support").expect("strip-support");
        assert!(
            s.normal_p50 < 1e-4,
            "{name}: surface nj {:.2e}",
            s.normal_p50
        );
        assert!(
            m.slope_p50 > 4.0 * s.slope_p50,
            "{name}: strip-support kink must be the mesh's: {s:?} vs {m:?}"
        );
    }

    // Surface-view ratchets (1.25x the measurement). PK-601
    // (2026-09-29): the curvature-variation fair interior lowered every
    // cap reading, so all bars were re-measured then. The cube's
    // cap-rim slope fell from 2.30e-3 (PK-599's bi-9 guide-border
    // relaxation) to 2.12e-4, and the variable pentagon's cap-rim nj
    // from 1.25e-2 to 4.34e-3, with certified rim bounds unchanged. The
    // blunt pentagon's cap-rim reading (4.97e-3) is still the
    // instrument's curvature x pixel-footprint floor on its r = 0.05
    // strips, not its certified rim g1 (2.4e-8).
    let bars: &[(&str, &str, f64, f64)] = &[
        ("pentagon_flagship", "cap-rim", 7.17e-4, 5.18e-3),
        ("pentagon_flagship", "cap-cap", 1.90e-4, 1.93e-4),
        ("pentagon5_blunt_g2", "cap-rim", 4.97e-3, 3.30e-2),
        ("pentagon5_blunt_g2", "cap-cap", 2.35e-4, 5.78e-4),
        ("cube_corner", "cap-rim", 1.48e-4, 2.12e-4),
        ("pentagon5_variable_g2", "cap-rim", 4.34e-3, 3.71e-3),
        ("pentagon5_variable_g2", "cap-cap", 5.69e-4, 6.84e-4),
        ("cube_corner", "cap-cap", 7.57e-5, 2.98e-4),
    ];
    for &(name, class, nj, slope) in bars {
        let (s, _) = get(name).class(class).expect(class);
        assert!(
            s.normal_p50 <= 1.25 * nj && s.slope_p50 <= 1.25 * slope,
            "{name} {class}: surface nj {:.3e} (bar {:.3e}), slope {:.3e} (bar {:.3e})",
            s.normal_p50,
            1.25 * nj,
            s.slope_p50,
            1.25 * slope
        );
    }

    // Cap-cap seams are G1-exact (analytic g1 <= 2e-5 on all three,
    // `setback_rim_family_profile`): the surface view must read them
    // clean, and the display mesh owns most of their visible kink.
    for name in ["pentagon_flagship", "pentagon5_blunt_g2", "cube_corner"] {
        let (s, m) = get(name).class("cap-cap").expect("cap-cap");
        assert!(
            s.normal_p50 < 1e-3,
            "{name}: cap-cap surface nj {:.2e}",
            s.normal_p50
        );
        assert!(m.slope_p50 > 2.0 * s.slope_p50, "{name}: {s:?} vs {m:?}");
    }
}
