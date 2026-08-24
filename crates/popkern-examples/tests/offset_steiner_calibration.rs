//! PK-323 row 468a: Steiner/tube volume bands on convex analytic
//! bodies (plus the torus below its reach), computed KERNEL-FREE.
//!
//! For each body the Minkowski-dilation volume is integrated on a
//! grid from the tessellation-based SDF and compared to the closed
//! form `V + A·d + M·d² + (4π/3)(χ/2)·d³`. Agreement calibrates BOTH
//! sides at once: the mean-curvature constants (derived by hand — a
//! wrong `M` shows up here immediately) and the SDF machinery every
//! offset-serving gate relies on. No offset needs to SERVE for this
//! battery: it is the acceptance substrate later stages measure
//! against.

use std::f64::consts::PI;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{filleted_cube, torus_lens};
use popkern_mesh::offset_oracle::{grid_offset_volume, steiner_offset_volume};
use popkern_mesh::point_query::PointMeshIndex;
use popkern_mesh::{TessellationOptions, tessellate};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn index_of(brep: &BRep) -> PointMeshIndex {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    PointMeshIndex::build(&mesh)
}

struct SteinerCase {
    name: &'static str,
    brep: BRep,
    volume: f64,
    area: f64,
    /// `∫H dA + Σ (θₑ/2)·ℓₑ`.
    mean_curvature: f64,
    euler: f64,
    distance: f64,
    band: f64,
}

#[test]
fn grid_offset_volume_matches_steiner_on_convex_analytic_bodies() {
    let d = 0.25;
    let fc = filleted_cube(tol()).expect("filleted cube");
    let fc_steiner = fc.steiner.expect("filleted cube constants");
    let tl = torus_lens(tol()).expect("torus lens");
    let tl_steiner = tl.steiner.expect("torus constants");
    let cases = vec![
        SteinerCase {
            name: "unit cube",
            brep: BRep::unit_cube(tol()).unwrap(),
            volume: 1.0,
            area: 6.0,
            // 12 edges, exterior angle π/2: Σ (θ/2)·ℓ = 3π.
            mean_curvature: 3.0 * PI,
            euler: 2.0,
            distance: d,
            band: 0.01,
        },
        SteinerCase {
            name: "cylinder r=1 h=2",
            brep: BRep::cylinder(1.0, 2.0, tol()).unwrap(),
            volume: 2.0 * PI,
            area: 6.0 * PI,
            // Lateral ∫H = πh; two rim circles at θ = π/2: π²r.
            mean_curvature: 2.0 * PI + PI * PI,
            euler: 2.0,
            distance: d,
            band: 0.01,
        },
        SteinerCase {
            name: "cone r=1 h=2",
            brep: BRep::cone(1.0, 2.0, tol()).unwrap(),
            volume: 2.0 * PI / 3.0,
            area: PI + PI * 5.0_f64.sqrt(),
            // Lateral ∫H = πh (κ₂ = cos α / ρ integrates to πℓcos α);
            // base rim at exterior angle π − atan(h/r); apex
            // contributes only to the cubic (Gauss) term.
            mean_curvature: 2.0 * PI + PI * (PI - 2.0_f64.atan()),
            euler: 2.0,
            distance: d,
            band: 0.01,
        },
        SteinerCase {
            name: "sphere r=1",
            brep: BRep::sphere(1.0, tol()).unwrap(),
            volume: 4.0 * PI / 3.0,
            area: 4.0 * PI,
            mean_curvature: 4.0 * PI,
            euler: 2.0,
            distance: d,
            band: 0.01,
        },
        SteinerCase {
            name: "filleted cube (exact rounded cube)",
            brep: fc.brep,
            volume: fc_steiner.volume,
            area: fc_steiner.surface_area,
            mean_curvature: fc_steiner.mean_curvature_integral,
            euler: fc_steiner.euler_characteristic,
            distance: d,
            band: 0.01,
        },
        SteinerCase {
            name: "torus lens (tube formula, below reach)",
            brep: tl.brep,
            volume: tl_steiner.volume,
            area: tl_steiner.surface_area,
            mean_curvature: tl_steiner.mean_curvature_integral,
            euler: tl_steiner.euler_characteristic,
            // Hole radius is 0.15: stay below the reach.
            distance: 0.10,
            band: 0.015,
        },
    ];
    for c in cases {
        let index = index_of(&c.brep);
        let grid = grid_offset_volume(&index, c.distance, 72);
        let steiner =
            steiner_offset_volume(c.volume, c.area, c.mean_curvature, c.euler, c.distance);
        let rel = (grid - steiner).abs() / steiner;
        assert!(
            rel <= c.band,
            "{}: grid dilation volume {grid:.6} vs Steiner {steiner:.6} (rel {rel:.5}, band {})",
            c.name,
            c.band
        );
    }
}
