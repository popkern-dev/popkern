//! Per-direction adaptive lattice density for NURBS charts
//! (extracted from `boundary_first.rs` under the structural line
//! budget). See [`directional_stations`].

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::TessellationOptions;

/// The interior lattice's station counts `(n_u, n_v)` given the grid
/// count `n` from the chart-diagonal probe. NURBS charts use
/// [`directional_stations`]: the diagonal probe reads a coiled chart
/// by phase luck. A torus takes the larger of the two per direction: a
/// band spanning a full minor turn (an imported pipe elbow) reads the
/// diagonal by the same luck (sagging 5x the chord height), while the
/// grid stays the floor native tori were tuned on (PK-641). Other
/// analytic charts keep the grid.
pub(super) fn lattice_counts(
    surface: &SurfaceGeometry,
    u_span: (f64, f64),
    v_span: (f64, f64),
    n: usize,
    options: &TessellationOptions,
) -> (usize, usize) {
    match surface {
        SurfaceGeometry::Nurbs(_) => directional_stations(surface, u_span, v_span, options),
        SurfaceGeometry::Torus { .. } => {
            let (n_u, n_v) = directional_stations(surface, u_span, v_span, options);
            (n_u.max(n), n_v.max(n))
        }
        _ => (n, n),
    }
}

/// Interior-lattice station counts per chart direction, refined until
/// the worst iso-line chord sagitta meets `chord_height` and (with an
/// angular facet bound) the worst iso-line tangent turn meets
/// `max_facet_angle`.
///
/// The single chart-diagonal probe that sizes the square grid reads a
/// COILED chart by phase luck: on a 6.5-turn helical tube wall the
/// diagonal's midpoint happened to land near its own chord, the grid
/// stayed at n=4 (metric cell 2.8 over a 0.25-wide band), every
/// interior point died to clearance, and the boundary-only Delaunay
/// spanned the tube with diameter chords — the PK-441 ribbon whose
/// fidelity decline sent the wall to the uniform-grid path (9% volume
/// deficit, rim T-junction cracks on cut walls). Sagitta is measured
/// along three cross-stations per direction so a chart curved in only
/// part of its span still refines; counts cap at 512 per direction
/// (the ruled path's station cap) with a combined budget guard.
pub(super) fn directional_stations(
    surface: &SurfaceGeometry,
    u_span: (f64, f64),
    v_span: (f64, f64),
    options: &TessellationOptions,
) -> (usize, usize) {
    let (u0, u1) = u_span;
    let (v0, v1) = v_span;
    let cap = 512usize;
    let max_turn = options.angular_bound();
    let refine = |dir_u: bool| -> usize {
        let mut n = options.min_subdivisions.max(2);
        loop {
            let mut worst = 0.0_f64;
            for cross in [0.25, 0.5, 0.75] {
                let at = |f: f64| -> Point3 {
                    if dir_u {
                        surface.point_at(u0 + (u1 - u0) * f, v0 + (v1 - v0) * cross)
                    } else {
                        surface.point_at(u0 + (u1 - u0) * cross, v0 + (v1 - v0) * f)
                    }
                };
                let mut prev = at(0.0);
                for k in 0..n {
                    let next = at((k + 1) as f64 / n as f64);
                    let mid = at((k as f64 + 0.5) / n as f64);
                    let chord_mid = Point3::new(
                        0.5 * (prev.x + next.x),
                        0.5 * (prev.y + next.y),
                        0.5 * (prev.z + next.z),
                    );
                    worst = worst.max(mid.distance_to(&chord_mid));
                    if let Some(max_turn) = max_turn {
                        // Iso-line tangent turn in chord-height units
                        // so one threshold decides both criteria.
                        let turn = crate::facet_angle::curve_turn(prev, mid, next);
                        worst = worst.max(options.chord_height * turn / max_turn);
                    }
                    prev = next;
                }
            }
            if worst <= options.chord_height || n >= cap {
                return n.min(cap);
            }
            n *= 2;
        }
    };
    let (mut n_u, mut n_v) = (refine(true), refine(false));
    // Combined budget: a chart that maxes BOTH directions would ask
    // for 512×512 interior points; halve the larger until the lattice
    // stays within a sane point budget.
    while n_u * n_v > 65_536 {
        if n_u >= n_v {
            n_u /= 2;
        } else {
            n_v /= 2;
        }
    }
    (n_u.max(2), n_v.max(2))
}
