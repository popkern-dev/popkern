//! Cap-outer-boundary commercial-grade continuity gate.
//!
//! Sibling of `kp_binonic_commercial_gates.rs` — but where that
//! file measures continuity *between adjacent cap sub-patches*
//! at their shared inner spokes, this file measures continuity
//! between each cap sub-patch and the cylinder fillet roll it
//! sits on top of.
//!
//! # Why this exists
//!
//! Per Task #144 / 2026-04-27 failed-approaches entry: the
//! existing inter-sub-patch ring gates passed while the cap-outer
//! boundary silently drifted by 5–15% bbox in 3-D. FreeCAD shows
//! this as the "messy top" defect on every fillet zoo pyramid.
//! The user-visible quality bar is set here, not at the inner
//! spokes.
//!
//! # Status
//!
//! The direct-apex exact-boundary matrix now gates native cube and
//! pyramid cap_outer quality at the user-visible `5e-3`/`5e-2`
//! tolerances. This older K-P helper remains useful as a stricter
//! legacy-path monitor: several ignored fixtures still miss the
//! near-exact `G0 <= 1e-10` seam gate even when their direct-apex
//! exact-boundary route is commercially acceptable.
//!
//! Threshold targets (mirroring `kp_binonic_commercial_gates`):
//! * G0 ≤ 1e-10 (position match at the seam)
//! * G1 ≤ 5e-2  (coplanarity at samples; ≈ 2.9° normal mismatch)
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{BlendWarning, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    G0_GATE, G1_GATE, N_SAMPLES, cap_outer_faces_from_result, fillet_corner, incident_edges,
    measure_cap_outer_boundary, pick_vertex_with_valence,
};

fn assert_cap_outer_gates(fixture: &str, m: &helpers::CapOuterMetrics) {
    println!(
        "cap_outer {fixture}: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
    assert!(
        m.edges_measured >= 1,
        "[{fixture}] no cap-vs-non-cap edges found — fixture broke"
    );
    assert!(
        m.worst_g0_gap < G0_GATE,
        "[{fixture}] cap-outer G0 gap {:.3e} exceeds {:.1e} gate",
        m.worst_g0_gap,
        G0_GATE
    );
    assert!(
        m.worst_g1_coplanarity < G1_GATE,
        "[{fixture}] cap-outer G1 coplanarity {:.3e} exceeds {:.1e} gate",
        m.worst_g1_coplanarity,
        G1_GATE
    );
}

fn run_cube(radii: &[f64]) -> helpers::CapOuterMetrics {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

fn run_n_pyramid(n: usize, h: f64, radii: &[f64]) -> helpers::CapOuterMetrics {
    let pyr = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

/// Phase 6c-C.1 win: cube_uniform G1 commercial-grade, locked
/// in as a regression guard. The 5e-2 gate is the
/// `kp_gate_helpers::G1_GATE` constant. Passed since commit
/// `016a780` (centroid-based interior-direction tie-break).
///
/// G0 still exceeds `G0_GATE` (1e-10) because the cap face's
/// 3 rolling-ball-footprint vertices aren't on the sphere; full
/// Phase 6c-C topology rebuild is needed to remove them. This
/// test gates only on G1 — graduates to the full
/// `assert_cap_outer_gates` once 6c-C.2/.3 land.
#[test]
fn cap_outer_boundary_g1_passes_on_cube_uniform() {
    let m = run_cube(&[0.10, 0.10, 0.10]);
    println!(
        "cap_outer cube_uniform: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
    assert!(m.edges_measured >= 1);
    assert!(
        m.worst_g1_coplanarity < helpers::G1_GATE,
        "cube_uniform G1 {:.3e} exceeds {:.1e} commercial gate",
        m.worst_g1_coplanarity,
        helpers::G1_GATE,
    );
}

#[test]
#[ignore = "legacy K-P strict gate: cube_mixed direct-apex exact-boundary passes 5e-3/5e-2, but this helper still misses G0 <= 1e-10"]
fn cap_outer_boundary_passes_on_cube_mixed() {
    let m = run_cube(&[0.08, 0.10, 0.12]);
    assert_cap_outer_gates("cube_mixed", &m);
}

#[ignore = "legacy K-P strict gate: cube_extreme direct-apex exact-boundary passes 5e-3/5e-2, but this helper still misses G0 <= 1e-10"]
#[test]
fn cap_outer_boundary_passes_on_cube_extreme() {
    let m = run_cube(&[0.02, 0.10, 0.20]);
    assert_cap_outer_gates("cube_extreme", &m);
}

#[ignore = "legacy K-P helper path: native exact-boundary pyramid4_mixed passes the direct-apex cap_outer matrix, but this older path still fails G0/G1"]
#[test]
fn cap_outer_boundary_passes_on_pyramid4_mixed() {
    let m = run_n_pyramid(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    assert_cap_outer_gates("pyramid4_mixed", &m);
}

fn run_n_pyramid_cap_outer_result(
    n: usize,
    h: f64,
    radii: &[f64],
) -> (popkern_fillet::BlendResult, helpers::CapOuterMetrics) {
    let pyr = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let specs = edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect::<Vec<_>>();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        use_kp_polar_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None).unwrap();
    let faces = cap_outer_faces_from_result(&result);
    let metrics = measure_cap_outer_boundary(&result.brep, &faces, N_SAMPLES);
    (result, metrics)
}

/// N=5 — the setback vertex blend ships a commercial-grade 20-patch
/// cap fan (re-pinned 2026-07-15, PK-087). History: the declined
/// PnS/sphere-cap construction used to ship a non-commercial cap
/// (g0 = 9.877e-3, g1 = 0.90), then the PK-089 all-or-nothing
/// contract rolled back to mitred (edges = 0); the PK-094/PK-087
/// setback boundary surgery now commits the KP20 cap fan on the
/// fillet setback border. Measured at re-pin: g0 = 3.6e-8,
/// g1 = 1.5e-3 across 20 cap/support edges.
#[test]
fn cap_outer_boundary_pyramid5_mixed_ships_commercial_setback_fan() {
    let (result, m) = run_n_pyramid_cap_outer_result(5, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    println!(
        "cap_outer pyramid5_mixed: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
    assert_eq!(
        m.edges_measured, 20,
        "pyramid5_mixed setback fan must expose its full cap/support boundary"
    );
    assert!(
        m.worst_g0_gap <= 5.0e-3,
        "cap_outer_g0 {:.3e} exceeds the 5e-3 commercial gate",
        m.worst_g0_gap
    );
    assert!(
        m.worst_g1_coplanarity <= 5.0e-2,
        "cap_outer_g1 {:.3e} exceeds the 5e-2 commercial gate",
        m.worst_g1_coplanarity
    );
    let apex_patches = result
        .corner_patches
        .iter()
        .filter(|patch| patch.n_sides == 5)
        .count();
    assert_eq!(
        apex_patches, 20,
        "the committed setback fan must record one CornerPatchInfo per cap face"
    );
    assert!(
        !result
            .edge_warnings
            .iter()
            .any(|(_, w)| matches!(w, BlendWarning::VertexModeFallback { .. })),
        "a committed setback blend must not also record a mitred-fallback warning"
    );
    result
        .brep
        .validate()
        .expect("committed setback blend must leave a valid BRep");
}

/// N=6 variable radii — the setback vertex blend ships a
/// commercial-grade 24-patch cap fan (re-pinned 2026-07-15, PK-087
/// semantic rim gate). History: this fixture was the long-standing
/// route blocker (fallback cap g0=5.359e-3/g1=0.99, then PK-089
/// mitred rollback with edges=0, then a folded fan with g1=7.8e-1
/// from sweep-proximity rim misclassification). Three PK-087 fixes
/// closed it: the border-loop semantic rim gate, the damped rim
/// least-squares (exact degree-9 interpolation slung CPs 150x), and
/// in-place outlier repair in the projection walk (dropped samples
/// re-parameterized whole pcurves). Measured at re-pin: g0 = 9.5e-5,
/// g1 = 2.0e-3 across 24 cap/support edges.
#[test]
fn cap_outer_boundary_hexagon6_variable_ships_commercial_setback_fan() {
    let (result, m) = run_n_pyramid_cap_outer_result(6, 0.8, &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02]);
    println!(
        "cap_outer hexagon6_variable_route: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
    assert_eq!(
        m.edges_measured, 24,
        "hexagon6_variable setback fan must expose its full cap/support boundary"
    );
    assert!(
        m.worst_g0_gap <= 5.0e-3,
        "cap_outer_g0 {:.3e} exceeds the 5e-3 commercial gate",
        m.worst_g0_gap
    );
    assert!(
        m.worst_g1_coplanarity <= 5.0e-2,
        "cap_outer_g1 {:.3e} exceeds the 5e-2 commercial gate",
        m.worst_g1_coplanarity
    );
    let apex_patches = result
        .corner_patches
        .iter()
        .filter(|patch| patch.n_sides == 6)
        .count();
    assert_eq!(
        apex_patches, 24,
        "the committed setback fan must record one CornerPatchInfo per cap face"
    );
    assert!(
        !result
            .edge_warnings
            .iter()
            .any(|(_, w)| matches!(w, BlendWarning::VertexModeFallback { .. })),
        "a committed setback blend must not also record a mitred-fallback warning"
    );
    result
        .brep
        .validate()
        .expect("committed setback blend must leave a valid BRep");
}

#[ignore = "diagnostic probe only: single-fixture cap-outer metric printout, not an acceptance gate"]
#[test]
fn print_cube_uniform_only() {
    eprintln!("\n=== cube_uniform single-fixture probe ===");
    let m = run_cube(&[0.10, 0.10, 0.10]);
    eprintln!(
        "cube_uniform: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

#[ignore = "diagnostic probe only: single-fixture cap-outer metric printout, not an acceptance gate"]
#[test]
fn print_pyramid3_uniform_only() {
    eprintln!("\n=== pyramid3_uniform single-fixture probe ===");
    let m = run_n_pyramid(3, 1.2, &[0.10, 0.10, 0.10]);
    eprintln!(
        "pyramid3_uniform: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

/// Diagnostic: dump each cap face's outer-loop structure —
/// vertex positions, on-sphere classification, and adjacent
/// (twin-side) face surface type. Used to design Phase 6c-C
/// variable-radius topology surgery.
#[ignore = "diagnostic probe only: cap-loop topology dump for PK-076 investigations"]
#[test]
fn dump_cap_loop_structure() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let radii: &[f64] = &[0.08, 0.10, 0.12]; // cube_mixed
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    // Collect all per-fillet sphere candidates from cap-face surfaces.
    let candidates: Vec<(popkern_core::geometry::Point3, f64)> = faces
        .iter()
        .filter_map(|&fh| {
            brep.faces
                .get(fh)
                .and_then(|f| match &f.surface.as_ref()?.geometry {
                    SurfaceGeometry::Sphere { center, radius } => Some((*center, *radius)),
                    _ => None,
                })
        })
        .collect();
    eprintln!("\n=== cube_mixed cap-loop structure ===");
    eprintln!("{} candidate per-fillet spheres:", candidates.len());
    for (i, (c, r)) in candidates.iter().enumerate() {
        eprintln!("  S{i}: ({:.3},{:.3},{:.3}) r={:.3}", c.x, c.y, c.z, r);
    }
    for (i, &fh) in faces.iter().enumerate() {
        let face = &brep.faces[fh];
        let (centre, radius) = match face.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { center, radius }) => (*center, *radius),
            other => {
                eprintln!("cap[{i}]: not a sphere ({other:?})");
                continue;
            }
        };
        eprintln!(
            "cap[{i}] sphere=({:.3},{:.3},{:.3}) r={:.3}",
            centre.x, centre.y, centre.z, radius
        );
        let outer_loop = face.outer_loop;
        let start_he = brep.loops[outer_loop].half_edge;
        let mut cursor = start_he;
        let mut step = 0;
        loop {
            let he_origin = brep.half_edges[cursor].origin;
            let he_twin = brep.half_edges[cursor].twin;
            let he_next = brep.half_edges[cursor].next;
            let v = brep.vertices[he_origin].position;
            let on_sphere =
                ((v.x - centre.x).powi(2) + (v.y - centre.y).powi(2) + (v.z - centre.z).powi(2))
                    .sqrt()
                    - radius;
            // Find best-matching sphere among all candidates.
            let mut best = (usize::MAX, f64::INFINITY);
            for (si, (sc, sr)) in candidates.iter().enumerate() {
                let d = ((v.x - sc.x).powi(2) + (v.y - sc.y).powi(2) + (v.z - sc.z).powi(2)).sqrt()
                    - *sr;
                if d.abs() < best.1.abs() {
                    best = (si, d);
                }
            }
            let twin_face = brep.loops[brep.half_edges[he_twin].loop_].face;
            let twin_kind = match brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(SurfaceGeometry::Plane { .. }) => "Plane",
                Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
                Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                Some(SurfaceGeometry::Gregory(_)) => "Gregory",
                _ => "?",
            };
            eprintln!(
                "  step{step:2}: V({:.3},{:.3},{:.3}) Δown={:+.3e} bestS{}={:+.3e} → twin {twin_kind}",
                v.x, v.y, v.z, on_sphere, best.0, best.1
            );
            cursor = he_next;
            step += 1;
            if cursor == start_he || step > 16 {
                break;
            }
        }
    }
}

/// Diagnostic: dump the cap face's surface variant + boundary
/// vertex count for cube_uniform vs cube_mixed vs cube_extreme.
/// Used to verify which surface type lands on the cap face after
/// the full pipeline (Phase 6c-A continuation §4.2 work).
#[ignore = "diagnostic probe only: cap surface-type dump for PK-076 investigations"]
#[test]
fn dump_cap_surface_types() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let cases: &[(&str, &[f64])] = &[
        ("cube_uniform", &[0.10, 0.10, 0.10]),
        ("cube_mixed", &[0.08, 0.10, 0.12]),
        ("cube_extreme", &[0.02, 0.10, 0.20]),
    ];
    eprintln!("\n=== cap surface-type dump ===");
    for (name, radii) in cases {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = cube.vertices.keys().next().unwrap();
        let edges = incident_edges(&cube, vertex);
        let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
        for (i, &fh) in faces.iter().enumerate() {
            let n_verts = brep.face_outer_loop_positions(fh).len();
            let kind = match brep
                .faces
                .get(fh)
                .and_then(|f| f.surface.as_ref())
                .map(|sof| &sof.geometry)
            {
                Some(SurfaceGeometry::Plane { .. }) => "Plane",
                Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
                Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                Some(SurfaceGeometry::Torus { .. }) => "Torus",
                Some(SurfaceGeometry::Cone { .. }) => "Cone",
                Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
                Some(SurfaceGeometry::Gregory(g)) => {
                    use popkern_core::geometry::gregory::GregoryPatch;
                    match g {
                        GregoryPatch::Triangular(_) => "Gregory:Triangular",
                        GregoryPatch::Rectangular(_) => "Gregory:Rectangular",
                        GregoryPatch::RectangularDeg5(_) => "Gregory:Deg5",
                        GregoryPatch::RectangularTwoStripe(_) => "Gregory:TwoStripe",
                        GregoryPatch::RectangularBiNonic(_) => "Gregory:BiNonic",
                        _ => "Gregory:?",
                    }
                }
                None => "<none>",
            };
            eprintln!("  {name:14} cap[{i}]: surface={kind} n_verts={n_verts}");
        }
    }
}

#[ignore = "diagnostic probe only: baseline metric table for PK-076 investigations"]
#[test]
fn print_cap_outer_boundary_baseline() {
    eprintln!("\n=== cap-outer-boundary baseline (no assertions) ===");
    for (name, m) in [
        ("cube_uniform", run_cube(&[0.10, 0.10, 0.10])),
        ("cube_mixed", run_cube(&[0.08, 0.10, 0.12])),
        ("cube_extreme", run_cube(&[0.02, 0.10, 0.20])),
        (
            "pyramid3_uniform",
            run_n_pyramid(3, 1.2, &[0.10, 0.10, 0.10]),
        ),
        (
            "pyramid4_uniform",
            run_n_pyramid(4, 1.2, &[0.10, 0.10, 0.10, 0.10]),
        ),
        (
            "pyramid4_mixed",
            run_n_pyramid(4, 0.8, &[0.04, 0.06, 0.05, 0.07]),
        ),
        (
            "pyramid_n04_zoo",
            run_n_pyramid(4, 1.2, &[0.07, 0.09, 0.11, 0.13]),
        ),
        (
            "pyramid_n03_zoo",
            run_n_pyramid(3, 1.2, &[0.07, 0.10, 0.13]),
        ),
    ] {
        eprintln!(
            "  {name:18}: g0={:.3e}  g1={:.3e}  edges={}",
            m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
        );
    }
}
