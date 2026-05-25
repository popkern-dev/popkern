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
//! # Status (2026-04-27)
//!
//! All fixtures `#[ignore]` because they're expected to fail
//! today — the gate is being introduced *before* the fix to
//! quantify the problem precisely. Once Task #144 closes, these
//! flip to passing and become regression guards.
//!
//! Threshold targets (mirroring `kp_binonic_commercial_gates`):
//! * G0 ≤ 1e-10 (position match at the seam)
//! * G1 ≤ 5e-2  (coplanarity at samples; ≈ 2.9° normal mismatch)

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    G0_GATE, G1_GATE, N_SAMPLES, fillet_corner, incident_edges, measure_cap_outer_boundary,
    pick_vertex_with_valence,
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
fn cap_outer_boundary_passes_on_cube_mixed() {
    let m = run_cube(&[0.08, 0.10, 0.12]);
    assert_cap_outer_gates("cube_mixed", &m);
}

#[ignore = "PK-080: native direct-apex needs outside-outside support continuation before extreme cube can meet G0/G1 gates"]
#[test]
fn cap_outer_boundary_passes_on_cube_extreme() {
    let m = run_cube(&[0.02, 0.10, 0.20]);
    assert_cap_outer_gates("cube_extreme", &m);
}

#[ignore = "PK-078: pyramid4 mixed cap-outer seam currently fails G0/G1 commercial gates"]
#[test]
fn cap_outer_boundary_passes_on_pyramid4_mixed() {
    let m = run_n_pyramid(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    assert_cap_outer_gates("pyramid4_mixed", &m);
}

/// N=5 — PnS path, not Gregory caps. Helper only counts cap edges
/// where one side is a Gregory face; PnS produces Sphere-surface
/// caps, so this fixture should report 0 edges measured. Kept here
/// so the matrix is documented; will be skipped or replaced when
/// the PnS sibling helper lands.
#[ignore = "PK-079: replace Gregory-only helper with PnS/sphere-cap cap-outer coverage"]
#[test]
fn cap_outer_boundary_pyramid5_mixed_path_check() {
    let m = run_n_pyramid(5, 0.8, &[0.04, 0.06, 0.05, 0.07, 0.05]);
    println!(
        "cap_outer pyramid5_mixed: g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

/// Probe-style: print numbers across the canonical fixtures with
/// no assertions so we can see the baseline without picking
/// thresholds first. Useful while developing the fix. Includes
/// `pyramid_n04_zoo` — the same geometry rendered by the
/// `fillet_zoo` showcase binary, where FreeCAD shows the
/// "messy top" defect.
///
/// G1 is reported but is *not* meaningful while G0 is large —
/// we're sampling each side's surface at 3-D points that are
/// 5-20cm apart, so the normals are essentially uncorrelated.
/// Once G0 drops below ~1e-3, G1 becomes informative.
/// Single-fixture probe for the cube_uniform path — used to trace
/// the per-vertex insert_vertex_blend invocation count without the
/// noise of the full 8-fixture matrix.
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
