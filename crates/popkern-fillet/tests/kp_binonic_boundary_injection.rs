//! Phase 6c-D Stage 1 — boundary-injection probe.
//!
//! Validates that, with `POPKERN_SPHERE_FIT_BOUNDARY=1` set AND
//! `use_kp_binonic_gregory = true`, the K-P bi-nonic ring builder
//! consumes the sphere-fit-substituted hex-side boundary CPs
//! (`BlendBoundary::hex_boundary_cps_deg9`) instead of falling
//! back to `k_p_boundary` Hermite synthesis.
//!
//! # Stage 1 status: plumbing only
//!
//! Stage 1 *injects* the new boundary but does **not** re-fit the
//! patch's interior CPs to match. Empirically, the resulting
//! bi-nonic patch self-intersects (the row0/col9 boundary moved
//! but the strict-interior 64 CPs — still Hermite-blend output —
//! produced a fold). The promotion ladder rejects the fold via
//! `self_intersects_with_tolerance` and falls back to a lower
//! rung. Net effect: with injection on, fewer fixtures land
//! binonic than without.
//!
//! This is *expected*. The empirical fall-off proves the injection
//! is reaching the K-P guide builder — without injection the
//! patches don't fold. Stage 2 (re-fit strict-interior CPs against
//! the injected boundary) closes this; once it lands, binonic
//! lands again *and* the cap-outer boundary tracks the analytic
//! arc to the polynomial-fit floor.
//!
//! # Matrix probe
//!
//! `boundary_injection_matrix_probe` (ignored by default) prints
//! a side-by-side comparison of inject=false vs inject=true on
//! cube_mixed and pyramid4_mixed. Run with
//!
//! ```text
//! cargo test -p popkern-fillet --test kp_binonic_boundary_injection \
//!   -- --ignored --nocapture
//! ```
//!
//! to read the binonic-landing-rate and cap-outer continuity numbers.
//!
//! # Env-var safety
//!
//! Each test that sets `POPKERN_SPHERE_FIT_BOUNDARY` uses an
//! `EnvGuard` so the prior value is restored even if the test
//! panics — same-process other tests aren't perturbed.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    N_SAMPLES, collect_patches, fillet_corner, incident_edges, measure_cap_outer_boundary,
    pick_vertex_with_valence,
};

/// RAII guard restoring the prior value of an env-var on drop.
struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        unsafe {
            std::env::set_var(key, value);
        }
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.prior {
                Some(v) => std::env::set_var(self.key, v),
                None => std::env::remove_var(self.key),
            }
        }
    }
}

struct ProbeResult {
    n_corner_faces: usize,
    n_gregory_patches: usize,
    n_binonic: usize,
    cap_outer_g0: f64,
    cap_outer_g1: f64,
    edges_measured: usize,
}

impl ProbeResult {
    fn print(&self, fixture: &str) {
        eprintln!(
            "  {fixture:18}: faces={} gregory={} binonic={} \
             cap_outer_g0={:.3e} cap_outer_g1={:.3e} edges={}",
            self.n_corner_faces,
            self.n_gregory_patches,
            self.n_binonic,
            self.cap_outer_g0,
            self.cap_outer_g1,
            self.edges_measured,
        );
    }
}

fn probe_cube(radii: &[f64], inject: bool) -> ProbeResult {
    let _guard = inject.then(|| EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1"));
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    finalize(&brep, &faces)
}

fn probe_square_pyramid(radii: &[f64], inject: bool) -> ProbeResult {
    let _guard = inject.then(|| EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1"));
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    finalize(&brep, &faces)
}

fn finalize(brep: &BRep, faces: &[popkern_core::topology::FaceHandle]) -> ProbeResult {
    let patches = collect_patches(brep, faces);
    let metrics = measure_cap_outer_boundary(brep, faces, N_SAMPLES);
    let n_binonic = patches
        .iter()
        .filter(|p| matches!(p, GregoryPatch::RectangularBiNonic(_)))
        .count();
    ProbeResult {
        n_corner_faces: faces.len(),
        n_gregory_patches: patches.len(),
        n_binonic,
        cap_outer_g0: metrics.worst_g0_gap,
        cap_outer_g1: metrics.worst_g1_coplanarity,
        edges_measured: metrics.edges_measured,
    }
}

/// Side-by-side comparison: cap-outer continuity and binonic
/// landing rate with vs without boundary injection. Diagnostic
/// only — Stage 1 expectation is that injection regresses the
/// binonic landing rate (folded interior, see module header).
/// Stage 2's interior re-fit will turn this into a real gate.
#[ignore]
#[test]
fn boundary_injection_matrix_probe() {
    eprintln!("\n=== Phase 6c-D Stage 1 boundary-injection matrix ===");
    eprintln!("(use_kp_binonic_gregory=true; ±POPKERN_SPHERE_FIT_BOUNDARY)");
    eprintln!("--- inject=false (baseline: Hermite-synth boundary) ---");
    let cm_off = probe_cube(&[0.08, 0.10, 0.12], false);
    cm_off.print("cube_mixed");
    let p4m_off = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], false);
    p4m_off.print("pyramid4_mixed");
    eprintln!("--- inject=true (Stage 1: substituted boundary) ---");
    let cm_on = probe_cube(&[0.08, 0.10, 0.12], true);
    cm_on.print("cube_mixed");
    let p4m_on = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], true);
    p4m_on.print("pyramid4_mixed");
}

/// Stage 1 wiring sanity: verify the *baseline* (no injection)
/// path still ships binonic on cube_mixed. Guards against the
/// Stage 1 plumbing accidentally regressing the inject-OFF path —
/// the `Option<&[...]>` parameter must thread `None` through
/// cleanly when the boundary's `hex_boundary_cps_deg9` is empty.
#[test]
fn cube_mixed_baseline_binonic_unaffected_by_stage1_plumbing() {
    let res = probe_cube(&[0.08, 0.10, 0.12], false);
    res.print("cube_mixed[no-inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 1 plumbing broke the inject=OFF baseline — cube_mixed should still ship binonic"
    );
}

/// Same baseline guard for pyramid4_mixed (the fixture that
/// already passes commercial gates per the H.7 matrix). With no
/// injection the K-P binonic dispatch must continue to land.
#[test]
fn pyramid4_mixed_baseline_binonic_unaffected_by_stage1_plumbing() {
    let res = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], false);
    res.print("pyramid4_mixed[no-inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 1 plumbing broke the inject=OFF baseline — pyramid4_mixed should still ship binonic"
    );
}

/// Phase 6c-D Stage 3 lock-in: under injection, both cube_mixed
/// and pyramid4_mixed must ship full binonic rings. The
/// row-9-based mid-edge formula (`row1[j] = (8/9)·row0[j] +
/// (1/9)·row9[j]`) + Coons interior fill produces fold-free
/// patches by construction, and the quality regression guard is
/// skipped under injection (the spoke-G1 metric is determined by
/// boundary geometry, not the L-M optimum, so the 1.2× threshold
/// doesn't apply).
///
/// Note: cap-outer continuity is **not** asserted yet. Stage 3
/// proves the dispatch ladder accepts injection-mode binonic; the
/// 3-D cap-vs-cylinder fidelity needs a follow-up (likely pcurve-
/// refit under the new boundary or wider edge-match tolerance in
/// `wire_substituted_boundary_to_edges`).
#[test]
fn cube_mixed_binonic_ships_under_injection() {
    let res = probe_cube(&[0.08, 0.10, 0.12], true);
    res.print("cube_mixed[inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 3 regression: cube_mixed under injection should ship at least one binonic patch"
    );
}

#[test]
fn pyramid4_mixed_binonic_ships_under_injection() {
    let res = probe_square_pyramid(&[0.04, 0.06, 0.05, 0.07], true);
    res.print("pyramid4_mixed[inject]");
    assert!(
        res.n_binonic >= 1,
        "Stage 3 regression: pyramid4_mixed under injection should ship at least one binonic patch"
    );
}

/// Phase 6c-D Stage 4 diagnostic — evaluate the Gregory cap
/// surface DIRECTLY at u=0 / v=1 boundaries (its row 0 / col 9)
/// and check whether those match the analytic FilletEnd /
/// FaceEdge curves the injection wrote into the BRep edges.
///
/// This bypasses the pcurve refit (which was producing UVs
/// outside [0, 1]² in the per-edge breakdown). If Gregory(0, v)
/// matches the analytic curve on the row 0 boundary, the cap
/// surface IS tracking the analytic geometry — the cap-outer
/// metric failure must be in the pcurve refit, not the surface.
#[ignore]
#[test]
fn cap_surface_at_row0_tracks_analytic_boundary() {
    use popkern_core::geometry::gregory::GregoryPatch as GP;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let _guard = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii: Vec<f64> = vec![0.04, 0.06, 0.05, 0.07];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);

    eprintln!("\n=== pyramid4_mixed cap-surface direct eval at row 0 ===");
    for (i, &fh) in faces.iter().enumerate() {
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            eprintln!("  cap[{i}]: no surface");
            continue;
        };
        let SurfaceGeometry::Gregory(GP::RectangularBiNonic(g)) = &sof.geometry else {
            eprintln!("  cap[{i}]: not RectangularBiNonic");
            continue;
        };
        // Sample row 0 (u=0, v varies) at v = 0.0, 0.5, 1.0.
        let p_v0 = g.point_at(0.0, 0.0);
        let p_mid = g.point_at(0.0, 0.5);
        let p_v1 = g.point_at(0.0, 1.0);
        eprintln!(
            "  cap[{i}] row0:  v=0 ({:.4},{:.4},{:.4})  v=.5 ({:.4},{:.4},{:.4})  v=1 ({:.4},{:.4},{:.4})",
            p_v0.x, p_v0.y, p_v0.z, p_mid.x, p_mid.y, p_mid.z, p_v1.x, p_v1.y, p_v1.z,
        );
        // Sample col 9 (v=1, u varies) at u = 0.0, 0.5, 1.0.
        let p_u0 = g.point_at(0.0, 1.0);
        let p_mid2 = g.point_at(0.5, 1.0);
        let p_u1 = g.point_at(1.0, 1.0);
        eprintln!(
            "  cap[{i}] col9:  u=0 ({:.4},{:.4},{:.4})  u=.5 ({:.4},{:.4},{:.4})  u=1 ({:.4},{:.4},{:.4})",
            p_u0.x, p_u0.y, p_u0.z, p_mid2.x, p_mid2.y, p_mid2.z, p_u1.x, p_u1.y, p_u1.z,
        );
    }

    // For each cap face, find its outer-loop edges that border a
    // Cylinder twin. Check whether the cap's row-0-or-col-9
    // surface evaluation tracks the edge's 3D curve.
    eprintln!("\n=== checking each cap-vs-cylinder edge ===");
    use popkern_core::geometry::curve::CurveEval;
    for &fh in &faces {
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Gregory(GP::RectangularBiNonic(g)) = &sof.geometry else {
            continue;
        };
        let outer = face.outer_loop;
        let start = brep.loops[outer].half_edge;
        let mut cur = start;
        let mut step = 0;
        loop {
            let he = &brep.half_edges[cur];
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            let twin = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref());
            let is_cyl = matches!(
                twin.map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            );
            if is_cyl {
                let edge = &brep.edges[he.edge];
                if let Some(coe) = edge.curve.as_ref() {
                    let t_lo = coe.t_start;
                    let t_hi = coe.t_end;
                    let t_mid = 0.5 * (t_lo + t_hi);
                    let p3 = coe.geometry.point_at(t_mid);
                    // Find closest u, v on the Gregory cap by sampling
                    // boundary (u=0 with v∈[0,1] or v=1 with u∈[0,1]).
                    let mut best = (f64::INFINITY, 0.0_f64, 0u8);
                    for k in 0..=20 {
                        let s = (k as f64) / 20.0;
                        // u=0 row
                        let p = g.point_at(0.0, s);
                        let d = p3.distance_to(&p);
                        if d < best.0 {
                            best = (d, s, 0);
                        }
                        // v=1 col
                        let p = g.point_at(s, 1.0);
                        let d = p3.distance_to(&p);
                        if d < best.0 {
                            best = (d, s, 1);
                        }
                    }
                    let bnd = if best.2 == 0 { "row0" } else { "col9" };
                    eprintln!(
                        "  cap face step{step}: edge_chord={:.3e} t_mid={t_mid:.3e}  cap_to_edge_min={:.3e} on {bnd}@{:.2}",
                        coe.geometry
                            .point_at(t_lo)
                            .distance_to(&coe.geometry.point_at(t_hi)),
                        best.0,
                        best.1,
                    );
                }
            }
            cur = he.next;
            step += 1;
            if cur == start || step > 16 {
                break;
            }
        }
    }
}

/// Phase 6c-D Stage 4.3 diagnostic — for each cap face, dump its
/// outer-loop edges (handle, vertex positions, curve3d midpoint)
/// alongside cap_outer's view of the same BRep. Reveals which
/// edges the wire step misses vs which edges cap_outer measures.
#[ignore]
#[test]
fn cap_loop_vs_cap_outer_edge_coverage() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let _guard = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let _guard2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii: Vec<f64> = vec![0.04, 0.06, 0.05, 0.07];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

    eprintln!(
        "\n=== cap face outer-loop edges (from wire's perspective) ===\n{} cap faces",
        faces.len()
    );
    let mut wire_edges: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
        std::collections::BTreeSet::new();
    for (face_idx, &fh) in faces.iter().enumerate() {
        let face = &brep.faces[fh];
        let outer = face.outer_loop;
        let start = brep.loops[outer].half_edge;
        eprintln!("\nCap face #{face_idx} ({fh:?}):");
        let mut cur = start;
        for step in 0..16 {
            let he = &brep.half_edges[cur];
            let edge = he.edge;
            wire_edges.insert(edge);
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            let twin_kind = brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Plane { .. } => "Plane",
                    SurfaceGeometry::Cylinder { .. } => "Cyl",
                    SurfaceGeometry::Sphere { .. } => "Sphere",
                    SurfaceGeometry::Gregory(_) => "Gregory",
                    SurfaceGeometry::Nurbs(_) => "Nurbs",
                    _ => "Other",
                })
                .unwrap_or("none");
            let v_origin = brep.vertices[he.origin].position;
            let curve_mid = brep.edges[edge].curve.as_ref().map(|coe| {
                let t_mid = 0.5 * (coe.t_start + coe.t_end);
                coe.geometry.point_at(t_mid)
            });
            eprintln!(
                "  step{step}: he={cur:?} edge={edge:?} twin={twin_kind:7} (face {twin_face:?}) v_origin=({:.3},{:.3},{:.3}) curve_mid={:?}",
                v_origin.x,
                v_origin.y,
                v_origin.z,
                curve_mid.map(|p| (p.x, p.y, p.z))
            );
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }

    eprintln!(
        "\n=== cap_outer's view (BRep edge iteration) ===\nwire visited {} unique edges",
        wire_edges.len()
    );
    let mut cap_outer_edges: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
        std::collections::BTreeSet::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(_coe) = edge.curve.as_ref() else {
            continue;
        };
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        if ha.pcurve.is_none() || hb.pcurve.is_none() {
            continue;
        }
        cap_outer_edges.insert(eh);
    }
    eprintln!("cap_outer measures {} unique edges", cap_outer_edges.len());
    let in_wire_not_outer: Vec<_> = wire_edges.difference(&cap_outer_edges).collect();
    let in_outer_not_wire: Vec<_> = cap_outer_edges.difference(&wire_edges).collect();
    eprintln!(
        "wire-only edges (not measured by cap_outer): {} = {:?}",
        in_wire_not_outer.len(),
        in_wire_not_outer
    );
    eprintln!(
        "cap_outer-only edges (not visited by wire): {} = {:?}",
        in_outer_not_wire.len(),
        in_outer_not_wire
    );
}

/// Phase 6c-D Stage 4 diagnostic: break down `cap_outer_g0` by
/// per-edge contribution. Reports for each cap-vs-non-cap edge:
/// the twin-face surface type, mid-edge G0 gap, and (for cap-side
/// pcurve only) the cap surface UV at midpoint. Helps localise
/// whether the worst gap is on a row 0 (cap-tangent arc), col 9
/// (slant bridge), or some other edge.
#[ignore]
#[test]
fn cap_outer_per_edge_breakdown_pyramid4_mixed() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let _guard = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii: Vec<f64> = vec![0.04, 0.06, 0.05, 0.07];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

    eprintln!("\n=== pyramid4_mixed cap-outer per-edge breakdown ===");
    eprintln!(
        "{} cap faces; {} BRep edges total",
        faces.len(),
        brep.edges.len()
    );

    let mut idx = 0;
    for (_eh, edge) in brep.edges.iter() {
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let Some(pa) = ha.pcurve.as_ref() else {
            continue;
        };
        let Some(pb) = hb.pcurve.as_ref() else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let kind = |g: &SurfaceGeometry| -> &'static str {
            match g {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
            }
        };
        let kind_a = kind(&sa.geometry);
        let kind_b = kind(&sb.geometry);
        let cap_kind = if a_is_cap { kind_a } else { kind_b };
        let twin_kind = if a_is_cap { kind_b } else { kind_a };

        // Sample at edge midpoint.
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let t_mid = 0.5 * (t_lo + t_hi);
        let p3_mid = coe.geometry.point_at(t_mid);
        let p3_lo = coe.geometry.point_at(t_lo);
        let p3_hi = coe.geometry.point_at(t_hi);
        let edge_len_chord = p3_lo.distance_to(&p3_hi);

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let pa_at = |t: f64| {
            let uv = pa.point_at(t);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let pb_at = |t: f64| {
            let uv = pb.point_at(t);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let _ = (
            p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
                <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo)),
            p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
                <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo)),
        );
        // Sample at pcurve midpoint (direction is symmetric for s=0.5).
        let t_a = pa_lo + 0.5 * (pa_hi - pa_lo);
        let t_b = pb_lo + 0.5 * (pb_hi - pb_lo);
        let p_a = pa_at(t_a);
        let p_b = pb_at(t_b);
        let g0 = p_a.distance_to(&p_b);

        // Distance from edge curve midpoint to each surface
        // evaluation. Reveals whether the curve drifts from one
        // or both surfaces.
        let dev_a = p3_mid.distance_to(&p_a);
        let dev_b = p3_mid.distance_to(&p_b);

        idx += 1;
        eprintln!(
            "  edge#{idx} cap={cap_kind} twin={twin_kind} cap_is_a={a_is_cap}  chord={edge_len_chord:.3e}  \
             g0={g0:.3e}  dev_cap={dev_a:.3e}  dev_twin={dev_b:.3e}",
        );
        let uv_a = pa.point_at(t_a);
        let uv_b = pb.point_at(t_b);
        let (cap_uv, twin_uv) = if a_is_cap { (uv_a, uv_b) } else { (uv_b, uv_a) };
        eprintln!(
            "         cap_uv=({:.3},{:.3})  twin_uv=({:.3},{:.3})",
            cap_uv.x, cap_uv.y, twin_uv.x, twin_uv.y
        );
        eprintln!(
            "         pa_range=[{pa_lo:.3},{pa_hi:.3}]  pb_range=[{pb_lo:.3},{pb_hi:.3}]  \
             coe_range=[{t_lo:.3},{t_hi:.3}]"
        );
        // Sample analytic edge curve at lo, mid, hi to compare to surface evals
        eprintln!(
            "         curve3d: lo=({:.3},{:.3},{:.3}) mid=({:.3},{:.3},{:.3}) hi=({:.3},{:.3},{:.3})",
            p3_lo.x, p3_lo.y, p3_lo.z, p3_mid.x, p3_mid.y, p3_mid.z, p3_hi.x, p3_hi.y, p3_hi.z,
        );
        eprintln!(
            "         cap_eval=({:.3},{:.3},{:.3})  twin_eval=({:.3},{:.3},{:.3})",
            p_a.x, p_a.y, p_a.z, p_b.x, p_b.y, p_b.z,
        );
    }
}
