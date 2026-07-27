//! PK-192 — ignored diagnosis probes for the low-valence PnS caps.
//!
//! These are the campaign's measurement instruments, split out of
//! `pns_low_valence_gates.rs` (which keeps the pass/fail gates). All
//! are `#[ignore]`d: each doc comment carries its run command, and
//! several read process-wide env vars so they must run in their own
//! process.
//!
//! Every probe reports through `run_pyramid_apex_outcome`, the same
//! harness the gates measure with.

mod kp_gate_helpers;
#[path = "support/pns_apex_outcome.rs"]
mod pns_apex_outcome;

use pns_apex_outcome::run_pyramid_apex_outcome;

/// PK-192 diagnosis: what actually bounds a committed cap face?
///
/// The PnS ring is built from each sector's (fillet-end arc + face
/// bridge) pair on the assumption that those two edges span the
/// sector's share of the rim. They do not — the ring failed to close
/// by about half a fillet radius. This walks the committed cap faces'
/// own loops and names every bounding edge by the face on its far
/// side, which says what occupies the gap.
#[test]
#[ignore = "PK-192 diagnosis: run with --ignored --nocapture"]
fn probe_cap_face_loop_anatomy() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;

    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        let brep = &o.brep;
        println!("--- n={n} cap faces={} ---", o.cap_face_handles.len());
        for (i, &face) in o.cap_face_handles.iter().enumerate() {
            let Some(face_data) = brep.faces.get(face) else {
                continue;
            };
            let start = brep.loops[face_data.outer_loop].half_edge;
            let mut he = start;
            let mut idx = 0;
            loop {
                let d = &brep.half_edges[he];
                let other = brep
                    .half_edges
                    .get(d.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| l.face);
                let kind = other
                    .and_then(|f| brep.faces.get(f))
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| match &s.geometry {
                        SurfaceGeometry::Plane { .. } => "Plane",
                        SurfaceGeometry::Cylinder { .. } => "Cylinder",
                        SurfaceGeometry::Sphere { .. } => "Sphere",
                        SurfaceGeometry::Torus { .. } => "Torus",
                        SurfaceGeometry::Cone { .. } => "Cone",
                        SurfaceGeometry::Nurbs(_) => "Nurbs",
                        SurfaceGeometry::Gregory(_) => "Gregory",
                        SurfaceGeometry::PolyhedralSpline(_) => "PnS",
                    })
                    .unwrap_or("<none>");
                let is_cap = other
                    .map(|f| o.cap_face_handles.contains(&f))
                    .unwrap_or(false);
                let len = brep
                    .edges
                    .get(d.edge)
                    .and_then(|e| e.curve.as_ref())
                    .map(|c| {
                        let a = (&c.geometry as &dyn CurveEval).point_at(c.t_start);
                        let b = (&c.geometry as &dyn CurveEval).point_at(c.t_end);
                        a.distance_to(&b)
                    })
                    .unwrap_or(f64::NAN);
                println!(
                    "  cap[{i}] he{idx}: across={kind}{} chord={len:.5e}",
                    if is_cap { " (CAP-CAP seam)" } else { "" }
                );
                he = d.next;
                idx += 1;
                if he == start || idx > 12 {
                    break;
                }
            }
        }
    }
}

/// PK-192: is the support-side pcurve deviation caused by the PnS
/// split, or already in the shipped body?
///
/// The PnS dispatch is refused at `cap_outer_g0` 5.545e-2, whose
/// dominant evidence term is a 4.618e-2 support-curve gap. That is
/// only PnS's fault if the shipped body's own pcurves are tighter.
/// This walks every half-edge of the INCUMBENT result and measures the
/// same quantity: the pcurve image against the edge curve at matched
/// fractions.
#[test]
#[ignore = "PK-192 attribution probe: run with --ignored --nocapture"]
fn probe_incumbent_pcurve_same_fraction_deviation() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        let brep = &o.brep;
        let mut worst = 0.0_f64;
        let mut over_bar = 0usize;
        let mut measured = 0usize;
        for (_eh, edge) in brep.edges.iter() {
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            for he in [edge.half_edge] {
                let Some(hd) = brep.half_edges.get(he) else {
                    continue;
                };
                let Some(pc) = hd.pcurve.as_ref() else {
                    continue;
                };
                let Some(face) = brep.loops.get(hd.loop_).map(|l| l.face) else {
                    continue;
                };
                let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                    continue;
                };
                measured += 1;
                // A half-edge pcurve runs along its LOOP, which for one
                // of the two half-edges is opposite to the edge curve's
                // own direction. Sampling both forward reads a REVERSED
                // pcurve as a deviation of about 0.9 * chord, which is
                // exactly the signature the first cut of this probe
                // produced. Resolve orientation from the endpoints, as
                // `cap_outer_acceptance_measure` does, before measuring.
                let img_at = |t: f64| {
                    let uv = pc.point_at(t);
                    surface.geometry.point_at(uv.x, uv.y)
                };
                let c_lo = (&curve.geometry as &dyn CurveEval).point_at(curve.t_start);
                let c_hi = (&curve.geometry as &dyn CurveEval).point_at(curve.t_end);
                let forward = c_lo.distance_to(&img_at(pc.t_start))
                    + c_hi.distance_to(&img_at(pc.t_end))
                    <= c_lo.distance_to(&img_at(pc.t_end)) + c_hi.distance_to(&img_at(pc.t_start));
                let mut local = 0.0_f64;
                for k in 0..=20 {
                    let s_raw = 0.05 + 0.90 * f64::from(k) / 20.0;
                    let s = if forward { s_raw } else { 1.0 - s_raw };
                    let uv = pc.point_at(pc.t_start + s * (pc.t_end - pc.t_start));
                    let img = surface.geometry.point_at(uv.x, uv.y);
                    let cp = (&curve.geometry as &dyn CurveEval)
                        .point_at(curve.t_start + s_raw * (curve.t_end - curve.t_start));
                    local = local.max(img.distance_to(&cp));
                }
                worst = worst.max(local);
                if local > 5.0e-3 {
                    over_bar += 1;
                    let kind = match &surface.geometry {
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                        popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => {
                            "Cylinder"
                        }
                        popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                        popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
                        popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
                        popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_) => "Nurbs",
                        popkern_core::geometry::surface::SurfaceGeometry::Gregory(_) => "Gregory",
                        popkern_core::geometry::surface::SurfaceGeometry::PolyhedralSpline(_) => {
                            "PnS"
                        }
                    };
                    let cd = pc.point_at(pc.t_start);
                    let cd2 = pc.point_at(pc.t_end);
                    println!(
                        "[pk192-attrib]   OVER: face_kind={kind} dev={local:.3e} \
                         pcurve_uv_span=({:.4},{:.4})->({:.4},{:.4}) chord={:.3e}",
                        cd.x,
                        cd.y,
                        cd2.x,
                        cd2.y,
                        (&curve.geometry as &dyn CurveEval)
                            .point_at(curve.t_start)
                            .distance_to(
                                &(&curve.geometry as &dyn CurveEval).point_at(curve.t_end)
                            )
                    );
                }
            }
        }
        println!(
            "[pk192-attrib] n={n} pcurves_measured={measured} worst_same_fraction_dev={worst:.3e} over_5e-3={over_bar}"
        );
    }
}

/// PK-192: is the mid-pipeline support-face reading a genuine defect,
/// or an unfinished state that later passes repair?
///
/// PnS dispatch runs before the pcurve-completion passes and measures
/// 3.277 of pcurve deviation on a neighbouring cylinder. If the same
/// faces read clean at pipeline END, the acceptance gate is being
/// applied too early rather than the geometry being wrong — a
/// different fix entirely.
#[test]
#[ignore = "PK-192 timing probe: run with --ignored --nocapture"]
fn probe_support_face_fidelity_at_pipeline_end() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_fillet::pcurve_fidelity::measure_pcurve_fidelity;

    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        let brep = &o.brep;
        let cylinders: Vec<_> = brep
            .faces
            .iter()
            .filter(|(_, d)| {
                matches!(
                    d.surface.as_ref().map(|s| &s.geometry),
                    Some(SurfaceGeometry::Cylinder { .. })
                )
            })
            .map(|(f, _)| f)
            .collect();
        let all: Vec<_> = brep.faces.iter().map(|(f, _)| f).collect();
        println!(
            "[pk192-timing] n={n} cylinders={} cyl_fidelity={:?}",
            cylinders.len(),
            measure_pcurve_fidelity(brep, &cylinders).map(|f| f.worst)
        );
        println!(
            "[pk192-timing] n={n} all_faces_fidelity={:?}",
            measure_pcurve_fidelity(brep, &all).map(|f| f.worst)
        );
        // PK-203: the same quantity the acceptance gate refuses on,
        // measured with NO pcurve in between — does each edge's curve
        // lie on its neighbouring face's surface? Mid-pipeline this
        // reads 4.632e-2. If it is small here, that state is transient
        // and the gate is simply running too early.
        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::surface::SurfaceEval;
        let mut worst_off = 0.0_f64;
        for &face in &all {
            let Some(fd) = brep.faces.get(face) else {
                continue;
            };
            let Some(sf) = fd.surface.as_ref() else {
                continue;
            };
            let Some(ld) = brep.loops.get(fd.outer_loop) else {
                continue;
            };
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..32 {
                let Some(hd) = brep.half_edges.get(he) else {
                    break;
                };
                if let Some(curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()) {
                    for k in 0..=20 {
                        let t = curve.t_start
                            + (0.05 + 0.90 * f64::from(k) / 20.0) * (curve.t_end - curve.t_start);
                        let p = (&curve.geometry as &dyn CurveEval).point_at(t);
                        if let Some((u, v)) = sf.geometry.closest_point(&p) {
                            worst_off = worst_off.max(sf.geometry.point_at(u, v).distance_to(&p));
                        }
                    }
                }
                he = hd.next;
                if he == start {
                    break;
                }
            }
        }
        println!("[pk192-timing] n={n} FINAL worst_curve_off_surface={worst_off:.3e}");
    }
}

/// Rim versus spoke on the PnS cap ITSELF, mid-pipeline.
///
/// This is the measurement `probe_cap_rim_versus_spoke_deviation`
/// below could not make: that probe runs on the final body, where the
/// PnS cap was rolled back, so it reads the incumbent. The split here
/// is computed INSIDE the dispatch (`cap_rim_spoke_split` diagnostic,
/// emitted at the same moment the acceptance gate judges the cap) and
/// surfaces through the structured diagnostics channel.
///
/// Reads process-wide `POPKERN_PNS_FIRST_N3`, so it must run in its
/// own process and is ignored by default:
///
/// ```text
/// POPKERN_PNS_FIRST_N3=1 cargo test -p popkern-fillet \
///     --test pns_low_valence_probes probe_pns_mid_pipeline_rim_spoke_split \
///     -- --ignored --nocapture
/// ```
#[test]
#[ignore = "PK-192 rim/spoke probe: reads process-wide POPKERN_PNS_FIRST_N3; see doc comment for the run command"]
fn probe_pns_mid_pipeline_rim_spoke_split() {
    // Same read as `pns_promote::pns_first_at_low_valence` — the test
    // only labels its output with it; routing reads the real one.
    let released = std::env::var("POPKERN_PNS_FIRST_N3")
        .map(|v| v != "0")
        .unwrap_or(false);
    // The uniform-radius apex collapses to the mitred/bridge
    // construction and never reaches the dispatch (zero corner
    // patches, zero pns diagnostics) — printed for the record, not
    // asserted. The mixed-radius apex is the fixture every PK-192
    // number was measured on, and there the split must be recorded.
    for (n, radii, must_record) in [
        (3_usize, vec![0.05, 0.05, 0.05], false),
        (3, vec![0.04, 0.06, 0.05], true),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        println!(
            "[pk192-mid-split] n={n} radii={radii:?} released_for_pns={released} served_by_pns={}",
            o.served_by_pns()
        );
        for d in &o.pns_diagnostics {
            println!("[pk192-mid-split]     diag: {d}");
        }
        if released && must_record {
            assert!(
                o.pns_diagnostics
                    .iter()
                    .any(|d| d.contains("cap_rim_spoke_split")),
                "n={n}: with N=3 released for PnS the dispatch must record the \
                 rim/spoke split before its verdict; diagnostics were {:?}",
                o.pns_diagnostics
            );
        } else if !released {
            println!(
                "[pk192-mid-split] POPKERN_PNS_FIRST_N3 not set — N=3 stays on the \
                 direct-apex route and the PnS dispatch never runs. Set the env var \
                 per the doc comment; this run measured nothing."
            );
        }
    }
}

/// PK-197 criterion: zebra evidence across the internal cap seams of
/// the n=3 and n=4 mixed apexes — the stripe-phase form of the zebra
/// render, measured where a render would be eyeballed. A creased seam
/// reads a large stripe-phase jump; a G1 seam reads near zero.
///
/// Asserted (not just printed): with `harmonize_mismatched_spoke_nets`
/// + `cross_spoke_g1` the caps must stay below a phase jump of 0.25
/// stripes on every internal seam. Pre-fix, the n=3 crease (89.9 deg)
/// saturated the measure.
#[test]
fn zebra_refutes_the_low_valence_cap_crease() {
    use popkern_core::geometry::Vec3;
    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        assert!(o.metrics.seam_count > 0, "n={n}: no internal seams");
        let brep = &o.brep;
        let caps: std::collections::BTreeSet<_> = o.cap_face_handles.iter().copied().collect();
        let mut measured = 0usize;
        let mut worst = 0.0_f64;
        for (_eh, edge) in brep.edges.iter() {
            let Some(he) = brep.half_edges.get(edge.half_edge) else {
                continue;
            };
            let Some(twin) = brep.half_edges.get(he.twin) else {
                continue;
            };
            let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_))
            else {
                continue;
            };
            let (fa, fb) = (la.face, lb.face);
            if fa == fb || !caps.contains(&fa) || !caps.contains(&fb) {
                continue;
            }
            if let Some(kink) = popkern_fillet::zebra::zebra_boundary_kink(
                brep,
                fa,
                fb,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(1.0, 0.0, 0.0),
                0.02,
                24,
                1e-6,
            ) {
                measured += 1;
                worst = worst.max(kink);
            }
        }
        println!(
            "[pk197-zebra] n={n} seams_measured={measured} worst_stripe_phase_jump={worst:.4}"
        );
        assert!(measured > 0, "n={n}: zebra measured no internal seams");
        assert!(
            worst < 0.25,
            "n={n}: internal cap seam shows a zebra stripe-phase jump of {worst:.4} \
             (creased); the cross-spoke + harmonization fix should have refuted this"
        );
    }
}

/// Rim versus spoke pcurve deviation on the SHIPPED cap.
///
/// Note what this does and does not measure. It runs on the final
/// body, where PnS never committed — so these are the *incumbent*
/// cap's numbers at pipeline end, not the PnS cap's mid-pipeline
/// ones. It was written to split PnS's 0.0796 into rim and spoke and
/// cannot do that; separating those requires instrumenting inside the
/// dispatch, before the rollback — which is what
/// `probe_pns_mid_pipeline_rim_spoke_split` above now does.
///
/// What it does establish is that the shipped cap's own pcurves are
/// far under the 5.000e-3 bar on both rim and spoke at pipeline end —
/// another instance of the pattern this ticket keeps meeting: the
/// geometry is fine by the time the pipeline finishes, and the PnS
/// gate judges it long before that.
#[test]
#[ignore = "PK-192 attribution: run with --ignored --nocapture"]
fn probe_cap_rim_versus_spoke_deviation() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    for (n, radii) in [
        (3_usize, vec![0.04, 0.06, 0.05]),
        (4, vec![0.04, 0.06, 0.05, 0.055]),
    ] {
        let o = run_pyramid_apex_outcome(n, &radii);
        let brep = &o.brep;
        let caps: std::collections::BTreeSet<_> = o.cap_face_handles.iter().copied().collect();
        let (mut rim, mut spoke) = (0.0_f64, 0.0_f64);
        for &face in &o.cap_face_handles {
            let Some(fd) = brep.faces.get(face) else {
                continue;
            };
            let Some(sf) = fd.surface.as_ref() else {
                continue;
            };
            let Some(ld) = brep.loops.get(fd.outer_loop) else {
                continue;
            };
            let start = ld.half_edge;
            let mut he = start;
            for _ in 0..32 {
                let Some(hd) = brep.half_edges.get(he) else {
                    break;
                };
                let neighbour_is_cap = brep
                    .half_edges
                    .get(hd.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| caps.contains(&l.face))
                    .unwrap_or(false);
                if let (Some(pc), Some(curve)) = (
                    hd.pcurve.as_ref(),
                    brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref()),
                ) {
                    let img = |t: f64| {
                        let uv = pc.point_at(t);
                        sf.geometry.point_at(uv.x, uv.y)
                    };
                    let c = |t: f64| (&curve.geometry as &dyn CurveEval).point_at(t);
                    let (lo, hi) = (c(curve.t_start), c(curve.t_end));
                    let fwd = lo.distance_to(&img(pc.t_start)) + hi.distance_to(&img(pc.t_end))
                        <= lo.distance_to(&img(pc.t_end)) + hi.distance_to(&img(pc.t_start));
                    let mut worst = 0.0_f64;
                    for k in 0..=20 {
                        let f = 0.05 + 0.90 * f64::from(k) / 20.0;
                        let s = if fwd { f } else { 1.0 - f };
                        worst = worst.max(
                            img(pc.t_start + s * (pc.t_end - pc.t_start))
                                .distance_to(&c(curve.t_start + f * (curve.t_end - curve.t_start))),
                        );
                    }
                    if neighbour_is_cap {
                        spoke = spoke.max(worst);
                    } else {
                        rim = rim.max(worst);
                    }
                }
                he = hd.next;
                if he == start {
                    break;
                }
            }
        }
        println!(
            "[pk192-split] FINAL/incumbent n={n} cap_rim={rim:.4e} cap_spoke={spoke:.4e} bar=5.000e-3"
        );
    }
}
