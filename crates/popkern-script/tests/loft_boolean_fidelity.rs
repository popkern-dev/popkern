//! PK-443: chart-wrapping boolean loops on closed lofted and swept
//! walls, and the loft construction fidelity they exposed.
//!
//! The propeller-blade fixture found two loft construction lies and a
//! whole boolean substrate gap:
//! - rational sections lofted with weights STAMPED to 1.0 — every arc
//!   span became a parabolic sheet 6.1% of the arc radius off, while
//!   the cap rim edges carried the exact arcs;
//! - span-wise lofting parameterized each span family by its own
//!   centroid chord-length, so adjacent spans disagreed along their
//!   shared boundary between sections (9.1e-5 at the girdle junctions);
//! - a chart-wrapping intersection loop on the closed wall walked as
//!   an interior ISLAND and the union dropped or mangled the pieces.

use popkern_script::Plane;
use popkern_script::body_sweep::HelixParams;
use popkern_script::prelude::*;

fn blade_section(m: &mut Model, z: f64, c: f64, t: f64, pitch_deg: f64, label: &str) -> Body {
    let pitch = pitch_deg * std::f64::consts::PI / 180.0;
    let (cs, sn) = (pitch.cos(), pitch.sin());
    let rot = |x: f64, y: f64| (x * cs - y * sn, x * sn + y * cs);
    let t2 = 0.5 * t;
    let c2 = 0.5 * c - t2;
    let p0 = rot(-c2, -t2);
    let p1 = rot(c2, -t2);
    let v1 = rot(c2 + t2, 0.0);
    let p2 = rot(c2, t2);
    let p3 = rot(-c2, t2);
    let v2 = rot(-c2 - t2, 0.0);
    m.profile(Plane::xy_at(z))
        .path(p0.0, p0.1)
        .line_to(p1.0, p1.1)
        .arc_to(v1, p2)
        .line_to(p3.0, p3.1)
        .arc_to(v2, p0)
        .close()
        .sketch(label)
        .expect("section sketch")
}

fn blade(m: &mut Model, tag: &str) -> Body {
    let s0 = blade_section(m, 0.08, 0.16, 0.050, 55.0, &format!("{tag}_root"));
    let s1 = blade_section(m, 0.22, 0.22, 0.032, 42.0, &format!("{tag}_mid1"));
    let s2 = blade_section(m, 0.38, 0.19, 0.020, 32.0, &format!("{tag}_mid2"));
    let s3 = blade_section(m, 0.52, 0.09, 0.010, 25.0, &format!("{tag}_tip"));
    m.loft(&[s0, s1, s2, s3]).expect("loft")
}

/// A stadium loft (rational arc spans, twisted sections) must carry NO
/// curve-surface lies: every edge curve lies on its adjacent surfaces.
/// This is the loft-fidelity gate — weight-dropping or per-span
/// station parameterization both fail it from birth.
#[test]
fn stadium_loft_carries_no_curve_surface_lies() {
    let mut m = Model::new();
    let blade = blade(&mut m, "blade");
    assert!(blade.is_alive(&m));
    let brep = m.brep_of(blade).expect("brep");
    let errors = brep.validate_curve_surface_all();
    assert!(
        errors.is_empty(),
        "lofted blade carries curve-surface lies: {errors:?}"
    );
}

/// Blade ∪ hub: the girdle (a chart-wrapping loop on the blade's
/// closed lofted wall) splits the wall into bands, the union serves,
/// every face tessellates, and the mesh is wound outward (the model is
/// star-shaped about the origin, so every face's signed volume is
/// positive).
#[test]
fn blade_hub_union_serves_with_full_outward_mesh() {
    let mut m = Model::new();
    let hub = m
        .profile(Plane::XZ)
        .path(0.0, -0.22)
        .line_to(0.14, -0.22)
        .line_to(0.14, 0.10)
        .line_to(0.025, 0.30)
        .line_to(0.0, 0.30)
        .close()
        .sketch("hub_section")
        .expect("hub sketch")
        .revolve(
            &mut m,
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            RevolveEnd::FullCircle,
            "hub",
        )
        .expect("hub revolve");
    let blade = blade(&mut m, "blade");
    let x_axis = AxisDesc::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [1.0, 0.0, 0.0],
    };
    let blade = blade.rotate(&mut m, x_axis, -90.0).expect("rotate");
    let joined = hub.union(&mut m, blade).expect("union");
    let brep = m.brep_of(joined).expect("brep");
    assert!(
        brep.faces.len() > 6,
        "hub plus blade faces expected, got {}",
        brep.faces.len()
    );

    let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut face_vol: std::collections::HashMap<popkern_core::topology::FaceHandle, f64> =
        std::collections::HashMap::new();
    for (ti, tri) in mesh.triangles.iter().enumerate() {
        let [a, b, c] = *tri;
        let (pa, pb, pc) = (
            mesh.vertices[a].position,
            mesh.vertices[b].position,
            mesh.vertices[c].position,
        );
        let cross = (
            pb.y * pc.z - pb.z * pc.y,
            pb.z * pc.x - pb.x * pc.z,
            pb.x * pc.y - pb.y * pc.x,
        );
        *face_vol.entry(mesh.face_map[ti]).or_default() +=
            (pa.x * cross.0 + pa.y * cross.1 + pa.z * cross.2) / 6.0;
    }
    assert_eq!(face_vol.len(), brep.faces.len(), "every face tessellates");
    for (fh, v) in &face_vol {
        assert!(*v > 0.0, "face {fh:?} wound inward (signed vol {v:+.6})");
    }
}

/// A lofted 16-gon-section blade unioned with the revolved hub
/// (PK-442): every facet is its own face today, the facet × hub
/// records march and imprint cleanly, and the union serves
/// watertight. (The single-creased-wall variant of this class lives
/// at core level: `surface_surface_crease_tests` pins the march
/// through C0 creases with geometric-crease-gated splitting.)
#[test]
fn polygon_section_blade_union_hub_serves_watertight() {
    fn polygon_section(m: &mut Model, z: f64, r: f64, label: &str) -> Body {
        let n = 16usize;
        let p = m.profile(Plane::xy_at(z));
        let mut path = p.path(r, 0.0);
        for k in 1..n {
            let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            path = path.line_to(r * a.cos(), r * a.sin());
        }
        path.close().sketch(label).expect("section")
    }
    let mut m = Model::new();
    let hub = m
        .profile(Plane::XZ)
        .path(0.0, -0.22)
        .line_to(0.14, -0.22)
        .line_to(0.14, 0.10)
        .line_to(0.025, 0.30)
        .line_to(0.0, 0.30)
        .close()
        .sketch("hub_section")
        .expect("hub sketch")
        .revolve(
            &mut m,
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 0.0, 1.0],
            },
            RevolveEnd::FullCircle,
            "hub",
        )
        .expect("hub revolve");
    let s0 = polygon_section(&mut m, 0.08, 0.055, "b0");
    let s1 = polygon_section(&mut m, 0.25, 0.075, "b1");
    let s2 = polygon_section(&mut m, 0.45, 0.045, "b2");
    let blade = m.loft(&[s0, s1, s2]).expect("loft");
    let x_axis = AxisDesc::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [1.0, 0.0, 0.0],
    };
    let blade = blade.rotate(&mut m, x_axis, -90.0).expect("rotate");
    let joined = hub.union(&mut m, blade).expect("union");
    let brep = m.brep_of(joined).expect("brep");
    assert!(
        brep.faces.len() > 10,
        "hub plus faceted blade expected, got {} faces",
        brep.faces.len()
    );
    assert!(brep.validate().is_ok(), "{:?}", brep.validate().err());
    m.assert_watertight()
        .expect("polygon blade ∪ hub watertight mesh");
}

/// PK-445 mesh serve gate: watertight (position-welded open-edge
/// census — the official topological gate alone passes zero-volume
/// ribbons) AND faithful (signed volume within 2% of the analytic
/// helical-tube volume π·r²·L; the ribbon read 49%, the uniform-grid
/// fallback 9% low with T-junction cracks).
fn assert_meshes_watertight_and_faithful(
    brep: &popkern_core::topology::BRep,
    expected_volume: f64,
    label: &str,
) {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 5.0e-4,
        min_subdivisions: 4,
        max_subdivisions: 512,
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(brep, &options).expect("mesh");
    let mut vol6 = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let b = mesh.vertices[t[1]].position;
        let c = mesh.vertices[t[2]].position;
        vol6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    let volume = vol6 / 6.0;
    let mut edge_uses: std::collections::HashMap<(u64, u64), i64> =
        std::collections::HashMap::new();
    let key = |p: popkern_core::geometry::Point3| -> u64 {
        let q = |x: f64| ((x * 1.0e7).round() as i64) as u64;
        q(p.x)
            .wrapping_mul(0x9E3779B97F4A7C15)
            .wrapping_add(q(p.y).wrapping_mul(0xC2B2AE3D27D4EB4F))
            .wrapping_add(q(p.z))
    };
    for t in &mesh.triangles {
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let a = key(mesh.vertices[t[i]].position);
            let b = key(mesh.vertices[t[j]].position);
            *edge_uses.entry((a.min(b), a.max(b))).or_insert(0) += if a < b { 1 } else { -1 };
        }
    }
    let open: usize = edge_uses.values().filter(|&&c| c != 0).count();
    assert_eq!(open, 0, "{label}: mesh has {open} open edges");
    let rel = (volume - expected_volume).abs() / expected_volume;
    assert!(
        rel < 0.02,
        "{label}: mesh volume {volume:.6} vs analytic {expected_volume:.6} ({:.2}% off)",
        rel * 100.0
    );
}

/// Analytic volume of the 6-turn coil after `cuts` end grinds:
/// π·r²·L (helix arc length L = turns·√((2πR)² + pitch²)) minus one
/// ground sliver per cut. Each grind plane sits 0.005 past the wire
/// centre's extreme, so the sliver is the tube material beyond it —
/// integrated as circular-segment cross-sections along the helix
/// (the ~4.6° pitch tilt of the cross-section plane changes the
/// integral at the 0.3% level of the sliver itself, far inside the
/// mesh gate's band).
fn coil_analytic_volume(cuts: usize) -> f64 {
    let big_r = 0.3_f64;
    let pitch = 0.15_f64;
    let wire_r = 0.04_f64;
    let speed = ((std::f64::consts::TAU * big_r).powi(2) + pitch * pitch).sqrt();
    let full = std::f64::consts::PI * wire_r * wire_r * 6.0 * speed;
    // Sliver: wire centre height z = pitch·t (t in turns); the plane
    // at −0.005 removes the segment below it while z + 0.005 < r.
    let seg = |d: f64| -> f64 {
        // Area of the disc part at signed centre distance d beyond
        // the cut line (d ≤ r).
        wire_r * wire_r * (d / wire_r).acos() - d * (wire_r * wire_r - d * d).sqrt()
    };
    let t_end = (wire_r - 0.005) / pitch;
    let n = 4096;
    let mut sliver = 0.0;
    for k in 0..n {
        let t = t_end * (k as f64 + 0.5) / n as f64;
        sliver += seg(pitch * t + 0.005) * speed * t_end / n as f64;
    }
    full - cuts as f64 * sliver
}

/// Spring end grinding: the plate's section loop crosses the helical
/// tube wall's seam mid-way — the march wraps the geometrically-closed
/// chart, the record chain closes through the seam and rim junctions,
/// and the subtraction serves a valid body (PK-441: chart-continuous
/// trim loops, exact corner records, robust classification). The
/// PK-445 mesh gate below serves at fidelity: the cut wall meshes
/// boundary-first on a per-direction curvature-adaptive lattice
/// (the chart-diagonal probe read the coiled chart by phase luck and
/// left the interior empty — the ribbon class), with breadth-first
/// canonical edge sampling (the depth-first cap left half the seam
/// as one segment whose protecting disk evacuated half the lattice).
#[test]
fn spring_end_grind_subtraction_serves() {
    let mut m = Model::new();
    let z_axis = AxisSpec::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    };
    let spring = m
        .profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("wire sketch")
        .helical_sweep(&mut m, z_axis, HelixParams::new(0.3, 0.15, 6.0))
        .expect("coil");
    let plate = m
        .profile(Plane::xy_at(-0.105))
        .rect(1.0, 1.0)
        .sketch("grind_plate")
        .expect("plate sketch")
        .extrude(&mut m, End::Blind(Scalar::from(0.1)))
        .expect("plate");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let ground = spring.subtract(&mut m, plate).expect("grind");
    let brep = m.brep_of(ground).expect("brep");
    assert!(
        brep.validate().is_ok(),
        "ground spring invalid: {:?}",
        brep.validate().err()
    );
    assert!(
        brep.validate_curve_surface().is_ok(),
        "ground spring curve-surface lie: {:?}",
        brep.validate_curve_surface().err()
    );
    // PK-444: the grind serves through the PRIMARY cellular route —
    // the graph-merge fallback registry must stay silent.
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        fallbacks.is_empty(),
        "grind fell back to graph-merge: {fallbacks:?}"
    );
    // PK-445: the ground coil meshes watertight AND faithful.
    assert_meshes_watertight_and_faithful(&brep, coil_analytic_volume(1), "single grind");
}

/// Double-ended spring grinding: the SECOND cut's plane grazes the
/// coil's top rim within 3.4e-5 of the seam-rim vertex — a triple
/// degeneracy (cut plane × rim circle × seam) that fragments the
/// records into a 5.6e-4 corner micro-circuit. Serving it exercises
/// the whole PK-441 chain: consolidation extent-guarded welds, closed-
/// edge split direction resolution, pcurve preservation through graph
/// splits, knot-span ray seeding, medial fragment sampling, and the
/// output loop's chart-continuity normalization. (Watertight meshing
/// of the cut wall is the chart-aware meshing campaign's scope — see
/// the single-grind pin's note.)
#[test]
fn spring_double_end_grind_serves() {
    let mut m = Model::new();
    let z_axis = AxisSpec::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    };
    let spring = m
        .profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("wire sketch")
        .helical_sweep(&mut m, z_axis, HelixParams::new(0.3, 0.15, 6.0))
        .expect("coil");
    let plate_bot = m
        .profile(Plane::xy_at(-0.105))
        .rect(1.0, 1.0)
        .sketch("grind_bot")
        .expect("plate sketch")
        .extrude(&mut m, End::Blind(Scalar::from(0.1)))
        .expect("plate");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let ground = spring.subtract(&mut m, plate_bot).expect("bottom grind");
    let plate_top = m
        .profile(Plane::xy_at(0.905))
        .rect(1.0, 1.0)
        .sketch("grind_top")
        .expect("plate sketch")
        .extrude(&mut m, End::Blind(Scalar::from(0.1)))
        .expect("plate");
    let ground2 = ground.subtract(&mut m, plate_top).expect("top grind");
    let brep = m.brep_of(ground2).expect("brep");
    assert!(
        brep.validate().is_ok(),
        "double-ground spring invalid: {:?}",
        brep.validate().err()
    );
    assert!(
        brep.validate_curve_surface().is_ok(),
        "double-ground spring curve-surface lie: {:?}",
        brep.validate_curve_surface().err()
    );
    // PK-444: BOTH grinds serve through the PRIMARY cellular route —
    // the graph-merge fallback registry must stay silent through the
    // full double-cut chain.
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        fallbacks.is_empty(),
        "a grind fell back to graph-merge: {fallbacks:?}"
    );
    // PK-445: the double-ground coil meshes watertight AND faithful.
    assert_meshes_watertight_and_faithful(&brep, coil_analytic_volume(2), "double grind");
}

/// Analytic volume of the 6.5-turn coil after a DEEP cut: everything
/// above `plane_z` removed (a band of partial circular-segment
/// sections, then whole turns). Same tilt caveat as
/// [`coil_analytic_volume`] — the ~4.6° section tilt moves the
/// integral well inside the mesh gate's 2% band.
fn coil_deep_cut_volume(plane_z: f64) -> f64 {
    let big_r = 0.3_f64;
    let pitch = 0.15_f64;
    let wire_r = 0.04_f64;
    let turns = 6.5_f64;
    let speed = ((std::f64::consts::TAU * big_r).powi(2) + pitch * pitch).sqrt();
    let full = std::f64::consts::PI * wire_r * wire_r * turns * speed;
    let seg = |d: f64| -> f64 {
        wire_r * wire_r * (d / wire_r).acos() - d * (wire_r * wire_r - d * d).sqrt()
    };
    let s_lo = (plane_z - wire_r) / pitch;
    let s_hi = ((plane_z + wire_r) / pitch).min(turns);
    let n = 4096;
    let mut removed = 0.0;
    for k in 0..n {
        let s = s_lo + (s_hi - s_lo) * (k as f64 + 0.5) / n as f64;
        removed += seg(plane_z - pitch * s) * speed * (s_hi - s_lo) / n as f64;
    }
    removed += (turns - s_hi).max(0.0) * speed * std::f64::consts::PI * wire_r * wire_r;
    full - removed
}

/// STATION-ALIGNED deep grind (PK-448 → PK-449): a cut plane through
/// the helix centreline exactly at a construction station (z = 0.9 =
/// 6 × pitch on the 6.5-turn coil) used to lose the seam × plate
/// pierce point to the ADJACENT COIL TURN — the graph splitter's
/// fixed-count coarse scan ranked the wrong turn's basin, the pierce
/// classified OnFaceInterior, the seam never split at the crossing,
/// the chart walker declined its substrate, and the wall dropped to
/// a knife-edge sample (two phantom face/loop pairs, a curveless
/// self-twinned stub, an OCC-rejected STEP; an off-station cut at
/// z = 0.91 always served by phase luck). With the knot-seeded curve
/// inversion the class SERVES census-clean through primary cellular;
/// the PK-448 unrepairable-stub decline stays armed as the fallback
/// boundary behind it.
#[test]
fn station_aligned_deep_grind_serves_census_clean() {
    let mut m = Model::new();
    let z_axis = AxisSpec::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    };
    let spring = m
        .profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("wire sketch")
        .helical_sweep(&mut m, z_axis, HelixParams::new(0.3, 0.15, 6.5))
        .expect("coil");
    let top = m
        .profile(Plane::xy_at(0.9))
        .rect(1.0, 1.0)
        .sketch("grind_top")
        .expect("plate sketch")
        .extrude(&mut m, End::Blind(Scalar::from(0.2)))
        .expect("plate");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let ground = spring
        .subtract(&mut m, top)
        .expect("station-aligned deep grind serves");
    let brep = m.brep_of(ground).expect("brep");
    assert!(
        brep.validate().is_ok(),
        "station-aligned grind invalid: {:?}",
        brep.validate().err()
    );
    assert!(
        brep.validate_curve_surface().is_ok(),
        "station-aligned grind curve-surface lie: {:?}",
        brep.validate_curve_surface().err()
    );
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        fallbacks.is_empty(),
        "station-aligned grind fell back to graph-merge: {fallbacks:?}"
    );
    // Census-clean: every edge exactly twice-used, no self-twins —
    // the serve contract PK-448 specified.
    for (eh, _) in brep.edges.iter() {
        let uses = brep
            .half_edges
            .iter()
            .filter(|(_, he)| he.edge == eh && brep.loops.contains_key(he.loop_))
            .count();
        assert_eq!(uses, 2, "edge {eh:?} not exactly twice-used");
    }
    assert!(
        brep.half_edges.iter().all(|(h, he)| he.twin != h),
        "served body carries a self-twinned half-edge"
    );
    assert_meshes_watertight_and_faithful(&brep, coil_deep_cut_volume(0.9), "station grind");
}
