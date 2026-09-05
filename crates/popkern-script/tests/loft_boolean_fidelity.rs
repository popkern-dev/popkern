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

/// Spring end grinding: the plate's section loop crosses the helical
/// tube wall's seam mid-way — the march wraps the geometrically-closed
/// chart, the record chain closes through the seam and rim junctions,
/// and the subtraction serves a valid body (PK-441: chart-continuous
/// trim loops, exact corner records, robust classification). Mesh
/// watertightness at full FIDELITY is the chart-aware meshing
/// campaign's frontier: boundary-first honestly declines the cut
/// wall (its boundary-only Delaunay is a ribbon of tube-diameter
/// chords) and the uniform-grid fallback leaves rim T-junction
/// chains beyond the crack zipper's band.
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
}
