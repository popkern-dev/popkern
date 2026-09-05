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

/// Spring end grinding: the plate's section loop crosses the helical
/// tube wall's seam mid-way — the march wraps the geometrically-closed
/// chart, the record chain closes through the seam and rim junctions,
/// and the subtraction serves a valid body. (Mesh watertightness of
/// the cut helical wall is PK-441's remaining frontier — this pins the
/// BOOLEAN serve, which OCC validates.)
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
}
