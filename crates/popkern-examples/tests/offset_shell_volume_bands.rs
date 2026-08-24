//! Volume-band acceptance for exact curved offset/shell serving
//! (PK-322): shelled and offset cone/cylinder fixtures must serve
//! validate_all-clean with tessellated volume inside the analytic band.
//! The volume oracle exists precisely because inverted or misplaced
//! skins can pass every structural validator (the
//! wrong-shape-served-clean class) — position and volume measurements
//! are the honest gates.

use std::f64::consts::PI;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellOptions, offset::offset_solid, shell, thicken::thicken};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Tessellated absolute volume via the divergence theorem.
fn mesh_volume(brep: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

fn assert_volume_band(measured: f64, expected: f64, rel_band: f64, label: &str) {
    let rel = (measured - expected).abs() / expected;
    assert!(
        rel <= rel_band,
        "{label}: volume {measured:.6} vs analytic {expected:.6} \
         (relative error {rel:.4}, band {rel_band})"
    );
}

fn find_plane_face(brep: &BRep) -> popkern_core::topology::FaceHandle {
    use popkern_core::geometry::surface::SurfaceGeometry;
    brep.faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .map(|(h, _)| h)
        .unwrap()
}

#[test]
fn shelled_cone_volume_in_analytic_band() {
    // shell(cone r=1 h=2, inward 0.1, base removed): the wall solid is
    // the outer cone minus the inner cone (apex translated down the
    // axis by d/sin(half_angle), rim flush in the base plane).
    let (r, h, d) = (1.0_f64, 2.0_f64, 0.1_f64);
    let brep = BRep::cone(r, h, tol()).unwrap();
    let base = find_plane_face(&brep);
    let result = shell(
        &brep,
        d,
        ShellDirection::Inward,
        &[base],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    assert!(
        result.brep.validate_all().is_empty(),
        "{:?}",
        result.brep.validate_all()
    );

    let half_angle = (r / h).atan();
    let h_inner = h - d / half_angle.sin();
    let r_inner = h_inner * half_angle.tan();
    let expected = PI / 3.0 * (r * r * h - r_inner * r_inner * h_inner);
    assert_volume_band(mesh_volume(&result.brep), expected, 0.02, "shelled cone");
}

#[test]
fn offset_cone_volume_in_analytic_band() {
    // offset_solid(cone, +0.15): every face moves out 0.15 with sharp
    // (intersection-join) edges — the grown body is the cone whose apex
    // moved up d/sin(half_angle) and whose base plane moved down d.
    let (r, h, d) = (1.0_f64, 2.0_f64, 0.15_f64);
    let brep = BRep::cone(r, h, tol()).unwrap();
    let result = offset_solid(&brep, d, &ShellOptions::default(), tol()).unwrap();
    assert!(
        result.brep.validate_all().is_empty(),
        "{:?}",
        result.brep.validate_all()
    );

    let half_angle = (r / h).atan();
    let h_grown = h + d / half_angle.sin() + d;
    let r_grown = h_grown * half_angle.tan();
    let expected = PI / 3.0 * r_grown * r_grown * h_grown;
    assert_volume_band(mesh_volume(&result.brep), expected, 0.02, "offset cone");

    // The apex vertex must sit exactly at the translated apex.
    let apex_z = h + d / half_angle.sin();
    let apex_count = result
        .brep
        .vertices
        .values()
        .filter(|v| {
            v.position.x.abs() < 1e-9
                && v.position.y.abs() < 1e-9
                && (v.position.z - apex_z).abs() < 1e-9
        })
        .count();
    assert_eq!(apex_count, 1, "apex translated to z={apex_z}");
}

#[test]
fn thickened_cylinder_volume_in_analytic_band() {
    // thicken(cylinder r=1 h=2, 0.2): hollow wall between the outward
    // and inward half-offset skins.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
    assert!(
        result.brep.validate_all().is_empty(),
        "{:?}",
        result.brep.validate_all()
    );
    let expected = PI * (1.1 * 1.1 * 2.2 - 0.9 * 0.9 * 1.8);
    assert_volume_band(
        mesh_volume(&result.brep),
        expected,
        0.02,
        "thickened cylinder",
    );
}

#[test]
fn volume_gauge_matches_primitives_and_hollow_bodies() {
    // Calibration of the gauge itself, including the HOLLOW case: the
    // mesher's global winding pass used to force every mesh component
    // to enclose positive volume, so a hollow body's void skin was
    // deliberately flipped inside-out and its volume ADDED instead of
    // subtracting (thickened cube measured 2.240 = outer+inner). The
    // containment-parity side rule makes void skins subtract.
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let cone = BRep::cone(1.0, 2.0, tol()).unwrap();
    let cube = BRep::unit_cube(tol()).unwrap();
    assert_volume_band(mesh_volume(&cyl), PI * 2.0, 0.01, "solid cylinder");
    // 1.5% band doubles as the seam-crack pin: with the PK-332 grid
    // truncation the cone read -2.02% (crack + chord); chord error
    // alone is ~0.7%.
    assert_volume_band(mesh_volume(&cone), PI * 2.0 / 3.0, 0.015, "solid cone");
    assert_volume_band(mesh_volume(&cube), 1.0, 1e-9, "unit cube");
    let tc = thicken(&cube, 0.2, &ShellOptions::default(), tol()).unwrap();
    assert_volume_band(
        mesh_volume(&tc.brep),
        1.2_f64.powi(3) - 0.8_f64.powi(3),
        1e-9,
        "thickened cube wall",
    );
}

/// Open (once-used) mesh edges after quantized position welding (1e-7:
/// seam pairs evaluate at u=0 vs u=2π and land ULPs apart) — a
/// watertight tessellation has none.
fn open_mesh_edges(brep: &BRep) -> usize {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let q = |x: f64| (x * 1.0e7).round() as i64;
    let key = |i: usize| {
        let p = mesh.vertices[i].position;
        (q(p.x), q(p.y), q(p.z))
    };
    let mut counts = std::collections::HashMap::new();
    for t in &mesh.triangles {
        let w = [key(t[0]), key(t[1]), key(t[2])];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0_i32) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}

#[test]
fn primitives_tessellate_watertight() {
    // PK-333: the adaptive quadtree emitted hanging nodes at every
    // refinement transition (an 80-edge open ring at the cone's
    // near-apex refinement boundary), and the 2% grid pad crossed the
    // cone's natural v ≥ 0 bound into the mirror nappe (an 8-edge
    // beyond-apex micro-ring). Conforming closure + natural-bound
    // clamping make the meshes genuinely watertight; hollow bodies
    // stay watertight per shell.
    let cube = BRep::unit_cube(tol()).unwrap();
    for (name, brep) in [
        ("cone", BRep::cone(1.0, 2.0, tol()).unwrap()),
        ("cylinder", BRep::cylinder(1.0, 2.0, tol()).unwrap()),
        ("cube", cube.clone()),
        (
            "thickened cube",
            thicken(&cube, 0.2, &ShellOptions::default(), tol())
                .unwrap()
                .brep,
        ),
    ] {
        let open = open_mesh_edges(&brep);
        assert_eq!(open, 0, "{name}: {open} open mesh edges");
    }
}

#[test]
fn shelled_bodies_tessellate_watertight() {
    // PK-334: shelled bodies mix tessellation paths (boundary-first
    // walls, legacy-clip rim faces), whose boundary chains carry
    // DIFFERENT station sets on the shared rim edges — all on the
    // same curve after sewing. The zipper's curve-parameter station
    // merge closes them; before it, the shelled cylinder read 27 open
    // edges on the opening rim (detection evidence for the pin).
    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let cone = BRep::cone(1.0, 2.0, tol()).unwrap();
    let shelled = |brep: &BRep, z: f64| {
        let top = brep
            .faces
            .keys()
            .find(|&fh| {
                use popkern_core::geometry::surface::SurfaceGeometry;
                matches!(
                    brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
                    Some(SurfaceGeometry::Plane { .. })
                ) && {
                    let ps = brep.face_outer_loop_positions(fh);
                    !ps.is_empty() && ps.iter().all(|p| (p.z - z).abs() < 1e-9)
                }
            })
            .unwrap();
        shell(
            brep,
            0.1,
            ShellDirection::Inward,
            &[top],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .unwrap()
        .brep
    };
    // (name, body, welds must close BY CONSTRUCTION — zipper unused)
    // PK-335: slit dissolution lets the rim wall serve boundary-first,
    // so the cube/cylinder shells weld with zero zipper repairs. Cone
    // bodies keep zipper reliance: their full-u-wrap laterals decline
    // boundary-first (PK-336) and the legacy grid staggers against
    // canonical neighbours — welded by the PK-334 parameter merge.
    for (name, brep, construction_weld) in [
        ("shelled cube", shelled(&cube, 1.0), true),
        ("shelled cylinder", shelled(&cyl, 2.0), true),
        ("shelled cone (base removed)", shelled(&cone, 0.0), false),
    ] {
        let open = open_mesh_edges(&brep);
        assert_eq!(open, 0, "{name}: {open} open mesh edges");
        if construction_weld {
            let mesh =
                popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default())
                    .expect("tessellate");
            assert_eq!(
                mesh.zipped_t_junctions, 0,
                "{name}: welds must close by construction, not zipper repair"
            );
        }
    }
}
