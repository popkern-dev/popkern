//! PK-143 — identify → modify → unfillet cycle across shipped cap
//! families (stage 1: measurement matrix).
//!
//! A-level round-trip editability: identify any shipped blend from
//! its faces alone, modify its radius in place, and unfillet back to
//! the pre-blend topology. Families measured here; classes that
//! survive become hard gates, gaps become named pins.

#[path = "support/mesh_gauges.rs"]
mod mesh_gauges;

use popkern_blend::{
    FilletOptions, FilletRadius, FilletSpec, fillet_edges, identify_blend, modify_fillet_radius,
    unfillet_face,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

fn cylinder_faces(brep: &BRep) -> Vec<FaceHandle> {
    brep.faces
        .iter()
        .filter(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .map(|(h, _)| h)
        .collect()
}

/// Family 1: single-edge plane-plane fillet (Cylinder band) — the
/// baseline family identify/modify/unfillet were built for. The
/// full cycle is a HARD GATE.
#[test]
fn cylinder_family_round_trips() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let pre_faces = brep.faces.len();
    let edge = brep.edges.keys().next().unwrap();
    let result = fillet_edges(
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        }],
        &FilletOptions::default(),
        default_tol(),
        None,
    )
    .expect("baseline fillet");
    let cyl = *cylinder_faces(&result.brep)
        .first()
        .expect("fillet cylinder face");

    // Identify from the face alone.
    let info = identify_blend(&result.brep, cyl).expect("identify must recognize the fillet");
    let r = info.radius.expect("radius must be extracted");
    assert!((r - 0.1).abs() < 1e-9, "identified radius {r} != 0.1");

    // Modify in place.
    let modified =
        modify_fillet_radius(&result.brep, cyl, 0.15, default_tol()).expect("modify radius");
    let cyl2 = *cylinder_faces(&modified)
        .first()
        .expect("modified cylinder");
    let info2 = identify_blend(&modified, cyl2).expect("re-identify");
    let r2 = info2.radius.expect("radius");
    assert!((r2 - 0.15).abs() < 1e-9, "modified radius {r2} != 0.15");

    // Unfillet back to the sharp cube.
    let restored = unfillet_face(&modified, cyl2).expect("unfillet");
    assert_eq!(
        restored.faces.len(),
        pre_faces,
        "unfillet must restore the pre-blend face census"
    );
    let (structural, open, other) = mesh_gauges::brep_open_defects(&restored);
    let conformity = mesh_gauges::welded_open_edges(&restored);
    eprintln!(
        "[editability] cylinder family: restored defects structural={structural} open={open} \
         other={other}, mesh conformity reads {conformity}"
    );
    assert_eq!(
        (structural, open, other),
        (0, 0, 0),
        "restored body must be structurally closed"
    );
    // PK-157: unfillet heals pcurves on the restored topology, so the
    // restored body must tessellate weld-conforming (was 38 non-welded
    // mesh edges before the heal pass).
    assert_eq!(
        conformity, 0,
        "restored body must mesh weld-conforming; {conformity} non-welded mesh edges"
    );
}

/// Family 2: uniform trihedral corner (SameSide sphere cap). Measures
/// the identify/modify/unfillet behavior on the sphere-cap family.
#[test]
fn sphere_cap_family_cycle_measured() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    // Three edges meeting at one vertex, equal radii.
    let mut chosen = Vec::new();
    'outer: for (v, _) in brep.vertices.iter() {
        chosen.clear();
        for (e, d) in brep.edges.iter() {
            let he = &brep.half_edges[d.half_edge];
            let ends = [he.origin, brep.half_edges[he.twin].origin];
            if ends.contains(&v) {
                chosen.push(e);
            }
        }
        if chosen.len() == 3 {
            break 'outer;
        }
    }
    assert_eq!(chosen.len(), 3, "cube corner has 3 edges");
    let specs: Vec<FilletSpec> = chosen
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(
        &brep,
        &specs,
        &FilletOptions::default(),
        default_tol(),
        None,
    )
    .expect("corner fillet");
    // Find the sphere cap.
    let sphere = result.brep.faces.iter().find(|(_, f)| {
        f.surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }))
    });
    match sphere {
        Some((face, _)) => match identify_blend(&result.brep, face) {
            Some(info) => eprintln!(
                "[editability] sphere cap: identified kind {:?} radius {:?}",
                info.kind, info.radius
            ),
            None => eprintln!("[editability] sphere cap: NOT identified (gap)"),
        },
        None => eprintln!("[editability] sphere cap: fixture produced no Sphere face"),
    }
    // Corner leg cylinders: measure identify coverage.
    let mut identified = 0usize;
    let mut missed = 0usize;
    for cyl in cylinder_faces(&result.brep) {
        if identify_blend(&result.brep, cyl).is_some() {
            identified += 1;
        } else {
            missed += 1;
        }
    }
    eprintln!("[editability] corner legs: identified {identified}, missed {missed}");
}

/// Shell-orientation control: the same single-edge fillet on an
/// OUTWARD-built cube (prism constructor) — discriminates whether
/// the raw-output open reads are an inward-shell-specific class.
#[test]
fn outward_cube_single_edge_fillet_weld_measured() {
    use popkern_core::geometry::{Point3, Vec3};
    let outline = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let brep =
        BRep::prism(&outline, Vec3::new(0.0, 0.0, 1.0), 1.0, default_tol()).expect("outward cube");
    let pristine_opens = mesh_gauges::welded_open_edges(&brep);
    eprintln!("[editability] outward cube pristine opens = {pristine_opens}");
    // PK-149: prism constructor authors planar pcurves, so its raw
    // tessellation conforms (read 358 before the authoring landed).
    assert_eq!(
        pristine_opens, 0,
        "pristine prism tessellation must conform"
    );
    let edge = brep.edges.keys().next().unwrap();
    match fillet_edges(
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        }],
        &FilletOptions::default(),
        default_tol(),
        None,
    ) {
        Ok(result) => {
            let opens = mesh_gauges::welded_open_edges(&result.brep);
            eprintln!("[editability] outward cube fillet opens = {opens}");
            // The fillet pipeline's pcurve authoring heals the raw
            // prism body's non-conforming tessellation (input reads
            // 358) — gate the healed output.
            assert_eq!(opens, 0, "outward-cube fillet output must conform");
        }
        Err(e) => eprintln!("[editability] outward cube fillet DECLINED: {e:?}"),
    }
}

/// Family 3: sphere-cap MODIFY + UNFILLET legs (stage 2
/// measurement). The corner cap's radius modify must either
/// re-dispatch the cap or decline structured; unfillet of the cap
/// face must decline structured (a vertex blend has no single
/// restored edge).
#[test]
fn sphere_cap_modify_and_unfillet_measured() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let mut chosen = Vec::new();
    'outer: for (v, _) in brep.vertices.iter() {
        chosen.clear();
        for (e, d) in brep.edges.iter() {
            let he = &brep.half_edges[d.half_edge];
            let ends = [he.origin, brep.half_edges[he.twin].origin];
            if ends.contains(&v) {
                chosen.push(e);
            }
        }
        if chosen.len() == 3 {
            break 'outer;
        }
    }
    let specs: Vec<FilletSpec> = chosen
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(
        &brep,
        &specs,
        &FilletOptions::default(),
        default_tol(),
        None,
    )
    .expect("corner fillet");
    let sphere = result
        .brep
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }))
        })
        .map(|(h, _)| h)
        .expect("sphere cap");
    // Contract pin (PK-143 stage 2): vertex-blend editing is not
    // implemented; both legs must DECLINE with the structured
    // VertexBlend reason — never panic, never ship a mutated body.
    // Flip these to completion gates when vertex-blend modify /
    // unfillet land.
    let modify_err = modify_fillet_radius(&result.brep, sphere, 0.15, default_tol())
        .expect_err("vertex-blend modify is unimplemented and must decline");
    assert!(
        format!("{modify_err:?}").contains("VertexBlend"),
        "modify decline must name the VertexBlend class: {modify_err:?}"
    );
    let unfillet_err = unfillet_face(&result.brep, sphere)
        .expect_err("vertex-blend unfillet is unimplemented and must decline");
    assert!(
        format!("{unfillet_err:?}").contains("VertexBlend"),
        "unfillet decline must name the VertexBlend class: {unfillet_err:?}"
    );
    eprintln!("[editability] sphere-cap modify+unfillet decline structured (VertexBlend)");
}

/// Family 4: pyramid4 uniform apex (Gregory/partition cap family).
/// Identify coverage over ALL blend faces + modify/unfillet legs on
/// the cap, measured.
#[test]
fn pyramid4_apex_family_measured() {
    let brep = BRep::square_pyramid(2.0, 1.5, default_tol()).unwrap();
    // The four lateral edges meet at the apex.
    let apex = brep
        .vertices
        .iter()
        .find(|(_, v)| (v.position.z - 1.5).abs() < 1e-9)
        .map(|(h, _)| h)
        .expect("apex vertex");
    let mut lateral = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let ends = [he.origin, brep.half_edges[he.twin].origin];
        if ends.contains(&apex) {
            lateral.push(e);
        }
    }
    assert_eq!(lateral.len(), 4, "square pyramid apex has 4 edges");
    let specs: Vec<FilletSpec> = lateral
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = match fillet_edges(
        &brep,
        &specs,
        &FilletOptions::default(),
        default_tol(),
        None,
    ) {
        Ok(r) => r,
        Err(e) => {
            let msg = format!("{e:?}");
            eprintln!(
                "[editability] pyramid4 apex fillet DECLINED: {}",
                &msg[..msg.len().min(240)]
            );
            return;
        }
    };
    let mut identified = 0usize;
    let mut missed = Vec::new();
    let mut blend_faces = 0usize;
    for (fh, f) in result.brep.faces.iter() {
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        // Blend faces: anything that is not one of the pyramid's
        // original planar walls/base (original planes have no
        // tangent boundaries; skip planes entirely for the census).
        if matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        blend_faces += 1;
        if identify_blend(&result.brep, fh).is_some() {
            identified += 1;
        } else {
            missed.push((fh, format!("{:?}", std::mem::discriminant(&sof.geometry))));
        }
    }
    eprintln!(
        "[editability] pyramid4 apex: {blend_faces} non-plane blend faces, {identified} identified, missed {missed:?}"
    );
    assert!(
        blend_faces >= 5,
        "pyramid4 apex fillet should carry 4 leg strips + cap faces"
    );
    assert!(
        missed.is_empty(),
        "every pyramid4 blend face must identify: missed {missed:?}"
    );
}
