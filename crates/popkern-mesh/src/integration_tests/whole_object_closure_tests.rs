//! PK-095: whole-object tessellation of `fillet_edges` output — the
//! per-face diagnostic probe and (eventually) the closure gate.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn pentagon_blend() -> BRep {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("pipeline");
    result.brep.validate().expect("valid blend");
    result.brep
}

/// Count boundary edges of a face mesh after welding (1e-6).
fn face_mesh_boundary(mesh: &crate::FaceMesh) -> (usize, usize) {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut remap = Vec::new();
    for v in &mesh.vertices {
        let k = (
            (v.position.x * 1.0e6).round() as i64,
            (v.position.y * 1.0e6).round() as i64,
            (v.position.z * 1.0e6).round() as i64,
        );
        let next = key.len();
        remap.push(*key.entry(k).or_insert(next));
    }
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    let mut degen = 0;
    for t in &mesh.triangles {
        let w = [remap[t[0]], remap[t[1]], remap[t[2]]];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            degen += 1;
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    // A correct face mesh is an open disk: its boundary edges should
    // form exactly one loop tracing the trim boundary; count edges
    // with odd counts as a proxy for internal tearing.
    let torn = counts.values().filter(|&&c| c > 2).count();
    (degen, torn)
}

#[test]
#[ignore = "PK-095 diagnostic: per-face tessellation report for the pentagon blend"]
fn per_face_tessellation_report() {
    let brep = pentagon_blend();
    for (fh, face) in brep.faces.iter() {
        let kind = match face.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) => "Gregory",
            Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cylinder",
            _ => "Other",
        };
        let loops = crate::trim::build_trim_loops(
            &brep,
            fh,
            &face.surface.as_ref().unwrap().geometry,
            crate::TessellationOptions::default().chord_height,
        );
        let loop_info: Vec<String> = loops
            .iter()
            .map(|l| {
                format!(
                    "pts={} from_pcurves={} outer={}",
                    l.points.len(),
                    l.from_pcurves,
                    l.is_outer
                )
            })
            .collect();
        match crate::tessellate::tessellate_face(&brep, fh, &crate::TessellationOptions::default())
        {
            Ok(mesh) => {
                let (degen, torn) = face_mesh_boundary(&mesh);
                eprintln!(
                    "face {fh:?} {kind}: tris={} degen={degen} torn_edges={torn} loops=[{}]",
                    mesh.triangles.len(),
                    loop_info.join("; ")
                );
            }
            Err(e) => eprintln!("face {fh:?} {kind}: TESSELLATION FAILED: {e}"),
        }
    }
}
