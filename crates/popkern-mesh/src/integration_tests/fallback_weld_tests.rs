//! PK-299 pin: a face that DECLINES boundary-first must still weld
//! watertight against served boundary-first neighbours.
//!
//! The uniform-grid tessellator used to sample its trim polylines adaptively from
//! each pcurve — boundary vertices then sat on private polyline
//! chords up to the UV chord-height off the shared edge curve,
//! beyond the crack zipper's matching band, and every shared edge
//! read open. With canonical-station trim sampling the fallback
//! face's boundary is bit-identical to its neighbours' and the
//! zipper closes the clip-introduced T-junctions exactly.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::TessellationOptions;

#[test]
fn pcurve_stripped_wall_face_welds_watertight() {
    let tol = Tolerance::default();
    let brep = {
        let mut b = BRep::cylinder(1.0, 2.0, tol).expect("cylinder");
        // Strip the BOTTOM CAP's rim pcurve: the planar cap declines
        // boundary-first while the wall and top cap still serve, and
        // its rim is a CURVED shared edge — the class where private
        // chord sets open the weld.
        let cap = b
            .faces
            .iter()
            .find(|(fh, fd)| {
                let planar = matches!(
                    fd.surface.as_ref().map(|s| &s.geometry),
                    Some(popkern_core::geometry::SurfaceGeometry::Plane { .. })
                );
                planar && {
                    let start = b.loops[fd.outer_loop].half_edge;
                    let p = b.vertices[b.half_edges[start].origin].position;
                    let _ = fh;
                    p.z.abs() < 1e-9
                }
            })
            .map(|(fh, _)| fh)
            .expect("bottom cap");
        let start = b.loops[b.faces[cap].outer_loop].half_edge;
        b.half_edges[start].pcurve = None;
        b
    };

    // The fixture must actually exercise the fallback: the cap face
    // declines boundary-first for its missing pcurve.
    let wall = brep
        .faces
        .iter()
        .find(|(_, fd)| {
            let planar = matches!(
                fd.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. })
            );
            planar && {
                let start = brep.loops[fd.outer_loop].half_edge;
                brep.vertices[brep.half_edges[start].origin]
                    .position
                    .z
                    .abs()
                    < 1e-9
            }
        })
        .map(|(fh, _)| fh)
        .expect("bottom cap");
    let options = TessellationOptions::default();
    let edge_samples = crate::boundary_first::sample_edges(&brep, &options, None);
    let mut sink = Vec::new();
    assert!(
        crate::boundary_first::tessellate_face_boundary_first(
            &brep,
            wall,
            &edge_samples,
            &options,
            &mut sink,
        )
        .is_none(),
        "fixture must force a boundary-first decline on the cap"
    );

    let mesh = crate::tessellate(&brep, &options).expect("tessellation");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut next = 0usize;
    let mut id = |p: Point3| -> usize {
        let k = (
            (p.x * 1.0e6).round() as i64,
            (p.y * 1.0e6).round() as i64,
            (p.z * 1.0e6).round() as i64,
        );
        *key.entry(k).or_insert_with(|| {
            next += 1;
            next - 1
        })
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let mut opens: Vec<((usize, usize), usize)> = Vec::new();
    for (&k, &c) in &counts {
        if c == 1 {
            opens.push((k, c));
        }
    }
    let mut pos: Vec<Point3> = vec![Point3::new(0.0, 0.0, 0.0); key.len()];
    for t in &mesh.triangles {
        for &vi in t {
            let p = mesh.vertices[vi].position;
            let k = (
                (p.x * 1.0e6).round() as i64,
                (p.y * 1.0e6).round() as i64,
                (p.z * 1.0e6).round() as i64,
            );
            if let Some(&i) = key.get(&k) {
                pos[i] = p;
            }
        }
    }
    for ((a, b), _) in opens.iter().take(12) {
        eprintln!(
            "open ({:.4},{:.4},{:.4})-({:.4},{:.4},{:.4})",
            pos[*a].x, pos[*a].y, pos[*a].z, pos[*b].x, pos[*b].y, pos[*b].z
        );
    }
    let open = counts.values().filter(|&&c| c == 1).count();
    assert_eq!(
        open, 0,
        "fallback cap must weld watertight against served neighbours (PK-299), {open} open"
    );
}
