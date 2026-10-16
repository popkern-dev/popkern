//! A tapped hole a hand-typed axis off its face's normal serves
//! (PK-845).
//!
//! On the angled fitting's wedge the tap's exact axis is the inclined
//! face's normal; a user's typed axis rarely is. The plane x cone
//! intersection snapped to the cone's cross-section circle whenever the
//! tilt passed the angular tolerance, though that circle sits r·sin(a)
//! off the plane: at 1e-4 rad the mouth circle lay 3.3e-4 off the face
//! and the tap declined on curve-surface validation (at 1e-3 and 1e-2
//! it declined with curveless edges). The snap is now bounded by the
//! linear tolerance, as plane x cylinder already was.

use popkern_rhai::ScriptHost;
use std::collections::HashMap;

#[test]
fn a_tapped_hole_tilted_off_the_face_normal_serves() {
    let n_len = (9.0f64 * 9.0 + 24.0 * 24.0).sqrt();
    let (nx, nz) = (-9.0 / n_len, -24.0 / n_len);
    for tilt in [1e-5, 1e-4, 1e-3, 1e-2] {
        let (c, s) = (f64::cos(tilt), f64::sin(tilt));
        let (ax, az) = (nx * c - nz * s, nx * s + nz * c);
        let source = format!(
            "let m = model();\n\
             let b = m.profile(xz()).path(0.0, 0.0).line_to(44.0, 0.0).line_to(44.0, 5.0)\
             .line_to(24.0, 5.0).line_to(0.0, 14.0).close().sketch(\"wedge\").extrude(blind(24.0));\n\
             b = b.tapped_hole(axis_at([12.0, -12.0, 9.5], [{ax}, 0.0, {az}]), \"M6\", blind(9.0), \
             #{{ clearance: 0.2, label: \"tap\" }});\nm"
        );
        let mut run = ScriptHost::new()
            .run(&source, "tilted_tap.rhai")
            .unwrap_or_else(|e| panic!("tilt {tilt}: {}", e.render()));
        let brep = run.model.take().expect("model").to_brep().expect("brep");
        assert!(brep.validate().is_ok(), "tilt {tilt}: invalid");
        let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::display())
            .expect("mesh");
        let key =
            |p: &popkern_core::geometry::Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        let mut weld: HashMap<[u64; 3], usize> = HashMap::new();
        let ids: Vec<usize> = mesh
            .vertices
            .iter()
            .map(|v| {
                let n = weld.len();
                *weld.entry(key(&v.position)).or_insert(n)
            })
            .collect();
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (ids[t[k]], ids[t[(k + 1) % 3]]);
                *edges.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        assert_eq!(
            edges.values().filter(|&&n| n != 2).count(),
            0,
            "tilt {tilt}: watertight"
        );
    }
}
