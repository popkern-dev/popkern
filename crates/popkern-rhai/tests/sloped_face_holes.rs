//! Countersunk holes entering or leaving a sloped face serve wherever
//! they sit on it (PK-846).
//!
//! On a wedge whose top slopes z = 14 - x/3, an M5 countersunk
//! clearance hole entering the slope declined at x = 6.5, 8 and 10
//! ("Loop not geometrically closed: gap 2.19", the countersink cone's
//! closed ellipse on the slope came back open), and one entering the
//! flat underside declined once it exited through the slope. Both
//! families now serve across the slope; each body validates, meshes
//! watertight, and encloses what its quadrature and its mesh agree on,
//! with at least the bore's own volume removed.

use popkern_rhai::ScriptHost;
use std::collections::HashMap;

const WEDGE_VOLUME: f64 = 240.0 * 16.0;

fn holed(x: f64, from_slope: bool) -> popkern_core::topology::BRep {
    let z = 14.0 - x / 3.0;
    let (at, dir) = if from_slope {
        (format!("[{x}, -8.0, {z}]"), "[0.0, 0.0, -1.0]")
    } else {
        (format!("[{x}, -8.0, 0.0]"), "[0.0, 0.0, 1.0]")
    };
    let source = format!(
        "let m = model();\n\
         let b = m.profile(xz()).path(0.0, 0.0).line_to(24.0, 0.0).line_to(24.0, 6.0)\
         .line_to(0.0, 14.0).close().sketch(\"wedge\").extrude(blind(16.0));\n\
         b = b.clearance_hole(axis_at({at}, {dir}), \"M5\", through_all(), \
         #{{ countersink: 0.2, label: \"hole\" }});\nm"
    );
    let mut run = ScriptHost::new()
        .run(&source, "sloped_hole.rhai")
        .unwrap_or_else(|e| panic!("x {x} from_slope {from_slope}: {}", e.render()));
    run.model.take().expect("model").to_brep().expect("brep")
}

#[test]
fn countersunk_holes_serve_across_a_sloped_face() {
    for from_slope in [true, false] {
        for x in [4.0, 5.0, 6.5, 8.0, 10.0, 14.0, 20.0] {
            let label = format!("x {x} from_slope {from_slope}");
            let brep = holed(x, from_slope);
            let issues = brep.validate_all();
            assert!(issues.is_empty(), "{label}: {issues:?}");
            let quadrature = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
            let mesh = popkern_mesh::tessellate(
                &brep,
                &popkern_mesh::TessellationOptions {
                    chord_height: 0.002,
                    ..Default::default()
                },
            )
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
            let mut meshed = 0.0;
            for t in &mesh.triangles {
                let p = t.map(|i| mesh.vertices[i].position);
                meshed += (p[0].x * (p[1].y * p[2].z - p[2].y * p[1].z)
                    - p[0].y * (p[1].x * p[2].z - p[2].x * p[1].z)
                    + p[0].z * (p[1].x * p[2].y - p[2].x * p[1].y))
                    / 6.0;
                for k in 0..3 {
                    let (a, b) = (ids[t[k]], ids[t[(k + 1) % 3]]);
                    *edges.entry((a.min(b), a.max(b))).or_default() += 1;
                }
            }
            assert_eq!(
                edges.values().filter(|&&n| n != 2).count(),
                0,
                "{label}: watertight"
            );
            assert!(
                ((meshed - quadrature) / quadrature).abs() < 1e-3,
                "{label}: {meshed} vs {quadrature}"
            );
            // At least the bore (M5 clearance, diameter >= 5.3) through
            // the wedge's height at the axis.
            let bore = std::f64::consts::PI * 2.65 * 2.65 * (14.0 - x / 3.0);
            assert!(
                WEDGE_VOLUME - quadrature > bore,
                "{label}: removed {}",
                WEDGE_VOLUME - quadrature
            );
        }
    }
}
