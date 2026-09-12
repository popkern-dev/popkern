//! PK-480 diagnostic probe (temporary): where are the angled tap's
//! open mesh edges?

use popkern_core::geometry::{Point3, Vec3};
use popkern_part::Part;
use popkern_part::fastener::{HoleDepth, TapOptions, ThreadSpec};
use popkern_part::prelude::*;
use std::collections::HashMap;

#[test]
#[ignore = "PK-480 probe: run explicitly"]
fn locate_open_edges() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let plate: Part = Profile::on(Plane::XY)
        .rect(20.0, 20.0)
        .extrude(8.0)
        .expect("block");
    let direction = Vec3::new(0.35, 0.2, -1.0);
    let tapped = plate
        .tapped_hole(
            &spec,
            Point3::new(1.0, -0.5, 8.0),
            direction,
            HoleDepth::Blind(4.5),
            &TapOptions {
                clearance: 0.15,
                entry_chamfer: true,
            },
        )
        .expect("angled tapped hole");
    let mesh = tapped
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let key = (a.min(b), a.max(b));
            *counts.entry(key).or_insert(0) += 1;
        }
    }
    let o = Point3::new(1.0, -0.5, 8.0);
    let d = direction.normalized().expect("dir");
    let mut open = 0;
    for (&(a, b), &c) in &counts {
        if c == 1 {
            open += 1;
            let pa = mesh.vertices[a].position;
            let pb = mesh.vertices[b].position;
            let hole = |p: Point3| {
                let v = p - o;
                let z = v.dot(&d);
                let r = (v.dot(&v) - z * z).max(0.0).sqrt();
                (z, r)
            };
            let (za, ra) = hole(pa);
            let (zb, rb) = hole(pb);
            eprintln!(
                "[open-edge] ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) hole(z={za:.4},r={ra:.4})->(z={zb:.4},r={rb:.4})",
                pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
            );
        }
    }
    eprintln!("[open-edge] total {open}");
}
