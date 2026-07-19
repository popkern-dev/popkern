use popkern_core::geometry::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};

#[test]
fn chamfer_end_edge_forensics() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let out = chamfer_edges(
        &cube,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(0.1),
        }],
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("chamfer")
    .brep;
    for (eh, e) in out.edges.iter() {
        let he = e.half_edge;
        let d = &out.half_edges[he];
        let o = out.vertices[d.origin].position;
        let t = out.half_edges[d.twin].origin;
        let tp = out.vertices[t].position;
        let curve = e
            .curve
            .as_ref()
            .map(|c| {
                let a = c.geometry.point_at(c.t_start);
                let b = c.geometry.point_at(c.t_end);
                format!(
                    "t[{:.3},{:.3}] ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    c.t_start, c.t_end, a.x, a.y, a.z, b.x, b.y, b.z
                )
            })
            .unwrap_or("NO CURVE".into());
        eprintln!(
            "{eh:?}: v ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) curve {curve}",
            o.x, o.y, o.z, tp.x, tp.y, tp.z
        );
    }
}
