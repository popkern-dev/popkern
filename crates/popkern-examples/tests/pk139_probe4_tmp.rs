use popkern_core::geometry::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
fn bulge_rail_curve_forensics() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let out = fillet_edges(
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::ByPoints(vec![(0.0, 0.05), (0.5, 0.12), (1.0, 0.05)]),
        }],
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("bulge completes")
    .brep;
    let strip = out
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface.as_ref().is_some_and(|s| {
                matches!(
                    s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
                )
            })
        })
        .map(|(h, _)| h)
        .unwrap();
    // strip surface stats
    if let Some(s) = out.faces[strip].surface.as_ref()
        && let popkern_core::geometry::surface::SurfaceGeometry::Nurbs(n) = &s.geometry
    {
        eprintln!(
            "strip surface: {}x{} deg {}x{} rational={}",
            n.num_u,
            n.num_v,
            n.degree_u,
            n.degree_v,
            n.weights.iter().any(|w| (w - 1.0).abs() > 1e-12)
        );
    }
    // rail edges = edges 17v1/18v1 pattern: find edges on strip with big trace
    for (eh, e) in out.edges.iter() {
        let Some(c) = e.curve.as_ref() else { continue };
        let he = e.half_edge;
        let on_strip = out.loops.get(out.half_edges[he].loop_).map(|l| l.face) == Some(strip)
            || out
                .half_edges
                .get(out.half_edges[he].twin)
                .and_then(|t| out.loops.get(t.loop_))
                .map(|l| l.face)
                == Some(strip);
        if !on_strip {
            continue;
        }
        let kind = format!("{:?}", std::mem::discriminant(&c.geometry));
        let n_cp = match &c.geometry {
            popkern_core::geometry::CurveGeometry::Nurbs(n) => n.control_points.len(),
            _ => 0,
        };
        let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
        eprintln!(
            "edge {eh:?} kind {kind} n_cp {n_cp} t[{:.3},{:.3}] mid ({:.4},{:.4},{:.4})",
            c.t_start, c.t_end, mid.x, mid.y, mid.z
        );
    }
}
