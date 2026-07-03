//! PK-095 follow-up: base-end anatomy of the G2 strips — the "indent"
//! near the bottom of the fillet edge and the floor ledge.
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-095 diagnostic: strip base-end anatomy"]
fn strip_base_end_anatomy() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex || t == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    let options = popkern_fillet::FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G2,
        ..Default::default()
    };
    let result =
        popkern_fillet::fillet_edges(&pyramid, &specs, &options, tol, None).expect("blend");
    let brep = &result.brep;

    for (fh, face) in brep.faces.iter() {
        let Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) =
            face.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        let surface = &face.surface.as_ref().unwrap().geometry;
        // Strip end borders: u=0 and u=1 iso curves — z profile.
        for uiso in [0.0, 1.0] {
            let zs: Vec<f64> = (0..=8)
                .map(|i| surface.point_at(uiso, i as f64 / 8.0).z)
                .collect();
            let s: Vec<String> = zs.iter().map(|z| format!("{z:.4}")).collect();
            eprintln!("strip {fh:?} u={uiso}: z = [{}]", s.join(", "));
        }
        // Trim straightness near the base: the plane-plane contact
        // line is straight; measure the strip border's deviation from
        // the straight line through its ends, per v-side.
        for viso in [0.0, 1.0] {
            let a = surface.point_at(0.0, viso);
            let b = surface.point_at(1.0, viso);
            let mut worst = (0.0_f64, 0.0_f64);
            for i in 0..=32 {
                let f = i as f64 / 32.0;
                let p = surface.point_at(f, viso);
                let ab = b - a;
                let t = (p - a).dot(&ab) / ab.dot(&ab);
                let q = popkern_core::geometry::Point3::new(
                    a.x + ab.x * t,
                    a.y + ab.y * t,
                    a.z + ab.z * t,
                );
                let d = p.distance_to(&q);
                if d > worst.0 {
                    worst = (d, f);
                }
            }
            eprintln!(
                "strip {fh:?} v={viso}: border straightness worst {:.4} at u={:.2}",
                worst.0, worst.1
            );
        }
        // Where do the strip's base-end corner vertices sit vs the
        // base plane z=0?
        let start = brep.loops[face.outer_loop].half_edge;
        let mut he = start;
        loop {
            let data = &brep.half_edges[he];
            let p = brep.vertices[data.origin].position;
            if p.z.abs() < 0.05 {
                eprintln!("  base-end vertex at ({:.4},{:.4},{:.4})", p.x, p.y, p.z);
            }
            he = data.next;
            if he == start {
                break;
            }
        }
        break;
    }
}
