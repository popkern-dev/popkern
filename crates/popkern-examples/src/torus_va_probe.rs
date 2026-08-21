fn main() {
    let tol = popkern_core::tolerance::Tolerance::default();
    let t = popkern_core::topology::BRep::torus(1.0, 0.35, tol).expect("torus");
    let va = t.validate_all();
    println!("raw torus validate_all: {}", va.len());
    for f in va.iter().take(5) {
        println!("  {f:?}");
    }
    {
        let opts = popkern_mesh::TessellationOptions::default();
        if let Ok(mesh) = popkern_mesh::tessellate(&t, &opts) {
            let mut vol6 = 0.0_f64;
            for tri in &mesh.triangles {
                let a = mesh.vertices[tri[0]].position;
                let bb = mesh.vertices[tri[1]].position;
                let c = mesh.vertices[tri[2]].position;
                vol6 += a.x * (bb.y * c.z - bb.z * c.y) - a.y * (bb.x * c.z - bb.z * c.x)
                    + a.z * (bb.x * c.y - bb.y * c.x);
            }
            let analytic = 2.0 * std::f64::consts::PI.powi(2) * 1.0 * 0.35_f64.powi(2);
            println!(
                "pristine torus mesh volume {:.6} vs analytic {:.6} ({:+.2}%)",
                vol6.abs() / 6.0,
                analytic,
                (vol6.abs() / 6.0 - analytic) / analytic * 100.0
            );
        }
    }
    {
        // PK-309: SSI branch census for the rotated torus pair + witness probe.
        use popkern_core::geometry::Transform3;
        use popkern_core::geometry::curve::CurveEval;
        let tb = popkern_core::topology::BRep::torus(1.0, 0.35, tol)
            .expect("torus")
            .transform(
                &Transform3::rotation_x(0.5).compose(&Transform3::translation(0.7, 0.4, 0.15)),
            )
            .expect("transform");
        let sa = t
            .faces
            .iter()
            .next()
            .and_then(|(_, f)| f.surface.clone())
            .expect("surf");
        let sb = tb
            .faces
            .iter()
            .next()
            .and_then(|(_, f)| f.surface.clone())
            .expect("surf");
        let tau = std::f64::consts::TAU;
        let r = popkern_core::geometry::intersection::intersect_surfaces(
            &sa.geometry,
            (0.0, tau),
            (0.0, tau),
            &sb.geometry,
            (0.0, tau),
            (0.0, tau),
            &tol,
        );
        let witness = popkern_core::geometry::Point3::new(1.295304, 0.056062, 0.185951);
        match r {
            Ok(res) => {
                println!("SSI curves: {}", res.intersections.len());
                for (i, isect) in res.intersections.iter().enumerate() {
                    if let popkern_core::geometry::intersection::SurfaceSurfaceIntersection::Curve(
                        c,
                    ) = isect
                    {
                        let (t0, t1) = c.curve_3d.parameter_range();
                        let mut best = f64::MAX;
                        for k in 0..=256 {
                            let tt = t0 + (t1 - t0) * k as f64 / 256.0;
                            best = best.min(c.curve_3d.point_at(tt).distance_to(&witness));
                        }
                        let s0 = c.curve_3d.point_at(t0);
                        let s1 = c.curve_3d.point_at(t1);
                        println!(
                            "  curve {i}: witness-dist {best:.3e} span ({:.2},{:.2},{:.2})->({:.2},{:.2},{:.2})",
                            s0.x, s0.y, s0.z, s1.x, s1.y, s1.z
                        );
                    }
                }
            }
            Err(e) => println!("SSI err: {e:?}"),
        }
    }
    use popkern_core::geometry::curve::CurveEval;
    for (eh, ed) in t.edges.iter().take(4) {
        if let Some(c) = ed.curve.as_ref() {
            let pm = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            let pa = c.geometry.point_at(c.t_start);
            let pb = c.geometry.point_at(c.t_end);
            println!(
                "{eh:?} kind={} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) mid ({:.3},{:.3},{:.3})",
                match &c.geometry {
                    popkern_core::geometry::curve::CurveGeometry::Arc { .. } => "Arc",
                    popkern_core::geometry::curve::CurveGeometry::Line { .. } => "Line",
                    popkern_core::geometry::curve::CurveGeometry::Circle { .. } => "Circle",
                    _ => "other",
                },
                pa.x,
                pa.y,
                pa.z,
                pb.x,
                pb.y,
                pb.z,
                pm.x,
                pm.y,
                pm.z
            );
        } else {
            println!("{eh:?} NO CURVE");
        }
    }
}
