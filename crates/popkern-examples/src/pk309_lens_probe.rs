//! PK-309 lens probe: two tori at adjustable near-tangency, SSI
//! branch census vs boolean record census. Transient diagnostic bin.
//!
//! The rotated torus × torus union leaves a ~3e-3 osculating-branch
//! lens uncut. This probe measures, per tangency depth: how many SSI
//! branches exist near the osculation, how far apart they run, and
//! how many survive the boolean's clip + consolidate into records —
//! the discriminator data for the tangential-crossing certificate.

use popkern_core::geometry::Transform3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn main() {
    let tol = Tolerance::default();
    // Depth sweep: translate the second torus along x so the tubes
    // graze at varying depth. Tube centers R=1, r=0.35; two coplanar
    // tori side by side touch (externally, outer-to-outer) when the
    // center distance is 2(R+r) = 2.7. Sample depths straddling
    // tangency.
    for depth in [0.0005_f64, 0.002, 0.008, 0.03, 0.1] {
        let d = 2.0 * (1.0 + 0.35) - depth;
        let ta = BRep::torus(1.0, 0.35, tol).expect("torus a");
        let tb = BRep::torus(1.0, 0.35, tol)
            .expect("torus b")
            .transform(&Transform3::translation(d, 0.0, 0.0))
            .expect("translate");
        let sa = ta
            .faces
            .iter()
            .next()
            .and_then(|(_, f)| f.surface.clone())
            .expect("surf a");
        let sb = tb
            .faces
            .iter()
            .next()
            .and_then(|(_, f)| f.surface.clone())
            .expect("surf b");
        let tau = std::f64::consts::TAU;
        let res = popkern_core::geometry::intersection::intersect_surfaces(
            &sa.geometry,
            (0.0, tau),
            (0.0, tau),
            &sb.geometry,
            (0.0, tau),
            (0.0, tau),
            &tol,
        );
        match res {
            Ok(r) => {
                let mut spans = Vec::new();
                for isect in &r.intersections {
                    if let popkern_core::geometry::intersection::SurfaceSurfaceIntersection::Curve(
                        c,
                    ) = isect
                    {
                        let (t0, t1) = c.curve_3d.parameter_range();
                        let s = c.curve_3d.point_at(t0);
                        let e = c.curve_3d.point_at(t1);
                        let m = c.curve_3d.point_at(0.5 * (t0 + t1));
                        spans.push((s, m, e));
                    }
                }
                println!("depth {depth:.4}: {} SSI curves", spans.len());
                for (i, (s, m, e)) in spans.iter().enumerate() {
                    println!(
                        "  c{i}: ({:.3},{:.3},{:.3}) mid ({:.3},{:.3},{:.3}) -> ({:.3},{:.3},{:.3})",
                        s.x, s.y, s.z, m.x, m.y, m.z, e.x, e.y, e.z
                    );
                }
            }
            Err(e) => println!("depth {depth:.4}: SSI err {e:?}"),
        }
        // Boolean record census for the same pair.
        match popkern_bool::boolean(popkern_bool::BoolOperation::Union, ta, tb, tol, None) {
            Ok(b) => {
                let va = b.validate_all().len();
                let vol = {
                    let opts = popkern_mesh::TessellationOptions::default();
                    popkern_mesh::tessellate(&b, &opts).ok().map(|mesh| {
                        let mut vol6 = 0.0_f64;
                        for tri in &mesh.triangles {
                            let a = mesh.vertices[tri[0]].position;
                            let bb = mesh.vertices[tri[1]].position;
                            let c = mesh.vertices[tri[2]].position;
                            vol6 += a.x * (bb.y * c.z - bb.z * c.y)
                                - a.y * (bb.x * c.z - bb.z * c.x)
                                + a.z * (bb.x * c.y - bb.y * c.x);
                        }
                        vol6.abs() / 6.0
                    })
                };
                println!(
                    "  union: {} faces validate={:?} validate_all={va} mesh_vol {vol:?}",
                    b.faces.len(),
                    b.validate().err().map(|e| e.to_string())
                );
            }
            Err(e) => println!("  union ERR: {e}"),
        }
    }
}
