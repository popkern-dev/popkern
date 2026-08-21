//! PK-304/PK-305 verification renders: torus × cylinder and
//! torus × torus unions. Transient inspection aid.

#[allow(dead_code)]
mod helpers;

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn main() {
    let tol = Tolerance::default();
    let torus = BRep::torus(1.0, 0.35, tol).expect("torus builds");
    let cyl = BRep::cylinder(0.5, 3.0, tol)
        .expect("cylinder builds")
        .transform(&Transform3::translation(0.8, 0.0, -1.5))
        .expect("translate");
    match boolean(BoolOperation::Union, torus.clone(), cyl, tol, None) {
        Ok(b) => {
            let boundary_hes = b
                .half_edges
                .iter()
                .filter(|(_, he)| !b.loops.contains_key(he.loop_))
                .count();
            println!(
                "torus x cylinder: {} faces validate={:?} boundary_hes={} validate_all={}",
                b.faces.len(),
                b.validate().err().map(|e| e.to_string()),
                boundary_hes,
                b.validate_all().len()
            );
            for f in b.validate_all() {
                println!("  [va] {f:?}");
                if let popkern_core::error::ValidationError::CurveSurfaceDeviation { edge, .. } = &f
                    && let Some(ed) = b.edges.get(*edge)
                {
                    {
                        use popkern_core::geometry::curve::CurveEval;
                        if let Some(c) = ed.curve.as_ref() {
                            let pm = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                            let pa = c.geometry.point_at(c.t_start);
                            let pb = c.geometry.point_at(c.t_end);
                            let he = ed.half_edge;
                            let f1 = b.loops[b.half_edges[he].loop_].face;
                            let f2 = b.loops[b.half_edges[b.half_edges[he].twin].loop_].face;
                            let sk = |fh| {
                                b.faces
                                    .get(fh)
                                    .and_then(|fd| fd.surface.as_ref())
                                    .map(|sf| format!("{:?}", std::mem::discriminant(&sf.geometry)))
                            };
                            println!(
                                "      ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) mid ({:.3},{:.3},{:.3}) faces {f1:?}{:?} {f2:?}{:?}",
                                pa.x,
                                pa.y,
                                pa.z,
                                pb.x,
                                pb.y,
                                pb.z,
                                pm.x,
                                pm.y,
                                pm.z,
                                sk(f1),
                                sk(f2)
                            );
                        }
                    }
                }
            }
            // Dump open loops: walk each loop, report face + endpoint gap.
            for (lh, lp) in &b.loops {
                let start = lp.half_edge;
                if !b.half_edges.contains_key(start) {
                    continue;
                }
                let mut cur = start;
                let mut n = 0;
                loop {
                    let he = &b.half_edges[cur];
                    let next = he.next;
                    if !b.half_edges.contains_key(next) {
                        break;
                    }
                    n += 1;
                    if n > 500 {
                        break;
                    }
                    cur = next;
                    if cur == start {
                        break;
                    }
                }
                let he0 = &b.half_edges[start];
                let first = b.vertices[he0.origin].position;
                let last_he = &b.half_edges[cur];
                let _ = last_he;
                // recompute walk end position
                let mut cur2 = start;
                let mut endpos = first;
                for _ in 0..n {
                    let he = &b.half_edges[cur2];
                    let twin = &b.half_edges[he.twin];
                    endpos = b.vertices[twin.origin].position;
                    cur2 = he.next;
                }
                let gap = endpos.distance_to(&first);
                if gap > 1e-6 {
                    println!(
                        "OPEN loop {lh:?} on face {:?}: {} hes, start ({:.3},{:.3},{:.3}) end ({:.3},{:.3},{:.3}) gap {gap:.3}",
                        lp.face, n, first.x, first.y, first.z, endpos.x, endpos.y, endpos.z
                    );
                    let mut cur3 = start;
                    for _ in 0..n.min(12) {
                        let he = &b.half_edges[cur3];
                        let twin = &b.half_edges[he.twin];
                        let p0 = b.vertices[he.origin].position;
                        let p1 = b.vertices[twin.origin].position;
                        println!(
                            "    he {:?} edge {:?} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                            cur3, he.edge, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                        );
                        cur3 = he.next;
                    }
                }
            }
            for (eh, ed) in &b.edges {
                use popkern_core::geometry::curve::CurveEval;
                let he = &b.half_edges[ed.half_edge];
                let f1 = b.loops[he.loop_].face;
                let f2 = b.loops[b.half_edges[he.twin].loop_].face;
                let (kind, pa, pb) = match ed.curve.as_ref() {
                    Some(c) => (
                        match &c.geometry {
                            popkern_core::geometry::curve::CurveGeometry::Arc { .. } => "Arc",
                            popkern_core::geometry::curve::CurveGeometry::Line { .. } => "Line",
                            popkern_core::geometry::curve::CurveGeometry::Circle { .. } => "Circle",
                            popkern_core::geometry::curve::CurveGeometry::Nurbs(_) => "Nurbs",
                            _ => "other",
                        },
                        c.geometry.point_at(c.t_start),
                        c.geometry.point_at(c.t_end),
                    ),
                    None => (
                        "NONE",
                        b.vertices[he.origin].position,
                        b.vertices[b.half_edges[he.twin].origin].position,
                    ),
                };
                println!(
                    "  [edge] {eh:?} {kind} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) faces {f1:?}/{f2:?}",
                    pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
                );
            }
            for (heh, he) in &b.half_edges {
                let self_twin = he.twin == heh;
                if !self_twin {
                    continue;
                }
                let nxt = b.half_edges.get(he.next);
                let d_next = nxt
                    .and_then(|n| b.vertices.get(n.origin))
                    .zip(b.vertices.get(he.origin))
                    .map(|(a, o)| a.position.distance_to(&o.position));
                println!(
                    "  [selftwin] {heh:?} edge {:?} next-orig-dist {d_next:?} in_loop {:?} outer={:?}",
                    he.edge,
                    he.loop_,
                    b.loops
                        .get(he.loop_)
                        .map(|l| b.faces[l.face].outer_loop == he.loop_)
                );
            }
            for (fh, fd) in &b.faces {
                let sk = fd
                    .surface
                    .as_ref()
                    .map(|s| match &s.geometry {
                        popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
                        popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => "Cyl",
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                        _ => "other",
                    })
                    .unwrap_or("none");
                println!("  [face] {fh:?} {sk} loops 1+{}", fd.inner_loops.len());
            }
            helpers::export(&b, "PK-304 torus x cylinder union", "pk304_torus_cyl.step");
            helpers::render_multi_angle(&b, "pk304_torus_cyl");
        }
        Err(e) => println!("torus x cylinder ERR: {e}"),
    }
    let torus_b = BRep::torus(1.0, 0.35, tol)
        .expect("torus builds")
        .transform(&Transform3::rotation_x(0.5).compose(&Transform3::translation(0.7, 0.4, 0.15)))
        .expect("transform");
    match boolean(BoolOperation::Union, torus, torus_b, tol, None) {
        Ok(b) => {
            let boundary_hes = b
                .half_edges
                .iter()
                .filter(|(_, he)| !b.loops.contains_key(he.loop_))
                .count();
            println!(
                "torus x torus: {} faces validate={:?} boundary_hes={} validate_all={}",
                b.faces.len(),
                b.validate().err().map(|e| e.to_string()),
                boundary_hes,
                b.validate_all().len()
            );
            helpers::export(&b, "PK-305 torus x torus union", "pk305_torus_torus.step");
            helpers::render_multi_angle(&b, "pk305_torus_torus");
        }
        Err(e) => println!("torus x torus ERR: {e}"),
    }
}
