//! Debug tests for cylinder-through-cube hole subtraction.
//!
//! The stress test revealed that cube-cylinder subtraction at tight tolerance
//! (1e-6) produces invalid topology, while at coarse tolerance (1e-2) it works.
//! These tests trace the issue.

use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn default_tol() -> Tolerance {
    Tolerance::default() // 1e-6
}

fn coarse_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-2)
}

fn medium_tol() -> Tolerance {
    Tolerance::default().with_linear(1e-4)
}

fn make_cube_with_hole(tol: Tolerance) -> Result<BRep, crate::error::BoolError> {
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.2, 2.0, tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    crate::subtraction(cube, cyl, tol, None)
}

#[test]
fn hole_at_coarse_tol_works() {
    let result = make_cube_with_hole(coarse_tol()).expect("coarse tol should work");
    let v = result.vertices.len();
    let e = result.edges.len();
    let f = result.faces.len();
    eprintln!("[coarse] V={v} E={e} F={f}");
    assert!(f >= 6, "should have at least 6 faces");
    let inner_loops: usize = result.faces.values().map(|f| f.inner_loops.len()).sum();
    eprintln!("[coarse] inner_loops={inner_loops}");
}

#[test]
fn hole_at_medium_tol() {
    let result = make_cube_with_hole(medium_tol());
    match result {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[medium] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[medium]   {err}");
            }
        }
        Err(e) => {
            eprintln!("[medium] FAILED: {e}");
        }
    }
}

#[test]
fn hole_at_default_tol() {
    let result = make_cube_with_hole(default_tol());
    match result {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[default] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[default]   {err}");
            }
            // Trace each loop to find the gap.
            use popkern_core::geometry::CurveEval;
            for (lh, ld) in &brep.loops {
                let start_he = ld.half_edge;
                let mut cursor = start_he;
                let mut first_pos = None;
                let mut end_pos = None;
                let mut count = 0;
                for _ in 0..100 {
                    let he = &brep.half_edges[cursor];
                    let pos = brep.vertices[he.origin].position;
                    if first_pos.is_none() {
                        first_pos = Some(pos);
                    }
                    let edge = &brep.edges[he.edge];
                    end_pos = if let Some(c) = &edge.curve {
                        let is_canonical = edge.half_edge == cursor;
                        Some(if is_canonical {
                            c.geometry.point_at(c.t_end)
                        } else {
                            c.geometry.point_at(c.t_start)
                        })
                    } else {
                        let next_he = &brep.half_edges[he.next];
                        Some(brep.vertices[next_he.origin].position)
                    };
                    count += 1;
                    cursor = he.next;
                    if cursor == start_he {
                        break;
                    }
                }
                if let (Some(fp), Some(ep)) = (first_pos, end_pos) {
                    let gap = fp.distance_to(&ep);
                    if gap > 1e-6 {
                        eprintln!("[gap] loop={lh:?} count={count} gap={gap:.6}");
                        eprintln!(
                            "  first=({:.4},{:.4},{:.4}) last_end=({:.4},{:.4},{:.4})",
                            fp.x, fp.y, fp.z, ep.x, ep.y, ep.z
                        );
                        // Also dump the last edge's curve info
                        cursor = start_he;
                        for _ in 0..count {
                            let he = &brep.half_edges[cursor];
                            let edge = &brep.edges[he.edge];
                            let has_curve = edge.curve.is_some();
                            let origin_pos = brep.vertices[he.origin].position;
                            eprintln!(
                                "  he origin=({:.4},{:.4},{:.4}) curve={has_curve}",
                                origin_pos.x, origin_pos.y, origin_pos.z
                            );
                            cursor = he.next;
                            if cursor == start_he {
                                break;
                            }
                        }
                    }
                }
            }

            // Trace vertex valence and face topology.
            for (vh, vertex) in &brep.vertices {
                let valence = brep
                    .half_edges
                    .iter()
                    .filter(|(_, he)| he.origin == vh)
                    .count();
                if valence < 2 {
                    let p = vertex.position;
                    eprintln!(
                        "[valence] vertex {vh:?} valence={valence} pos=({:.4},{:.4},{:.4})",
                        p.x, p.y, p.z
                    );
                    // Show which half-edges reference this vertex.
                    for (heh, he) in &brep.half_edges {
                        if he.origin == vh {
                            let edge = &brep.edges[he.edge];
                            let has_curve = edge.curve.is_some();
                            let twin_origin = brep.half_edges[he.twin].origin;
                            let is_self_loop = he.origin == twin_origin;
                            eprintln!(
                                "  he={heh:?} twin_origin={twin_origin:?} self_loop={is_self_loop} loop={:?} curve={has_curve}",
                                he.loop_
                            );
                        }
                    }
                }
            }
            // Trace face loops.
            for (fh, face) in &brep.faces {
                let outer_count = {
                    let mut c = 0;
                    let start = brep.loops[face.outer_loop].half_edge;
                    let mut cur = start;
                    for _ in 0..200 {
                        c += 1;
                        cur = brep.half_edges[cur].next;
                        if cur == start {
                            break;
                        }
                    }
                    c
                };
                let inner_counts: Vec<usize> = face
                    .inner_loops
                    .iter()
                    .map(|&lh| {
                        let mut c = 0;
                        let start = brep.loops[lh].half_edge;
                        let mut cur = start;
                        for _ in 0..200 {
                            c += 1;
                            cur = brep.half_edges[cur].next;
                            if cur == start {
                                break;
                            }
                        }
                        c
                    })
                    .collect();
                eprintln!(
                    "[face] {fh:?} outer={outer_count} inner_loops={} inner_counts={inner_counts:?}",
                    face.inner_loops.len()
                );
            }

            // Check intersection curve count
            let tol = default_tol();
            let cube = BRep::unit_cube(tol).unwrap();
            let cyl = BRep::cylinder(0.2, 2.0, tol)
                .unwrap()
                .transform(&Transform3::translation(0.5, 0.5, -0.5))
                .unwrap();
            let isect = crate::intersect::find_intersections(&cube, &cyl, &tol);
            match isect {
                Ok(result) => {
                    eprintln!(
                        "[default] intersection: {} curves, {} coincident overlaps",
                        result.records.len(),
                        result.coincidence.face_overlaps.len()
                    );
                    for (i, rec) in result.records.iter().enumerate() {
                        let (t0, t1) = rec.curve.curve_3d.parameter_range();
                        use popkern_core::geometry::CurveEval;
                        let p0 = rec.curve.curve_3d.point_at(t0);
                        let p1 = rec.curve.curve_3d.point_at(t1);
                        let closed = rec.curve.is_closed;
                        eprintln!(
                            "[default]   curve {i}: face_a={:?} face_b={:?} closed={closed}",
                            rec.face_a, rec.face_b
                        );
                        eprintln!(
                            "[default]     start=({:.3},{:.3},{:.3}) end=({:.3},{:.3},{:.3})",
                            p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                        );
                    }
                }
                Err(diag) => {
                    eprintln!("[default] intersection FAILED: {:?}", diag.error);
                }
            }
        }
        Err(e) => {
            eprintln!("[default] FAILED: {e}");
        }
    }
}

// ── Overlapping cubes subtraction (stress test 3) ────────────────────

#[test]
fn overlap_cubes_sub_debug() {
    let tol = default_tol();
    let a = BRep::unit_cube(tol).unwrap();
    let b = BRep::unit_cube(tol)
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.5, 0.0))
        .unwrap();

    // Check intersection details first.
    let isect = crate::intersect::find_intersections(&a, &b, &tol);
    match &isect {
        Ok(result) => {
            eprintln!(
                "[overlap] intersection: {} curves, {} hits, {} coincident",
                result.records.len(),
                result.edge_face_hits.len(),
                result.coincidence.face_overlaps.len()
            );
            for (i, rec) in result.records.iter().enumerate() {
                use popkern_core::geometry::CurveEval;
                let (t0, t1) = rec.curve.curve_3d.parameter_range();
                let p0 = rec.curve.curve_3d.point_at(t0);
                let p1 = rec.curve.curve_3d.point_at(t1);
                eprintln!(
                    "  curve {i}: a={:?} b={:?} tangent={} ({:.2},{:.2},{:.2})->({:.2},{:.2},{:.2})",
                    rec.face_a, rec.face_b, rec.is_tangent, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                );
            }
            for ov in &result.coincidence.face_overlaps {
                eprintln!(
                    "  coincident: a={:?} b={:?} orient={:?} full={} center=({:.2},{:.2},{:.2})",
                    ov.face_a,
                    ov.face_b,
                    ov.orientation,
                    ov.full_overlap,
                    ov.overlap_center.x,
                    ov.overlap_center.y,
                    ov.overlap_center.z,
                );
            }
        }
        Err(e) => eprintln!("[overlap] intersection FAILED: {e:?}"),
    }

    // Run subtraction.
    match crate::subtraction(a, b, tol, None) {
        Ok(brep) => {
            let v = brep.vertices.len();
            let e = brep.edges.len();
            let f = brep.faces.len();
            let errors = brep.validate_all();
            eprintln!("[overlap] V={v} E={e} F={f} errors={}", errors.len());
            for err in &errors {
                eprintln!("[overlap]   {err}");
            }
            for (fh, face) in &brep.faces {
                let outer_count = {
                    let mut c = 0;
                    let start = brep.loops[face.outer_loop].half_edge;
                    let mut cur = start;
                    for _ in 0..200 {
                        c += 1;
                        cur = brep.half_edges[cur].next;
                        if cur == start {
                            break;
                        }
                    }
                    c
                };
                if outer_count == 5 || outer_count <= 3 {
                    // Trace vertices of short loops.
                    let start = brep.loops[face.outer_loop].half_edge;
                    let mut cur = start;
                    for _ in 0..10 {
                        let he = &brep.half_edges[cur];
                        let p = brep.vertices[he.origin].position;
                        let twin_he = &brep.half_edges[he.twin];
                        let pe = brep.vertices[twin_he.origin].position;
                        let has_curve = brep.edges[he.edge].curve.is_some();
                        eprintln!(
                            "  vtx ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) curve={has_curve}",
                            p.x, p.y, p.z, pe.x, pe.y, pe.z
                        );
                        cur = he.next;
                        if cur == start {
                            break;
                        }
                    }
                }
                eprintln!("[overlap] face {fh:?} outer={outer_count}");
            }
        }
        Err(e) => eprintln!("[overlap] subtraction FAILED: {e}"),
    }
}
