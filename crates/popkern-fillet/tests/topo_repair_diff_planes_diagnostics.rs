use super::*;

#[ignore]
#[test]
fn t_136_b3_anchor_oncyl_status() {
    // Pyramid apex-only fillets can produce chains of off-cylinder
    // vertices in cylinder loops. This diagnostic counts whether
    // chord-pair neighbours are valid cylinder anchors.
    eprintln!("\n=== #136 anchor on-cyl status per candidate ===");
    for n in 4..=8 {
        let brep = fillet_apex_only(n, 0.8, 0.05);
        let mut clean = 0;
        let mut polluted = 0;
        for (_fh, face) in brep.faces.iter() {
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            if !matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) {
                continue;
            }
            let outer_loop = face.outer_loop;
            let Some(loop_data) = brep.loops.get(outer_loop) else {
                continue;
            };
            let start_he = loop_data.half_edge;
            let mut hes = Vec::new();
            let mut cursor = start_he;
            let mut safety = 64;
            loop {
                if safety == 0 {
                    break;
                }
                safety -= 1;
                hes.push(cursor);
                cursor = brep.half_edges[cursor].next;
                if cursor == start_he {
                    break;
                }
            }
            let nh = hes.len();
            for i in 0..nh {
                let he_out = hes[i];
                let he_in = hes[(i + nh - 1) % nh];
                let v = brep.half_edges[he_out].origin;
                let p = brep.vertices[v].position;
                if vertex_on_cylinder(&p, &sof.geometry) {
                    continue;
                }
                let anchor_in_v = brep.half_edges[he_in].origin;
                let twin_out = brep.half_edges[he_out].twin;
                let anchor_out_v = brep.half_edges[twin_out].origin;
                let p_in = brep.vertices[anchor_in_v].position;
                let p_out = brep.vertices[anchor_out_v].position;
                let in_on = vertex_on_cylinder(&p_in, &sof.geometry);
                let out_on = vertex_on_cylinder(&p_out, &sof.geometry);
                if in_on && out_on {
                    clean += 1;
                } else {
                    polluted += 1;
                }
            }
        }
        eprintln!(
            "  N={n}: clean (both anchors on-cyl) = {clean}, polluted (>=1 anchor off-cyl) = {polluted}"
        );
    }
}

#[ignore]
#[test]
fn t_136_b2_attempt_cap_sweep() {
    eprintln!("\n=== #136 attempt-cap sweep: tess-nm vs MAX_ATTEMPTS ===");
    eprintln!(
        "{:>3}  {:>14}  {:>14}  {:>14}  {:>14}",
        "N", "attempts=1", "attempts=2", "attempts=4", "attempts=8"
    );
    eprintln!("{:-<76}", "");
    for n in 4..=8 {
        let mut row = Vec::new();
        for k in [1usize, 2, 4, 8] {
            let brep = run_with_attempts(n, 0.8, 0.05, k);
            row.push(worst_tess_nm_divergence(&brep));
        }
        eprintln!(
            "{:>3}  {:>14.4e}  {:>14.4e}  {:>14.4e}  {:>14.4e}",
            n, row[0], row[1], row[2], row[3]
        );
    }
}

#[ignore]
#[test]
fn t_136_c_diff_planes_neighbour_drift() {
    eprintln!("\n=== #136 diff-planes neighbour drift: where does v_off snap break things? ===");
    let n = 7;
    let height = 0.8;
    let radius = 0.05;

    let pre = run_with_flags(n, height, radius, true, false);
    let pre_survey = survey_off_cyl(&pre);
    eprintln!(
        "After same-plane only (N={n}): cyl_faces={} off_cyl_v={} same_plane={} diff_planes={}",
        pre_survey.cylinder_faces,
        pre_survey.off_cyl_vertices,
        pre_survey.same_plane_chord_pairs,
        pre_survey.diff_plane_chord_pairs
    );

    let mut off_cyl_vs = std::collections::BTreeSet::<VertexHandle>::new();
    for (_fh, face) in pre.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let outer_loop = face.outer_loop;
        let Some(loop_data) = pre.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        let mut cursor = start_he;
        let mut safety = 64;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            let v = pre.half_edges[cursor].origin;
            let p = pre.vertices[v].position;
            if !vertex_on_cylinder(&p, &sof.geometry) {
                off_cyl_vs.insert(v);
            }
            cursor = pre.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }
    }
    eprintln!("Off-cyl vertices to inspect: {}", off_cyl_vs.len());

    let post = run_with_flags(n, height, radius, true, true);
    eprintln!(
        "Post diff-planes:  worst tess-nm = {:.3e}",
        worst_tess_nm_divergence(&post)
    );

    for v in off_cyl_vs {
        let p_pre = pre.vertices[v].position;
        let (v_post, p_post) = post
            .vertices
            .iter()
            .min_by(|a, b| {
                let da = a.1.position.distance_to(&p_pre);
                let db = b.1.position.distance_to(&p_pre);
                da.partial_cmp(&db).unwrap()
            })
            .map(|(h, vd)| (h, vd.position))
            .expect("non-empty post");

        let snap_dist = p_pre.distance_to(&p_post);
        eprintln!(
            "\n  v_pre={v:?} pos=({:.4},{:.4},{:.4}) -> v_post={v_post:?} pos=({:.4},{:.4},{:.4}) snap={:.4e} m",
            p_pre.x, p_pre.y, p_pre.z, p_post.x, p_post.y, p_post.z, snap_dist
        );

        let drift_pre = surfaces_at_vertex(&pre, v);
        let drift_post = surfaces_at_vertex(&post, v_post);

        let max_pre = drift_pre.iter().map(|t| t.2).fold(0.0, f64::max);
        let max_post = drift_post.iter().map(|t| t.2).fold(0.0, f64::max);
        eprintln!(
            "    incident surfaces (pre): n={}, max_drift={:.4e}",
            drift_pre.len(),
            max_pre
        );
        for (fh, kind, dist) in &drift_pre {
            eprintln!("      pre  {fh:?} {kind:>8}  drift={:.4e}", dist);
        }
        eprintln!(
            "    incident surfaces (post): n={}, max_drift={:.4e}",
            drift_post.len(),
            max_post
        );
        for (fh, kind, dist) in &drift_post {
            eprintln!("      post {fh:?} {kind:>8}  drift={:.4e}", dist);
        }
    }
}
