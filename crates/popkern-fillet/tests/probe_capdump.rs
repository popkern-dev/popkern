use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_fillet::fillet_edges;
use popkern_fillet::types::{FilletOptions, FilletRadius, FilletSpec};

fn incident(brep: &BRep, vertex: popkern_core::topology::VertexHandle) -> Vec<EdgeHandle> {
    let mut seen: Vec<EdgeHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        let edge = he.edge;
        if !seen.contains(&edge) {
            seen.push(edge);
        }
    }
    seen
}

#[test]
fn dump() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    }
    unsafe {
        std::env::set_var("POPKERN_A2_TRACE", "1");
    }
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let vertex = brep.vertices.keys().next().unwrap();
    let edges = incident(&brep, vertex);
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        use_a2_apex_construction: true,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&brep, &specs, &options, tol, None).expect("fillet_edges ok");
    for (_, w) in &result.edge_warnings {
        if let popkern_fillet::error::BlendWarning::A2FellBackToLegacy { reason, .. } = w {
            eprintln!("[probe] A2 fallback reason: {reason}");
        }
    }
    eprintln!(
        "[probe] post-pipeline: faces={} edges={}",
        result.brep.faces.len(),
        result.brep.edges.len()
    );
}

#[test]
fn dump_pyramid3_uniform() {
    use popkern_core::topology::BRep;
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    }
    unsafe {
        std::env::set_var("POPKERN_A2_TRACE", "1");
    }
    let tol = Tolerance::default();
    let brep = BRep::n_gonal_pyramid(3, 1.0, 0.8, tol).unwrap();
    // Find apex vertex (the one not on z=0 plane).
    let apex_vh = brep
        .vertices
        .iter()
        .filter(|(_, v)| v.position.z.abs() > 1e-6)
        .map(|(vh, _)| vh)
        .next()
        .unwrap();
    let edges_at_apex: Vec<popkern_core::topology::EdgeHandle> = brep
        .edges
        .iter()
        .filter(|(_, e)| {
            let primary = e.half_edge;
            let he = &brep.half_edges[primary];
            let twin_he = &brep.half_edges[he.twin];
            he.origin == apex_vh || twin_he.origin == apex_vh
        })
        .map(|(eh, _)| eh)
        .collect();
    let radii = [0.10_f64, 0.10, 0.10];
    let specs: Vec<FilletSpec> = edges_at_apex
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        use_a2_apex_construction: true,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&brep, &specs, &options, tol, None).expect("fillet_edges ok");
    for (_, w) in &result.edge_warnings {
        if let popkern_fillet::error::BlendWarning::A2FellBackToLegacy { reason, .. } = w {
            eprintln!("[probe-pyr3] A2 fallback reason: {reason}");
        }
    }
}

#[test]
fn dump_legacy() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    }
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let vertex = brep.vertices.keys().next().unwrap();
    let edges = incident(&brep, vertex);
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        use_a2_apex_construction: false,
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&brep, &specs, &options, tol, None).expect("fillet_edges ok");
    eprintln!(
        "[probe-legacy] post-pipeline: faces={} edges={}",
        result.brep.faces.len(),
        result.brep.edges.len()
    );
    for cp in &result.corner_patches {
        let face = cp.face;
        let Some(face_data) = result.brep.faces.get(face) else {
            continue;
        };
        let Some(loop_data) = result.brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = result.brep.half_edges.len() + 1;
        eprintln!("[probe-legacy] cap_face={face:?} method={:?}", cp.method);
        while let Some(he) = result.brep.half_edges.get(cursor) {
            let twin_face = result
                .brep
                .half_edges
                .get(he.twin)
                .and_then(|t| result.brep.loops.get(t.loop_))
                .map(|l| l.face);
            let p = result.brep.vertices.get(he.origin).map(|v| v.position);
            eprintln!(
                "[probe-legacy]   [{iters}] edge={:?} v={:?}@{:?} twin_face={twin_face:?}",
                he.edge, he.origin, p
            );
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
    }
}
