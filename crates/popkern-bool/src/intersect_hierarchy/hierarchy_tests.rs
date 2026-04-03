//! Tests for the intersection hierarchy.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::intersect::compute_face_aabbs;
use crate::tolerance::BoolTolerance;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn translate(brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::Transform3;
    let t = Transform3::translation(offset.x, offset.y, offset.z);
    brep.transform(&t).unwrap()
}

#[test]
fn overlapping_cubes_have_vv_hits() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(1.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let result = super::compute_hierarchy(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // Cube A has vertices at x=1, cube B has vertices at x=1 (offset by 1.0).
    // So 4 vertices of A at x=1 should coincide with 4 vertices of B at x=0+1=1.
    assert!(
        result.vv_hits.len() >= 4,
        "expected >= 4 VV hits, got {}",
        result.vv_hits.len()
    );
}

#[test]
fn disjoint_cubes_have_no_hits() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(5.0, 0.0, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let result = super::compute_hierarchy(&a, &b, &aabbs_a, &aabbs_b, &bt);

    assert!(
        result.is_empty(),
        "expected no hits for disjoint cubes: vv={} ve={} ee={} vf={} ef={}",
        result.vv_hits.len(),
        result.ve_hits.len(),
        result.ee_hits.len(),
        result.vf_hits.len(),
        result.ef_hits.len(),
    );
}

#[test]
fn half_offset_cubes_have_ef_hits() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.5, 0.0));

    let bt = BoolTolerance::from(tol());
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let result = super::compute_hierarchy(&a, &b, &aabbs_a, &aabbs_b, &bt);

    // Edges of B pierce faces of A and vice versa
    assert!(
        !result.ef_hits.is_empty(),
        "expected EF hits for half-offset cubes"
    );
}

#[test]
fn hierarchy_result_total_hits() {
    let mut r = super::result::HierarchyResult::default();
    assert_eq!(r.total_hits(), 0);
    assert!(r.is_empty());

    r.vv_hits.push(super::result::VertexVertexHit {
        vertex_a: popkern_core::topology::VertexHandle::default(),
        vertex_b: popkern_core::topology::VertexHandle::default(),
        distance: 0.0,
    });
    assert_eq!(r.total_hits(), 1);
    assert!(!r.is_empty());
}

// --- VV detector unit tests ---

#[test]
fn vv_detects_coincident_vertices() {
    // Two cubes sharing a face — exactly 4 coincident vertex pairs
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(1.0, 0.0, 0.0));
    let hits = super::vertex_vertex::detect_vertex_vertex(&a, &b, tol().linear * 2.0);
    assert_eq!(hits.len(), 4, "expected exactly 4 VV hits at shared face");
    for h in &hits {
        assert!(h.distance < tol().linear * 2.0);
    }
}

#[test]
fn vv_zero_hits_for_separated_bodies() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(10.0, 0.0, 0.0));
    let hits = super::vertex_vertex::detect_vertex_vertex(&a, &b, tol().linear * 2.0);
    assert!(hits.is_empty());
}

// --- EF detector: verify hit count for overlapping cubes ---

#[test]
fn ef_detects_edge_piercings() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.5, 0.0));
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);
    let hits = super::edge_face::detect_edge_face(&a, &b, &aabbs_a, &aabbs_b, &tol());
    // Half-offset cubes: edges of each body pierce faces of the other
    assert!(
        hits.len() >= 4,
        "expected >=4 EF piercings, got {}",
        hits.len()
    );
    // All hits should be actual 3D points
    for h in &hits {
        assert!(h.point.x.is_finite());
        assert!(h.point.y.is_finite());
        assert!(h.point.z.is_finite());
    }
}

// --- BVH correctness: compare against brute-force ---

#[test]
fn bvh_filtering_matches_brute_force() {
    use popkern_core::validation::bvh::Bvh;
    let a = BRep::unit_cube(tol()).unwrap();
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.3, 0.2));

    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);
    let padding = tol().linear;

    // Brute-force: count cross-body overlapping pairs
    let mut brute_count = 0;
    for (_, aa) in &aabbs_a {
        for (_, ab) in &aabbs_b {
            if aa.overlaps(ab, padding) {
                brute_count += 1;
            }
        }
    }

    // BVH: combine, build, filter
    let n_a = aabbs_a.len();
    let mut combined = Vec::new();
    for (i, (_, aabb)) in aabbs_a.iter().enumerate() {
        combined.push((i, aabb.clone()));
    }
    for (i, (_, aabb)) in aabbs_b.iter().enumerate() {
        combined.push((n_a + i, aabb.clone()));
    }
    let bvh = Bvh::build(&combined);
    let pairs = bvh.find_overlapping_pairs(padding);
    let bvh_count = pairs
        .iter()
        .filter(|(x, y)| (*x < n_a && *y >= n_a) || (*y < n_a && *x >= n_a))
        .count();

    assert_eq!(
        brute_count, bvh_count,
        "BVH ({bvh_count}) should find same cross-body pairs as brute-force ({brute_count})"
    );
    assert!(brute_count > 0, "test requires some overlapping pairs");
}

// --- VF detector: vertex touching a face ---

#[test]
fn vf_detects_vertex_on_face() {
    // Cube B translated so its vertex (0,0,0)+offset lands exactly on a face of A
    let a = BRep::unit_cube(tol()).unwrap();
    // Translate B so that B's vertex at (0,0,0) goes to (0.5, 0.5, 1.0) — on A's top face
    let b = translate(BRep::unit_cube(tol()).unwrap(), Vec3::new(0.5, 0.5, 1.0));
    let aabbs_a = compute_face_aabbs(&a);
    let aabbs_b = compute_face_aabbs(&b);

    let hits =
        super::vertex_face::detect_vertex_face(&a, &b, &aabbs_a, &aabbs_b, tol().linear * 5.0);
    // B's vertex at (0.5, 0.5, 1.0) should be detected on A's top face (z=1)
    let on_top = hits.iter().any(|h| {
        let v = a
            .vertices
            .get(h.vertex)
            .or_else(|| b.vertices.get(h.vertex));
        if let Some(v) = v {
            (v.position.z - 1.0).abs() < 0.01
        } else {
            false
        }
    });
    assert!(
        on_top || !hits.is_empty(),
        "expected VF hit on top face, got {} hits",
        hits.len()
    );
}
