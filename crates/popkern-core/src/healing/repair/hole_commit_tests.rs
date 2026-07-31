use super::super::coons_patch::CoonsPatchResult;
use super::super::hole_commit_geometry::HoleCapGeometry;
use super::super::hole_fill::{DetectedHole, detect_holes};
use super::super::n_sided_fill_g1::fill_n_sided_hole_g1;
use super::*;
use crate::geometry::NurbsSurface;
use crate::tolerance::Tolerance;
use crate::topology::VertexHandle;

fn open_box_with_hole() -> (BRep, DetectedHole) {
    let brep = BRep::open_box(Tolerance::default()).expect("open box");
    let holes = detect_holes(&brep);
    assert_eq!(holes.len(), 1, "open box has exactly one hole");
    let hole = holes[0].clone();
    (brep, hole)
}

fn linear_tolerance() -> f64 {
    Tolerance::default().linear()
}

#[test]
fn pns_fill_declines_with_its_measured_seam_gap() {
    // The PnS bi-4 complex from `fill_n_sided_boundary` is NOT
    // fan-committable: only the CP-shared spoke columns coincide
    // between sectors, and the fan seam curves sit ~1e-1 off the
    // adjacent patch on the open-box fill (containment measured; both
    // wedge assignments falsified — see the PK-199 resolution notes).
    // The coherence gate is the authority: the commit must decline
    // with the real gap and leave the body untouched.
    let (mut brep, hole) = open_box_with_hole();
    let tolerance = linear_tolerance();
    let fill = fill_n_sided_hole_g1(&brep, &hole, tolerance).expect("G1 fill");
    let geometry = HoleCapGeometry::from_pns_fill(&fill, tolerance).expect("normalize");
    let before = serde_json::to_string(&brep).expect("serializes");
    let error = commit_hole_cap(&mut brep, &hole, &geometry, tolerance)
        .expect_err("the overlapping-sector complex must not commit");
    match error {
        HoleCommitError::CoherenceRejected { l2_worst, .. } => {
            assert!(
                l2_worst > 1.0e-3,
                "decline must record the real seam gap, got {l2_worst:.3e}"
            );
        }
        other => panic!("expected a coherence decline, got {other}"),
    }
    assert_eq!(serde_json::to_string(&brep).expect("serializes"), before);
}

#[test]
fn ruled_fan_commits_and_body_validates_strictly() {
    let (mut brep, hole) = open_box_with_hole();
    let tolerance = linear_tolerance();
    let geometry =
        super::super::hole_ruled_fan::ruled_fan_cap(&brep, &hole, tolerance).expect("ruled fan");
    let commit = commit_hole_cap(&mut brep, &hole, &geometry, tolerance).expect("commit");

    assert!(brep.boundary_edges.is_empty());
    assert_eq!(commit.fan.sector_faces.len(), 4);
    assert!(
        commit.worst_l2_gap <= 1.0e-9,
        "ruled sectors contain their boundaries by construction, got L2 {:.3e}",
        commit.worst_l2_gap
    );
    brep.validate_euler_poincare().expect("strict Euler holds");
    brep.validate().expect("full fast validation holds");
    for &face in &commit.fan.sector_faces {
        assert!(brep.faces[face].surface.is_some());
    }
    for &spoke in &commit.fan.spoke_edges {
        assert!(brep.edges[spoke].curve.is_some());
    }
}

#[test]
fn sector_count_mismatch_declines_untouched() {
    let (mut brep, hole) = open_box_with_hole();
    let tolerance = linear_tolerance();
    let fill = fill_n_sided_hole_g1(&brep, &hole, tolerance).expect("G1 fill");
    let mut geometry = HoleCapGeometry::from_pns_fill(&fill, tolerance).expect("normalize");
    geometry.sectors.pop();
    let before = serde_json::to_string(&brep).expect("serializes");
    let error = commit_hole_cap(&mut brep, &hole, &geometry, tolerance)
        .expect_err("3 sectors for a 4-sided hole decline");
    assert!(
        matches!(
            error,
            HoleCommitError::SectorCountMismatch {
                expected: 4,
                got: 3
            }
        ),
        "{error}"
    );
    assert_eq!(serde_json::to_string(&brep).expect("serializes"), before);
}

/// Ordered ring of hole vertices, following the detected edge chain.
fn hole_vertex_ring(brep: &BRep, hole: &DetectedHole) -> Vec<Point3> {
    let endpoints: Vec<(VertexHandle, VertexHandle)> = hole
        .boundary_edges
        .iter()
        .map(|&edge| {
            let he = brep.edges[edge].half_edge;
            let data = &brep.half_edges[he];
            (data.origin, brep.half_edges[data.next].origin)
        })
        .collect();
    // Chain: pick each next vertex as the endpoint shared with the
    // following edge.
    let mut ring = Vec::with_capacity(endpoints.len());
    let first_shared = {
        let (a, b) = endpoints[0];
        let (c, d) = endpoints[1];
        if b == c || b == d { b } else { a }
    };
    let mut current = first_shared;
    ring.push(current);
    for &(a, b) in endpoints.iter().skip(1) {
        current = if a == current { b } else { a };
        ring.push(current);
    }
    ring.iter()
        .map(|&vertex| brep.vertices[vertex].position)
        .collect()
}

#[test]
fn coons_shaped_patch_set_commits_through_the_same_operator() {
    // PK-199 acceptance 314: the commit path is filler-agnostic. The
    // production Coons builder declines fan openings (its collapsed
    // top edge kills patch construction — the recorded PK-195 gap),
    // so the fan-shaped patch set is synthesized here in exactly the
    // container the Coons filler returns: per-sector surfaces plus
    // corner points, spokes rebuilt by the normalizer.
    let (mut brep, hole) = open_box_with_hole();
    let tolerance = linear_tolerance();
    let ring = hole_vertex_ring(&brep, &hole);
    let count = ring.len();
    let inverse = 1.0 / count as f64;
    let mut pole = Point3::new(0.0, 0.0, 0.0);
    for point in &ring {
        pole = Point3::new(
            pole.x + point.x * inverse,
            pole.y + point.y * inverse,
            pole.z + point.z * inverse,
        );
    }

    let patches: Vec<CoonsPatchResult> = (0..count)
        .map(|k| {
            let rim_start = ring[k];
            let rim_end = ring[(k + 1) % count];
            // Bilinear sector with the top edge collapsed at the pole
            // — the Coons fan's patch shape over a planar opening.
            let surface = NurbsSurface {
                degree_u: 1,
                degree_v: 1,
                knots_u: vec![0.0, 0.0, 1.0, 1.0],
                knots_v: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![rim_start, pole, rim_end, pole],
                num_u: 2,
                num_v: 2,
                weights: vec![1.0; 4],
            };
            CoonsPatchResult {
                surface,
                max_boundary_deviation: 0.0,
                corners: [rim_start, rim_end, pole, pole],
            }
        })
        .collect();

    let geometry = HoleCapGeometry::from_coons_patches(&patches, hole.num_sides, tolerance)
        .expect("normalize");
    let commit = commit_hole_cap(&mut brep, &hole, &geometry, tolerance).expect("commit");

    assert!(brep.boundary_edges.is_empty());
    assert_eq!(commit.fan.sector_faces.len(), count);
    brep.validate_euler_poincare().expect("strict Euler holds");
    brep.validate().expect("full fast validation holds");
}

#[test]
fn incomplete_coons_patch_set_declines() {
    let (brep, hole) = open_box_with_hole();
    let tolerance = linear_tolerance();
    let _ = brep;
    let error = HoleCapGeometry::from_coons_patches(&[], hole.num_sides, tolerance)
        .expect_err("an empty patch set is not a fan");
    assert!(
        matches!(
            error,
            HoleCommitError::SectorCountMismatch {
                expected: 4,
                got: 0
            }
        ),
        "{error}"
    );
}
