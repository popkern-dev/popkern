use super::*;
use crate::types::{ChamferKind, ChamferSpec, FilletOptions, FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::collections::HashSet;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube() -> BRep {
    BRep::unit_cube(tol()).unwrap()
}

fn cube_and_edge() -> (BRep, popkern_core::topology::EdgeHandle) {
    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    (brep, edge)
}

/// Apply a constant-radius fillet and return the (pre_faces, post_brep) pair.
fn fillet_cube(radius: f64) -> (HashSet<FaceHandle>, BRep) {
    let (brep, edge) = cube_and_edge();
    let pre: HashSet<_> = brep.faces.keys().collect();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(radius),
    }];
    let result = crate::fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None)
        .expect("fillet should succeed");
    (pre, result.brep)
}

fn new_face_handles(pre: &HashSet<FaceHandle>, post: &BRep) -> Vec<FaceHandle> {
    post.faces.keys().filter(|fh| !pre.contains(fh)).collect()
}

// ---------------------------------------------------------------------------
// identify_blend: fillet on a plane-plane cube edge
// ---------------------------------------------------------------------------

#[test]
fn identify_plane_plane_cube_fillet_as_fillet_with_exact_radius() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);
    assert!(
        !new_faces.is_empty(),
        "fillet should create at least one face"
    );

    // Find the first new face the identifier calls a fillet.
    let info = new_faces
        .iter()
        .find_map(|fh| identify_blend(&post, *fh))
        .expect("at least one new face should identify as a blend");

    assert_eq!(info.kind, BlendKind::Fillet);
    let r = info.radius.expect("analytic fillets have known radius");
    assert!(
        (r - 0.1).abs() < 1e-6,
        "radius should match the requested value: got {r}, expected 0.1"
    );
    assert!(
        info.confidence >= 0.7,
        "fillet should be identified with reasonable confidence: got {}",
        info.confidence
    );
}

#[test]
fn identify_blend_returns_support_face_pair_for_fillet() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);

    let info = new_faces
        .iter()
        .find_map(|fh| identify_blend(&post, *fh))
        .expect("should find a blend face");

    let (s1, s2) = info
        .support_faces
        .expect("fillet should expose its support pair");
    // The two supports must be distinct from each other and from
    // the blend face itself. Note: the supports may be re-handled
    // sub-faces of the original cube faces (popkern's fillet
    // operation can split or re-issue the support face handles
    // after trimming).
    assert_ne!(s1, s2, "support faces should be distinct");
    assert_ne!(s1, info.face);
    assert_ne!(s2, info.face);
    let _ = pre;
}

// ---------------------------------------------------------------------------
// identify_all_blends: batch path
// ---------------------------------------------------------------------------

#[test]
fn identify_all_blends_finds_new_fillet_on_cube() {
    let (pre, post) = fillet_cube(0.1);
    let blends = identify_all_blends(&post);

    // At least one of the new faces should be recognised as a blend.
    let any_new = blends.iter().any(|info| !pre.contains(&info.face));
    assert!(
        any_new,
        "identify_all_blends should recognise at least one new fillet face"
    );

    // Every recognised face should carry a non-empty radius or at
    // least one surfaces-identified BlendKind.
    for info in &blends {
        assert!(
            matches!(
                info.kind,
                BlendKind::Fillet | BlendKind::Chamfer | BlendKind::VertexBlend
            ),
            "unexpected blend kind {:?}",
            info.kind
        );
    }
}

#[test]
fn identify_all_blends_rejects_original_cube_faces() {
    // On the un-filleted cube, every face is an original plane
    // bounded by four sharp 90° neighbours. None of them should
    // identify as a blend.
    let brep = cube();
    let blends = identify_all_blends(&brep);
    assert!(
        blends.is_empty(),
        "an un-blended cube should have no blend faces, got {} blends",
        blends.len()
    );
}

// ---------------------------------------------------------------------------
// identify_blend: negative cases
// ---------------------------------------------------------------------------

#[test]
fn identify_blend_returns_none_for_unknown_face_handle() {
    let brep = cube();
    let bogus = FaceHandle::default();
    assert!(identify_blend(&brep, bogus).is_none());
}

#[test]
fn identify_blend_returns_none_for_original_cube_face() {
    let brep = cube();
    let face = brep.faces.keys().next().unwrap();
    let info = identify_blend(&brep, face);
    assert!(
        info.is_none(),
        "original cube face should not identify as a blend, got {info:?}"
    );
}

// ---------------------------------------------------------------------------
// Convenience wrappers
// ---------------------------------------------------------------------------

#[test]
fn blend_radius_returns_same_value_as_identify_blend() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);

    let fh = new_faces
        .iter()
        .find(|fh| identify_blend(&post, **fh).is_some())
        .copied()
        .expect("some new face should identify as a blend");

    let via_identify = identify_blend(&post, fh).unwrap().radius;
    let via_shortcut = blend_radius(&post, fh);
    assert_eq!(via_identify, via_shortcut);
}

#[test]
fn blend_support_faces_returns_pair_for_fillet() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);

    let fh = new_faces
        .iter()
        .find(|fh| identify_blend(&post, **fh).is_some())
        .copied()
        .expect("some new face should identify as a blend");

    let supports = blend_support_faces(&post, fh);
    assert!(supports.is_some());
    let (a, b) = supports.unwrap();
    assert_ne!(a, b);
}

#[test]
fn blend_neighbours_returns_non_empty_for_fillet() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);

    let fh = new_faces
        .iter()
        .find(|fh| identify_blend(&post, **fh).is_some())
        .copied()
        .expect("some new face should identify as a blend");

    let neighbours = blend_neighbours(&post, fh);
    assert!(
        !neighbours.is_empty(),
        "a recognised fillet should have at least one neighbour"
    );
}

#[test]
fn blend_neighbours_returns_empty_for_unknown_face() {
    let brep = cube();
    let bogus = FaceHandle::default();
    assert!(blend_neighbours(&brep, bogus).is_empty());
}

// ---------------------------------------------------------------------------
// Chamfer identification
// ---------------------------------------------------------------------------

#[test]
fn identify_chamfer_on_cube_edge_as_chamfer() {
    let (brep, edge) = cube_and_edge();
    let pre: HashSet<_> = brep.faces.keys().collect();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    let result = crate::chamfer_edges(
        &brep,
        &specs,
        &crate::types::ChamferOptions::default(),
        tol(),
        None,
    )
    .expect("chamfer should succeed");

    let new_faces: Vec<_> = result
        .brep
        .faces
        .keys()
        .filter(|fh| !pre.contains(fh))
        .collect();

    // At least one of the new faces must identify as a chamfer.
    let chamfer_count = new_faces
        .iter()
        .filter_map(|fh| identify_blend(&result.brep, *fh))
        .filter(|info| info.kind == BlendKind::Chamfer)
        .count();

    assert!(
        chamfer_count >= 1,
        "chamfer operation should produce at least one Chamfer face; \
         found {chamfer_count} chamfer identifications out of {} new faces",
        new_faces.len()
    );
}

#[test]
fn identified_chamfer_reports_positive_radius_estimate() {
    let (brep, edge) = cube_and_edge();
    let pre: HashSet<_> = brep.faces.keys().collect();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.1),
    }];
    let result = crate::chamfer_edges(
        &brep,
        &specs,
        &crate::types::ChamferOptions::default(),
        tol(),
        None,
    )
    .unwrap();

    let chamfer_info = result
        .brep
        .faces
        .keys()
        .filter(|fh| !pre.contains(fh))
        .find_map(|fh| identify_blend(&result.brep, fh))
        .expect("chamfer should be identified");

    if chamfer_info.kind == BlendKind::Chamfer {
        let r = chamfer_info
            .radius
            .expect("chamfer should expose a radius estimate");
        assert!(r > 0.0, "chamfer radius estimate must be positive, got {r}");
    }
}

// ---------------------------------------------------------------------------
// Determinism and idempotence
// ---------------------------------------------------------------------------

#[test]
fn identify_blend_is_deterministic_across_calls() {
    let (_, post) = fillet_cube(0.1);
    let first: Vec<_> = identify_all_blends(&post)
        .into_iter()
        .map(|info| (info.kind, info.radius))
        .collect();
    let second: Vec<_> = identify_all_blends(&post)
        .into_iter()
        .map(|info| (info.kind, info.radius))
        .collect();
    assert_eq!(first, second);
}

#[test]
fn identify_blend_does_not_mutate_input() {
    let (_, post) = fillet_cube(0.1);
    let faces_before = post.faces.len();
    let edges_before = post.edges.len();
    let vertices_before = post.vertices.len();

    let _ = identify_all_blends(&post);

    assert_eq!(post.faces.len(), faces_before);
    assert_eq!(post.edges.len(), edges_before);
    assert_eq!(post.vertices.len(), vertices_before);
}

// ---------------------------------------------------------------------------
// BlendInfo fields are accessible and documented
// ---------------------------------------------------------------------------

#[test]
fn blend_info_source_is_curvature_analysis_for_this_round() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);
    let info = new_faces
        .iter()
        .find_map(|fh| identify_blend(&post, *fh))
        .expect("should find a blend");
    // uses curvature analysis only — metadata path is future work.
    assert_eq!(info.source, IdentificationSource::CurvatureAnalysis);
}

#[test]
fn blend_info_confidence_is_bounded_unit_interval() {
    let (pre, post) = fillet_cube(0.1);
    let new_faces = new_face_handles(&pre, &post);
    for fh in new_faces {
        if let Some(info) = identify_blend(&post, fh) {
            assert!(
                info.confidence >= 0.0 && info.confidence <= 1.0,
                "confidence must be in [0, 1], got {}",
                info.confidence
            );
            // threshold: no identification under 0.5.
            assert!(info.confidence >= 0.5);
        }
    }
}
