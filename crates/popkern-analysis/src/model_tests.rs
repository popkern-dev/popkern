use crate::model::*;

/// A closed, outward-oriented unit tetrahedron with one face per
/// triangle.
fn tetrahedron() -> AnalysisModel {
    let vertices = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let triangles = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
    let faces = (0..4)
        .map(|i| FaceInfo {
            id: FaceId(i),
            solid: SolidId(0),
            name: None,
            surface: Surface::Freeform,
            role: None,
            role_source: RoleSource::Inferred,
            area: 0.5,
            exact_area: None,
            folded_area: None,
        })
        .collect();
    AnalysisModel {
        units: Units::Millimetre,
        chord_tolerance: 0.01,
        solids: vec![SolidMesh {
            id: SolidId(0),
            name: None,
            vertices,
            triangles,
            triangle_face: (0..4).map(FaceId).collect(),
            exact_volume: Some(1.0 / 6.0),
        }],
        faces,
        concave_edges: Vec::new(),
        provenance: Provenance::Step {
            file_sha256: "0".repeat(64),
        },
    }
}

#[test]
fn a_closed_tetrahedron_validates() {
    assert_eq!(tetrahedron().validate(), Ok(()));
}

#[test]
fn a_short_face_tag_list_is_rejected() {
    let mut m = tetrahedron();
    m.solids[0].triangle_face.pop();
    assert_eq!(
        m.validate(),
        Err(ModelError::TriangleFaceLength {
            solid: SolidId(0),
            triangles: 4,
            tagged: 3
        })
    );
}

#[test]
fn a_vertex_index_out_of_range_is_rejected() {
    let mut m = tetrahedron();
    m.solids[0].triangles[2] = [0, 3, 9];
    assert_eq!(
        m.validate(),
        Err(ModelError::VertexOutOfRange {
            solid: SolidId(0),
            triangle: 2,
            vertex: 9,
            vertex_count: 4
        })
    );
}

#[test]
fn a_dangling_face_tag_is_rejected() {
    let mut m = tetrahedron();
    m.solids[0].triangle_face[1] = FaceId(7);
    assert_eq!(
        m.validate(),
        Err(ModelError::DanglingFace {
            solid: SolidId(0),
            triangle: 1,
            face: FaceId(7)
        })
    );
}

#[test]
fn a_face_of_a_missing_solid_is_rejected() {
    let mut m = tetrahedron();
    m.faces[3].solid = SolidId(1);
    assert_eq!(
        m.validate(),
        Err(ModelError::FaceSolidMissing {
            face: FaceId(3),
            solid: SolidId(1)
        })
    );
}

#[test]
fn a_face_of_another_solid_is_rejected() {
    let mut m = tetrahedron();
    let mut second = m.solids[0].clone();
    second.id = SolidId(1);
    m.solids.push(second);
    let error = m.validate();
    assert_eq!(
        error,
        Err(ModelError::FaceOfOtherSolid {
            solid: SolidId(1),
            triangle: 0,
            face: FaceId(0)
        })
    );
}

#[test]
fn a_non_finite_vertex_is_rejected() {
    let mut m = tetrahedron();
    m.solids[0].vertices[2][1] = f64::NAN;
    assert_eq!(
        m.validate(),
        Err(ModelError::NonFiniteVertex {
            solid: SolidId(0),
            vertex: 2
        })
    );
}

#[test]
fn ids_out_of_index_order_are_rejected() {
    let mut m = tetrahedron();
    m.faces.swap(0, 1);
    assert_eq!(
        m.validate(),
        Err(ModelError::FaceIdOrder {
            index: 0,
            id: FaceId(1)
        })
    );
    let mut m = tetrahedron();
    m.solids[0].id = SolidId(3);
    assert!(matches!(m.validate(), Err(ModelError::SolidIdOrder { .. })));
}

#[test]
fn a_non_positive_chord_tolerance_is_rejected() {
    let mut m = tetrahedron();
    m.chord_tolerance = 0.0;
    assert_eq!(m.validate(), Err(ModelError::ChordTolerance(0.0)));
}
