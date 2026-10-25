use popkern_analysis::{FaceId, FaceInfo, RoleSource, SolidId, Surface};

use crate::adapter::AdapterError;
use crate::orient::{orient_faces, signed_volume};
use crate::weld::open_edges;

/// Unit tetrahedron, outward, one face per triangle; face 0 is z = 0.
fn tetrahedron() -> (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<FaceId>) {
    (
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
        (0..4).map(FaceId).collect(),
    )
}

fn faces(bottom: Surface) -> Vec<FaceInfo> {
    (0..4)
        .map(|i| FaceInfo {
            id: FaceId(i),
            solid: SolidId(0),
            name: None,
            surface: if i == 0 {
                bottom.clone()
            } else {
                Surface::Freeform
            },
            role: None,
            role_source: RoleSource::Inferred,
            area: 0.0,
            exact_area: None,
            folded_area: None,
        })
        .collect()
}

#[test]
fn an_inward_plane_normal_is_turned_outward() {
    let (v, t, tf) = tetrahedron();
    assert!((signed_volume(&v, &t) - 1.0 / 6.0).abs() < 1e-15);
    assert_eq!(open_edges(&t), 0);
    let mut f = faces(Surface::Plane {
        normal: [0.0, 0.0, 1.0],
        offset: 0.0,
    });
    assert_eq!(orient_faces(SolidId(0), &v, &t, &tf, &mut f), Ok(()));
    assert_eq!(
        f[0].surface,
        Surface::Plane {
            normal: [0.0, 0.0, -1.0],
            offset: -0.0
        }
    );
    assert_eq!(f[0].area, 0.5);
}

#[test]
fn an_inside_out_mesh_is_declined() {
    let (v, mut t, tf) = tetrahedron();
    for tri in &mut t {
        tri.swap(1, 2);
    }
    assert_eq!(open_edges(&t), 0, "reversal keeps it closed");
    let mut f = faces(Surface::Freeform);
    assert_eq!(
        orient_faces(SolidId(0), &v, &t, &tf, &mut f),
        Err(AdapterError::InsideOut { solid: SolidId(0) })
    );
}

#[test]
fn triangles_off_their_plane_are_declined() {
    let (v, t, tf) = tetrahedron();
    // face 0 lies in z = 0; claim it is the plane x = 0
    let mut f = faces(Surface::Plane {
        normal: [1.0, 0.0, 0.0],
        offset: 0.0,
    });
    assert!(matches!(
        orient_faces(SolidId(0), &v, &t, &tf, &mut f),
        Err(AdapterError::OffPlane {
            face: FaceId(0),
            ..
        })
    ));
}

#[test]
fn a_face_without_triangles_is_declined() {
    let (v, t, mut tf) = tetrahedron();
    tf[3] = FaceId(2);
    let mut f = faces(Surface::Freeform);
    assert_eq!(
        orient_faces(SolidId(0), &v, &t, &tf, &mut f),
        Err(AdapterError::UntessellatedFace { face: FaceId(3) })
    );
}

#[test]
fn open_and_misoriented_edges_are_counted() {
    let (_, t, _) = tetrahedron();
    assert_eq!(open_edges(&t[..3]), 3);
    let mut flipped = t.clone();
    flipped[3].swap(1, 2);
    assert_eq!(
        open_edges(&flipped),
        6,
        "one flipped triangle unbalances 3 edges twice"
    );
}
