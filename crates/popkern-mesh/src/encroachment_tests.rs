//! A face subset meshes by the whole-body pipeline (PK-900).

use popkern_core::geometry::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{TessellationOptions, tessellate, tessellate_face, tessellate_faces};

#[test]
fn a_face_subset_meshes_like_its_body_not_like_a_uniform_grid() {
    let body = BRep::cylinder(10.0, 10.0, Tolerance::default()).expect("cylinder");
    // A certification-level chord: fine, with a deep subdivision cap.
    let options = TessellationOptions {
        chord_height: 1e-3,
        max_subdivisions: 4096,
        ..TessellationOptions::default()
    };
    let wall = body
        .faces
        .keys()
        .find(|&f| {
            matches!(
                body.faces[f].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .expect("wall");
    let subset = tessellate_faces(&body, &[wall], &options).expect("subset");
    assert_eq!(subset.len(), 1);
    assert_eq!(subset[0].face, wall);
    let whole = tessellate(&body, &options).expect("whole");
    let in_whole = whole.face_map.iter().filter(|&&f| f == wall).count();
    let alone = subset[0].triangles.len();
    assert!(
        alone > 0 && alone <= in_whole + in_whole / 2,
        "subset {alone} vs whole body's wall {in_whole}"
    );
    // The subset honours the chord: every facet centroid sags at most
    // the chord height below the wall (radius 10 about z).
    for t in &subset[0].triangles {
        let c = t.iter().fold((0.0, 0.0), |(x, y), &i| {
            let p = subset[0].vertices[i].position;
            (x + p.x / 3.0, y + p.y / 3.0)
        });
        let sag = 10.0 - (c.0 * c.0 + c.1 * c.1).sqrt();
        assert!((-1e-9..=1.5e-3).contains(&sag), "sag {sag}");
    }
    // The single-face mesher is the uniform grid at this cap: the
    // reason a region cannot be meshed face by face with it.
    let grid = tessellate_face(&body, wall, &options).expect("grid");
    assert!(
        grid.triangles.len() >= 10 * alone,
        "grid {} vs subset {alone}",
        grid.triangles.len()
    );
}

#[test]
fn a_whole_body_subset_meshes_every_face() {
    let body = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    let faces: Vec<_> = body.faces.keys().collect();
    let options = TessellationOptions::default();
    let meshes = tessellate_faces(&body, &faces, &options).expect("subset");
    assert_eq!(meshes.len(), faces.len());
    assert!(meshes.iter().all(|m| !m.triangles.is_empty()));
}
