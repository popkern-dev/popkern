//! Test models of several solids.

use crate::mesh_test_helpers::model_of;
use crate::model::{FaceId, SolidId, SolidMesh};

/// A model of several solids, each with `faces` faces: solid `s` gets
/// `SolidId(s)` and faces `faces·s ..`, its meshes' face tags offset to
/// match.
pub(crate) fn model_of_many(meshes: Vec<SolidMesh>, faces: usize) -> crate::model::AnalysisModel {
    let mut model = model_of(meshes[0].clone(), faces);
    model.solids.clear();
    model.faces.clear();
    for (s, mut mesh) in meshes.into_iter().enumerate() {
        let mut one = model_of(mesh.clone(), faces);
        mesh.id = SolidId(s as u32);
        for f in &mut mesh.triangle_face {
            f.0 += (faces * s) as u32;
        }
        for (f, info) in one.faces.iter_mut().enumerate() {
            info.id = FaceId((faces * s + f) as u32);
            info.solid = SolidId(s as u32);
        }
        model.faces.extend(one.faces);
        model.solids.push(mesh);
    }
    model
}
