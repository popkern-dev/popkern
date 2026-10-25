//! A solid cut to a box as a model of its own: the domain of a
//! submodel, structural or thermal.
//!
//! The solid's faces that reach into the box keep their metadata (each
//! becomes one face of the cut model); each box face that cuts the solid
//! adds a planar cap face. The caller attaches conditions by the parent
//! face of each cut-model face, or by cap.

use super::box_cut::{CutError, Origin, cut_to_box};
use super::{AnalysisModel, FaceId, FaceInfo, RoleSource, SolidId, SolidMesh, Surface};

/// A box face: its normal axis and whether it is the upper one.
pub type BoxFace = (usize, bool);

/// A solid cut to a box, as a one-solid model.
#[derive(Debug, Clone, PartialEq)]
pub struct CutModel {
    /// The cut solid (its faces first the parent's that reach into the
    /// box, in order of first appearance, then the caps).
    pub model: AnalysisModel,
    /// The cut-model face of each parent face, by parent `FaceId`.
    pub face_map: Vec<Option<FaceId>>,
    /// The parent face of each cut-model face (`None` for caps).
    pub parent_face: Vec<Option<FaceId>>,
    /// Box faces that cut the solid, with their cap's face.
    pub caps: Vec<(BoxFace, FaceId)>,
    /// The box cut to (each bound moved off mesh vertices).
    pub lo: [f64; 3],
    /// Upper corner of that box.
    pub hi: [f64; 3],
}

impl CutModel {
    /// Whether `p` lies at least `margin` from every cap's plane.
    pub fn clear_of_caps(&self, p: [f64; 3], margin: f64) -> bool {
        clear_of_caps(&self.caps, self.lo, self.hi, p, margin)
    }
}

/// Whether `p` lies at least `margin` from the plane of every cap in
/// `caps` of the box [`lo`, `hi`].
pub fn clear_of_caps(
    caps: &[(BoxFace, FaceId)],
    lo: [f64; 3],
    hi: [f64; 3],
    p: [f64; 3],
    margin: f64,
) -> bool {
    caps.iter().all(|((axis, upper), _)| {
        let plane = if *upper { hi[*axis] } else { lo[*axis] };
        (p[*axis] - plane).abs() >= margin
    })
}

fn triangle_area(v: &[[f64; 3]], t: [u32; 3]) -> f64 {
    let [a, b, c] = t.map(|i| v[i as usize]);
    let (u, w) = (
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    let n = [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ];
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Solid `solid` of `model` cut to the box [`lo`, `hi`].
pub fn cut_model(
    model: &AnalysisModel,
    solid: SolidId,
    lo: [f64; 3],
    hi: [f64; 3],
) -> Result<CutModel, CutError> {
    let mesh = model.solids.get(solid.0 as usize).ok_or(CutError::Empty)?;
    let cut = cut_to_box(mesh, lo, hi)?;
    let mut face_map: Vec<Option<FaceId>> = vec![None; model.faces.len()];
    let mut parent_face: Vec<Option<FaceId>> = Vec::new();
    let mut faces: Vec<FaceInfo> = Vec::new();
    let mut caps: Vec<(BoxFace, FaceId)> = Vec::new();
    let mut triangle_face = Vec::with_capacity(cut.triangles.len());
    for o in &cut.origin {
        let id = match *o {
            Origin::Face(f) => match face_map.get(f.0 as usize).copied().flatten() {
                Some(id) => id,
                None => {
                    let id = FaceId(faces.len() as u32);
                    let info = model.faces.get(f.0 as usize).ok_or(CutError::Empty)?;
                    faces.push(FaceInfo {
                        id,
                        solid: SolidId(0),
                        area: 0.0,
                        ..info.clone()
                    });
                    parent_face.push(Some(f));
                    if let Some(slot) = face_map.get_mut(f.0 as usize) {
                        *slot = Some(id);
                    }
                    id
                }
            },
            Origin::Cap { axis, upper } => match caps.iter().find(|c| c.0 == (axis, upper)) {
                Some(c) => c.1,
                None => {
                    let id = FaceId(faces.len() as u32);
                    let mut normal = [0.0; 3];
                    normal[axis] = if upper { 1.0 } else { -1.0 };
                    let at = if upper { cut.hi[axis] } else { cut.lo[axis] };
                    faces.push(FaceInfo {
                        id,
                        solid: SolidId(0),
                        name: None,
                        surface: Surface::Plane {
                            normal,
                            offset: normal[axis] * at,
                        },
                        role: None,
                        role_source: RoleSource::Inferred,
                        area: 0.0,
                        exact_area: None,
                        folded_area: None,
                    });
                    parent_face.push(None);
                    caps.push(((axis, upper), id));
                    id
                }
            },
        };
        triangle_face.push(id);
    }
    for (t, f) in cut.triangles.iter().zip(&triangle_face) {
        if let Some(info) = faces.get_mut(f.0 as usize) {
            info.area += triangle_area(&cut.vertices, *t);
        }
    }
    let sub = AnalysisModel {
        units: model.units,
        chord_tolerance: model.chord_tolerance,
        solids: vec![SolidMesh {
            id: SolidId(0),
            name: mesh.name.clone(),
            vertices: cut.vertices,
            triangles: cut.triangles,
            triangle_face,
            exact_volume: None,
        }],
        faces,
        concave_edges: Vec::new(),
        provenance: model.provenance.clone(),
    };
    Ok(CutModel {
        model: sub,
        face_map,
        parent_face,
        caps,
        lo: cut.lo,
        hi: cut.hi,
    })
}
