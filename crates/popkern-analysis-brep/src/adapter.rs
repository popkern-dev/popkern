//! BRep → [`AnalysisModel`].
//!
//! Solids are numbered in input order (each BRep's solids in its slot
//! order); faces are numbered by solid, each solid's faces in the
//! BRep's face slot order. The tessellation is welded per solid, so a
//! solid's mesh is closed when its BRep is, and checked for open edges,
//! positive volume, and planar faces lying on their planes.

use std::collections::BTreeMap;

use popkern_analysis::{
    AnalysisModel, FaceId, FaceInfo, Provenance, RoleSource, SolidId, SolidMesh, Surface, Units,
};
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::{BRep, FaceHandle, SolidHandle};
use popkern_mesh::TessellationOptions;

use crate::features::{FeatureTolerances, detect_features};
use crate::orient::orient_faces;
use crate::weld::{WeldedMesh, weld};

/// Why a BRep could not become an analysis model.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AdapterError {
    /// No input BRep has a solid.
    #[error("no solid to analyse")]
    NoSolid,
    /// A face's shell or the shell's solid is missing.
    #[error("BRep {brep}: a face has no shell or solid")]
    OrphanFace {
        /// Index of the BRep in the input.
        brep: usize,
    },
    /// Tessellation failed.
    #[error("tessellation failed: {detail}")]
    Mesh {
        /// The mesher's message.
        detail: String,
    },
    /// A measurement (face area, volume) failed.
    #[error("measurement failed: {detail}")]
    Measure {
        /// The measurer's message.
        detail: String,
    },
    /// A solid's welded mesh has open edges.
    #[error("solid {solid:?} is not watertight: {open_edges} open edges")]
    NotWatertight {
        /// The solid.
        solid: SolidId,
        /// Open edges after welding.
        open_edges: usize,
    },
    /// More solids, faces or vertices than fit a `u32` id.
    #[error("too many entities for u32 ids")]
    TooLarge,
    /// A solid's mesh encloses non-positive volume (inward shell).
    #[error("solid {solid:?} is inside out")]
    InsideOut {
        /// The solid.
        solid: SolidId,
    },
    /// A face has no triangles after welding.
    #[error("face {face:?} has no triangles")]
    UntessellatedFace {
        /// The face.
        face: FaceId,
    },
    /// A planar face's triangles do not lie on its plane.
    #[error("face {face:?}: projected area {projected} differs from area {area}")]
    OffPlane {
        /// The face.
        face: FaceId,
        /// |n · Σ n_tri| / 2.
        projected: f64,
        /// Summed triangle area.
        area: f64,
    },
}

/// A centroid's projection onto its face's surface counts for the fold
/// check within this many chord heights (a sound tessellation lies
/// within one).
const FOLD_PROJECTION_CHORDS: f64 = 4.0;

fn point(p: popkern_core::geometry::Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

/// The surface of a face. A plane's normal sense is provisional here:
/// `orient::orient_faces` sets it outward from the tessellation.
fn surface_kind(geometry: &SurfaceGeometry) -> Surface {
    match geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let n = [normal.x, normal.y, normal.z];
            Surface::Plane {
                normal: n,
                offset: n[0] * origin.x + n[1] * origin.y + n[2] * origin.z,
            }
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Surface::Cylinder {
            origin: point(*origin),
            axis: [axis.x, axis.y, axis.z],
            radius: *radius,
        },
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Surface::Cone {
            apex: point(*apex),
            axis: [axis.x, axis.y, axis.z],
            half_angle: *half_angle,
        },
        // The pole axis orients the chart only; the analysis sphere is
        // its centre and radius.
        SurfaceGeometry::Sphere { center, radius, .. } => Surface::Sphere {
            centre: point(*center),
            radius: *radius,
        },
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => Surface::Torus {
            centre: point(*center),
            axis: [axis.x, axis.y, axis.z],
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        },
        _ => Surface::Freeform,
    }
}

/// The solid each face of `brep` bounds, in face slot order.
fn face_solids(brep: &BRep, index: usize) -> Result<Vec<(FaceHandle, SolidHandle)>, AdapterError> {
    brep.faces
        .iter()
        .map(|(handle, face)| {
            let shell = brep
                .shells
                .get(face.shell)
                .ok_or(AdapterError::OrphanFace { brep: index })?;
            if !brep.solids.contains_key(shell.solid) {
                return Err(AdapterError::OrphanFace { brep: index });
            }
            Ok((handle, shell.solid))
        })
        .collect()
}

fn id(n: usize) -> Result<u32, AdapterError> {
    u32::try_from(n).map_err(|_| AdapterError::TooLarge)
}

/// Builds the analysis model of `breps`, tessellated with `options`
/// (popkern live's display tessellation keeps results aligned with the
/// display buffers). `chord_tolerance` in the model is
/// `options.chord_height`. Concave sharp edges and inferred fillet
/// roles come from [`detect_features`] with `features`.
pub fn model_from_breps(
    breps: &[BRep],
    options: &TessellationOptions,
    features: &FeatureTolerances,
    provenance: Provenance,
) -> Result<AnalysisModel, AdapterError> {
    model_from_breps_named(breps, options, features, provenance, &|_, _| None)
        .map(|named| named.model)
}

/// An analysis model built from BReps a host can name, with each face's
/// source (popkern live, PK-948).
#[derive(Debug, Clone, PartialEq)]
pub struct NamedModel {
    /// The analysis model; each face's `name` is what the namer gave.
    pub model: AnalysisModel,
    /// By `FaceId`: the index of the input BRep and the face's handle in
    /// it, so a host maps its own face picks to analysis faces.
    pub face_sources: Vec<(usize, FaceHandle)>,
}

/// [`model_from_breps`] with stable face names: `name(brep, face)`
/// spells the face's popkern name (or answers `None`), and setups then
/// reference the face by [`FaceRef::Name`](popkern_analysis::FaceRef),
/// which survives a parametric edit where a raw id does not.
pub fn model_from_breps_named(
    breps: &[BRep],
    options: &TessellationOptions,
    features: &FeatureTolerances,
    provenance: Provenance,
    name: &dyn Fn(usize, FaceHandle) -> Option<String>,
) -> Result<NamedModel, AdapterError> {
    let mut face_sources = Vec::new();
    let mut solids = Vec::new();
    let mut faces = Vec::new();
    let mut concave_edges = Vec::new();
    for (index, brep) in breps.iter().enumerate() {
        let owners = face_solids(brep, index)?;
        let mesh = popkern_mesh::tessellate(brep, options).map_err(|e| AdapterError::Mesh {
            detail: e.to_string(),
        })?;
        for (solid_handle, _) in &brep.solids {
            let solid_id = SolidId(id(solids.len())?);
            let first_face = faces.len();
            let mut face_ids: BTreeMap<FaceHandle, FaceId> = BTreeMap::new();
            for &(face_handle, owner) in &owners {
                if owner != solid_handle {
                    continue;
                }
                let face_id = FaceId(id(faces.len())?);
                face_ids.insert(face_handle, face_id);
                let mut info = face_info(brep, face_handle, face_id, solid_id)?;
                info.name = name(index, face_handle);
                faces.push(info);
                face_sources.push((index, face_handle));
            }
            let tagged: Vec<(usize, FaceId)> = mesh
                .face_map
                .iter()
                .enumerate()
                .filter_map(|(t, h)| face_ids.get(h).map(|&f| (t, f)))
                .collect();
            let WeldedMesh {
                vertices,
                triangles,
                triangle_face,
                open_edges,
            } = weld(&mesh, &tagged)?;
            if open_edges > 0 {
                return Err(AdapterError::NotWatertight {
                    solid: solid_id,
                    open_edges,
                });
            }
            orient_faces(
                solid_id,
                &vertices,
                &triangles,
                &triangle_face,
                &mut faces[first_face..],
            )?;
            // Folds of curved faces against popkern's surfaces; chord
            // error keeps a sound tessellation within a chord height.
            let curved: BTreeMap<FaceId, SurfaceGeometry> = face_ids
                .iter()
                .filter_map(|(handle, id)| {
                    let surface = brep.faces.get(*handle)?.surface.as_ref()?;
                    (!matches!(surface.geometry, SurfaceGeometry::Plane { .. }))
                        .then(|| (*id, surface.geometry.clone()))
                })
                .collect();
            let folded = crate::fold_check::folded_areas(
                &curved,
                &vertices,
                &triangles,
                &triangle_face,
                FOLD_PROJECTION_CHORDS * options.chord_height,
            );
            for (face, area) in folded {
                if let Some(info) = faces.get_mut(face.0 as usize) {
                    info.folded_area = Some(area);
                }
            }
            // Exact normals of freeform faces from popkern's surfaces.
            let geometry: BTreeMap<FaceId, SurfaceGeometry> = face_ids
                .iter()
                .filter_map(|(handle, id)| {
                    let surface = brep.faces.get(*handle)?.surface.as_ref()?;
                    matches!(surface_kind(&surface.geometry), Surface::Freeform)
                        .then(|| (*id, surface.geometry.clone()))
                })
                .collect();
            let freeform = |face: FaceId, p: [f64; 3]| -> Option<[f64; 3]> {
                let g = geometry.get(&face)?;
                let (u, v) = project_point_to_uv(g, &Point3::new(p[0], p[1], p[2]))?;
                let n = g.normal_at(u, v);
                crate::geometry::unit([n.x, n.y, n.z])
            };
            let found = detect_features(
                &vertices,
                &triangles,
                &triangle_face,
                &mut faces[first_face..],
                features,
                &freeform,
            );
            concave_edges.extend(found.concave_edges);
            let exact_volume = Some(solid_volume(brep, solid_handle, &owners)?);
            solids.push(SolidMesh {
                id: solid_id,
                name: None,
                vertices,
                triangles,
                triangle_face,
                exact_volume,
            });
        }
    }
    if solids.is_empty() {
        return Err(AdapterError::NoSolid);
    }
    Ok(NamedModel {
        model: AnalysisModel {
            units: Units::Millimetre,
            chord_tolerance: options.chord_height,
            solids,
            faces,
            concave_edges,
            provenance,
        },
        face_sources,
    })
}

/// Exact volume of one solid of `brep`: `body_volume` of a copy that
/// keeps only that solid's faces (popkern has no per-solid extraction
/// at the pinned commit). A single-solid BRep is measured directly.
fn solid_volume(
    brep: &BRep,
    solid: SolidHandle,
    owners: &[(FaceHandle, SolidHandle)],
) -> Result<f64, AdapterError> {
    let measure = |b: &BRep| {
        popkern_measure::body_volume(b, &b.tolerance).map_err(|e| AdapterError::Measure {
            detail: e.to_string(),
        })
    };
    if brep.solids.len() == 1 {
        return measure(brep);
    }
    let mut view = brep.clone();
    for &(face, owner) in owners {
        if owner != solid {
            view.faces.remove(face);
        }
    }
    view.solids.retain(|handle, _| handle == solid);
    measure(&view)
}

fn face_info(
    brep: &BRep,
    handle: FaceHandle,
    id: FaceId,
    solid: SolidId,
) -> Result<FaceInfo, AdapterError> {
    let face = brep
        .faces
        .get(handle)
        .ok_or(AdapterError::OrphanFace { brep: 0 })?;
    let surface = face
        .surface
        .as_ref()
        .map_or(Surface::Freeform, |s| surface_kind(&s.geometry));
    // The kernel measures curved analytic faces without holes in closed
    // form from their parameter ranges (planes from their loops'
    // vertices only, which misses curved edges, so not those).
    let curved = matches!(
        face.surface.as_ref().map(|s| &s.geometry),
        Some(
            SurfaceGeometry::Cylinder { .. }
                | SurfaceGeometry::Cone { .. }
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
        )
    );
    let exact_area = (curved && face.inner_loops.is_empty())
        .then(|| popkern_measure::face_area(brep, handle, &brep.tolerance).ok())
        .flatten();
    Ok(FaceInfo {
        id,
        solid,
        name: None,
        surface,
        role: None,
        role_source: RoleSource::Inferred,
        // set from the tessellation by orient::orient_faces
        area: 0.0,
        exact_area,
        // set from the tessellation by fold_check::folded_areas
        folded_area: None,
    })
}
