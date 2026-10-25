//! The analysis model: the solver's only geometric input.
//!
//! Closed, outward-oriented triangle meshes per solid, each triangle
//! tagged with the face it came from, plus per-face metadata. Every
//! popkern-specific advantage has an optional slot: the STEP path fills
//! what it can, and the popkern path fills it exactly. See
//! `docs/architecture/analysis-library-boundary.md`.

pub mod box_cut;
mod cap_fill;
pub mod cut_model;
pub mod tessellation;

use serde::Serialize;

/// Identifies a solid; equal to its index in [`AnalysisModel::solids`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct SolidId(pub u32);

/// Identifies a face; equal to its index in [`AnalysisModel::faces`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct FaceId(pub u32);

/// Length unit of the model's coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Units {
    /// Millimetres (the only unit at Stage 1).
    Millimetre,
}

/// Where the geometry came from, recorded in every report.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub enum Provenance {
    /// Imported from a STEP file.
    Step {
        /// SHA-256 of the file bytes, lowercase hex.
        file_sha256: String,
    },
    /// Built in memory by popkern.
    Popkern {
        /// popkern document identifier.
        document: String,
        /// Document revision the geometry was taken from.
        revision: u64,
    },
}

/// The geometry the solver analyses.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisModel {
    /// Unit of every coordinate and length.
    pub units: Units,
    /// Chord tolerance of the tessellation, recorded in reports.
    pub chord_tolerance: f64,
    /// Solids, indexed by [`SolidId`].
    pub solids: Vec<SolidMesh>,
    /// Faces of all solids, indexed by [`FaceId`].
    pub faces: Vec<FaceInfo>,
    /// Concave sharp edges (stress raisers), from BRep dihedral angles.
    pub concave_edges: Vec<ConcaveEdge>,
    /// Source of the geometry.
    pub provenance: Provenance,
}

/// One solid as a closed, outward-oriented triangle mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidMesh {
    /// This solid's id (its index in the model).
    pub id: SolidId,
    /// Body name, when the source has one.
    pub name: Option<String>,
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles as vertex indices, counter-clockwise seen from outside.
    pub triangles: Vec<[u32; 3]>,
    /// Face of each triangle; same length as `triangles`.
    pub triangle_face: Vec<FaceId>,
    /// The source kernel's volume of the solid, when known: exact for
    /// all-planar solids, a fine tessellation otherwise (popkern-measure;
    /// within 0.6% of OpenCASCADE on the gallery, PK-807).
    pub exact_volume: Option<f64>,
}

/// Metadata of one BRep face.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceInfo {
    /// This face's id (its index in the model).
    pub id: FaceId,
    /// The solid this face bounds.
    pub solid: SolidId,
    /// Stable popkern face name; `None` from STEP.
    pub name: Option<String>,
    /// The face's underlying surface.
    pub surface: Surface,
    /// Design role of the face, when known.
    pub role: Option<FaceRole>,
    /// Whether `role` is exact or inferred.
    pub role_source: RoleSource,
    /// Face area in model units squared.
    pub area: f64,
    /// The kernel's closed-form area of a curved analytic face without
    /// holes, over the face's parameter box: exact for a rectangular
    /// domain, too large otherwise, so an upper bound. Checks the
    /// tessellation (`tessellation::mismatched_faces`).
    pub exact_area: Option<f64>,
    /// Area of the face's tessellation folded against its surface (the
    /// minority of triangles whose normal opposes the kernel's exact
    /// normal), where the kernel can tell (curved faces).
    pub folded_area: Option<f64>,
}

/// The surface a face lies on, exact from the BRep.
#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    /// The plane `normal · x = offset`, with unit `normal`.
    Plane {
        /// Unit normal.
        normal: [f64; 3],
        /// Signed distance of the plane from the origin along `normal`.
        offset: f64,
    },
    /// A circular cylinder.
    Cylinder {
        /// A point on the axis.
        origin: [f64; 3],
        /// Unit axis direction.
        axis: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// A circular cone.
    Cone {
        /// Apex.
        apex: [f64; 3],
        /// Unit axis direction.
        axis: [f64; 3],
        /// Half-angle in radians.
        half_angle: f64,
    },
    /// A sphere.
    Sphere {
        /// Centre.
        centre: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// A torus.
    Torus {
        /// Centre.
        centre: [f64; 3],
        /// Unit axis direction.
        axis: [f64; 3],
        /// Distance from the centre to the tube centre.
        major_radius: f64,
        /// Tube radius.
        minor_radius: f64,
    },
    /// Any other surface (B-spline, swept, offset, ...).
    Freeform,
}

/// Design role of a face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FaceRole {
    /// A fillet (blend) of the given radius.
    Fillet {
        /// Fillet radius.
        radius: f64,
    },
    /// A chamfer.
    Chamfer,
}

/// Whether a face role is known exactly or inferred from geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum RoleSource {
    /// From popkern construction history.
    Exact,
    /// Inferred from surface kind and tangency (STEP path).
    Inferred,
}

/// A concave sharp edge between two faces.
#[derive(Debug, Clone, PartialEq)]
pub struct ConcaveEdge {
    /// The two faces meeting at the edge.
    pub faces: (FaceId, FaceId),
    /// Dihedral angle in radians.
    pub dihedral: f64,
    /// The edge as a polyline.
    pub polyline: Vec<[f64; 3]>,
}

/// Why an [`AnalysisModel`] is malformed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ModelError {
    /// The chord tolerance is not a positive finite number.
    #[error("chord tolerance {0} is not positive and finite")]
    ChordTolerance(f64),
    /// A solid's id differs from its index.
    #[error("solid at index {index} has id {id:?}")]
    SolidIdOrder {
        /// Position in `solids`.
        index: usize,
        /// The id found there.
        id: SolidId,
    },
    /// A face's id differs from its index.
    #[error("face at index {index} has id {id:?}")]
    FaceIdOrder {
        /// Position in `faces`.
        index: usize,
        /// The id found there.
        id: FaceId,
    },
    /// A face bounds a solid that does not exist.
    #[error("face {face:?} refers to missing solid {solid:?}")]
    FaceSolidMissing {
        /// The face.
        face: FaceId,
        /// The missing solid.
        solid: SolidId,
    },
    /// `triangle_face` and `triangles` differ in length.
    #[error("solid {solid:?} has {triangles} triangles but {tagged} face tags")]
    TriangleFaceLength {
        /// The solid.
        solid: SolidId,
        /// Number of triangles.
        triangles: usize,
        /// Number of face tags.
        tagged: usize,
    },
    /// A triangle uses a vertex index past the end of `vertices`.
    #[error("solid {solid:?} triangle {triangle} uses vertex {vertex} of {vertex_count}")]
    VertexOutOfRange {
        /// The solid.
        solid: SolidId,
        /// Triangle index.
        triangle: usize,
        /// The bad vertex index.
        vertex: u32,
        /// Number of vertices.
        vertex_count: usize,
    },
    /// A triangle is tagged with a face that does not exist.
    #[error("solid {solid:?} triangle {triangle} is tagged with missing face {face:?}")]
    DanglingFace {
        /// The solid.
        solid: SolidId,
        /// Triangle index.
        triangle: usize,
        /// The missing face.
        face: FaceId,
    },
    /// A triangle is tagged with a face of another solid.
    #[error("solid {solid:?} triangle {triangle} is tagged with face {face:?} of another solid")]
    FaceOfOtherSolid {
        /// The solid.
        solid: SolidId,
        /// Triangle index.
        triangle: usize,
        /// The face, which bounds a different solid.
        face: FaceId,
    },
    /// A vertex coordinate is NaN or infinite.
    #[error("solid {solid:?} vertex {vertex} is not finite")]
    NonFiniteVertex {
        /// The solid.
        solid: SolidId,
        /// Vertex index.
        vertex: usize,
    },
}

impl SolidMesh {
    /// Signed volume enclosed by the mesh, Σ a·(b×c)/6 in triangle
    /// order: positive for a closed outward mesh.
    pub fn enclosed_volume(&self) -> f64 {
        let mut six = 0.0;
        for t in &self.triangles {
            let [a, b, c] = t.map(|i| self.vertices.get(i as usize).copied().unwrap_or([0.0; 3]));
            six += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
        six / 6.0
    }
}

impl AnalysisModel {
    /// Checks the model's structural invariants: ids equal indices,
    /// every triangle indexes real vertices and is tagged with a face
    /// of its own solid, and every coordinate is finite. Closedness and
    /// orientation are geometric checks made by the rasteriser.
    pub fn validate(&self) -> Result<(), ModelError> {
        if !(self.chord_tolerance.is_finite() && self.chord_tolerance > 0.0) {
            return Err(ModelError::ChordTolerance(self.chord_tolerance));
        }
        for (index, solid) in self.solids.iter().enumerate() {
            if solid.id.0 as usize != index {
                return Err(ModelError::SolidIdOrder {
                    index,
                    id: solid.id,
                });
            }
        }
        for (index, face) in self.faces.iter().enumerate() {
            if face.id.0 as usize != index {
                return Err(ModelError::FaceIdOrder { index, id: face.id });
            }
            if face.solid.0 as usize >= self.solids.len() {
                return Err(ModelError::FaceSolidMissing {
                    face: face.id,
                    solid: face.solid,
                });
            }
        }
        for solid in &self.solids {
            self.validate_solid(solid)?;
        }
        Ok(())
    }

    fn validate_solid(&self, solid: &SolidMesh) -> Result<(), ModelError> {
        let id = solid.id;
        if solid.triangle_face.len() != solid.triangles.len() {
            return Err(ModelError::TriangleFaceLength {
                solid: id,
                triangles: solid.triangles.len(),
                tagged: solid.triangle_face.len(),
            });
        }
        if let Some(vertex) = solid
            .vertices
            .iter()
            .position(|v| v.iter().any(|c| !c.is_finite()))
        {
            return Err(ModelError::NonFiniteVertex { solid: id, vertex });
        }
        let vertex_count = solid.vertices.len();
        for (triangle, (corners, &face)) in
            solid.triangles.iter().zip(&solid.triangle_face).enumerate()
        {
            if let Some(&vertex) = corners.iter().find(|&&v| v as usize >= vertex_count) {
                return Err(ModelError::VertexOutOfRange {
                    solid: id,
                    triangle,
                    vertex,
                    vertex_count,
                });
            }
            match self.faces.get(face.0 as usize) {
                None => {
                    return Err(ModelError::DanglingFace {
                        solid: id,
                        triangle,
                        face,
                    });
                }
                Some(info) if info.solid != id => {
                    return Err(ModelError::FaceOfOtherSolid {
                        solid: id,
                        triangle,
                        face,
                    });
                }
                Some(_) => {}
            }
        }
        Ok(())
    }
}
