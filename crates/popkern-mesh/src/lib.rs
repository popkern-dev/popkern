#![warn(missing_docs)]

//! popkern-mesh: Tessellation and mesh export for popkern BRep geometry.
//!
//! This crate provides adaptive tessellation of BRep solids into triangle meshes,
//! and export to common mesh formats (STL, OBJ, glTF).

mod mesh;

pub use mesh::{FaceMesh, MeshVertex, TriMesh};

use popkern_core::topology::BRep;
use std::io::Write;
use thiserror::Error;

/// Options controlling the quality of adaptive tessellation.
#[derive(Debug, Clone)]
pub struct TessellationOptions {
    /// Maximum chord-height deviation allowed between the mesh and the true surface.
    ///
    /// Smaller values produce finer meshes. Default: `0.01`.
    pub chord_height: f64,
    /// Minimum number of subdivisions along each parametric direction per face.
    ///
    /// Default: `2`.
    pub min_subdivisions: usize,
    /// Maximum number of subdivisions along each parametric direction per face.
    ///
    /// Default: `64`.
    pub max_subdivisions: usize,
}

impl Default for TessellationOptions {
    fn default() -> Self {
        Self {
            chord_height: 0.01,
            min_subdivisions: 2,
            max_subdivisions: 64,
        }
    }
}

/// Errors that can occur during tessellation or mesh export.
#[derive(Debug, Error)]
pub enum MeshError {
    /// Tessellation of a face or solid failed.
    #[error("tessellation failed: {0}")]
    TessellationFailed(String),
    /// A face has no surface geometry and cannot be tessellated.
    #[error("face has no surface geometry")]
    NoSurfaceGeometry,
    /// Serialisation of the exported mesh format failed.
    #[error("export failed: {0}")]
    ExportFailed(String),
    /// An I/O error occurred while writing mesh output.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Tessellate every face of `brep` into a combined [`TriMesh`].
///
/// The `options` parameter controls chord-height tolerance and subdivision limits.
pub fn tessellate(_brep: &BRep, _options: &TessellationOptions) -> Result<TriMesh, MeshError> {
    todo!()
}

/// Export a [`TriMesh`] to the STL binary format, writing to `writer`.
pub fn export_stl<W: Write>(_mesh: &TriMesh, _writer: &mut W) -> Result<(), MeshError> {
    todo!()
}

/// Export a [`TriMesh`] to the Wavefront OBJ format, writing to `writer`.
pub fn export_obj<W: Write>(_mesh: &TriMesh, _writer: &mut W) -> Result<(), MeshError> {
    todo!()
}

/// Export a [`TriMesh`] to the glTF 2.0 format, writing to `writer`.
pub fn export_gltf<W: Write>(_mesh: &TriMesh, _writer: &mut W) -> Result<(), MeshError> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tessellation_options_default_values() {
        let opts = TessellationOptions::default();
        assert!((opts.chord_height - 0.01).abs() < f64::EPSILON);
        assert_eq!(opts.min_subdivisions, 2);
        assert_eq!(opts.max_subdivisions, 64);
    }
}
