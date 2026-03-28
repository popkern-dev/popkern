#![warn(missing_docs)]

//! popkern-mesh: Tessellation and mesh export for popkern BRep geometry.
//!
//! This crate provides adaptive tessellation of BRep solids into triangle meshes,
//! and export to common mesh formats (STL, OBJ, glTF).

mod adaptive;
mod export;
mod mesh;
mod tessellate;
mod trim;
mod trim_clip;

pub use mesh::{FaceMesh, MeshVertex, TriMesh};

use popkern_core::topology::BRep;
use std::io::Write;
use std::path::Path;
use thiserror::Error;

/// Options controlling the quality of adaptive tessellation.
#[derive(Debug, Clone)]
pub struct TessellationOptions {
    /// Maximum chord-height deviation allowed between the mesh and the true surface.
    ///
    /// Smaller values produce finer meshes. Default: `0.005`.
    pub chord_height: f64,
    /// Minimum number of subdivisions along each parametric direction per face.
    ///
    /// Default: `4`.
    pub min_subdivisions: usize,
    /// Maximum number of subdivisions along each parametric direction per face.
    ///
    /// Default: `64`.
    pub max_subdivisions: usize,
}

impl Default for TessellationOptions {
    fn default() -> Self {
        Self {
            chord_height: 0.005,
            min_subdivisions: 4,
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
pub fn tessellate(brep: &BRep, options: &TessellationOptions) -> Result<TriMesh, MeshError> {
    tessellate::tessellate_brep(brep, options)
}

/// Export a [`TriMesh`] to an STL file at `path`.
///
/// When `binary` is `true` the file is written in STL binary format (80-byte
/// header + per-triangle 50-byte records).  When `false` the file is written
/// in human-readable STL ASCII format.
pub fn export_stl(mesh: &TriMesh, path: &Path, binary: bool) -> Result<(), MeshError> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    if binary {
        export::stl::write_stl_binary(mesh, &mut writer)
    } else {
        export::stl::write_stl_ascii(mesh, &mut writer)
    }
}

/// Export a [`TriMesh`] to the STL format, writing to `writer`.
///
/// When `binary` is `true` writes STL binary; otherwise writes STL ASCII.
pub fn export_stl_writer<W: Write>(
    mesh: &TriMesh,
    writer: &mut W,
    binary: bool,
) -> Result<(), MeshError> {
    if binary {
        export::stl::write_stl_binary(mesh, writer)
    } else {
        export::stl::write_stl_ascii(mesh, writer)
    }
}

/// Export a [`TriMesh`] to the Wavefront OBJ format, writing to `writer`.
///
/// Writes a Wavefront OBJ text stream containing all vertices, normals, UV
/// coordinates, and per-face triangle groups.  Face groups are named `face_0`,
/// `face_1`, etc., matching the order of [`TriMesh::face_meshes`].  All indices
/// in the output are 1-based as required by the OBJ specification.
pub fn export_obj<W: Write>(mesh: &TriMesh, writer: &mut W) -> Result<(), MeshError> {
    export::obj::write_obj(mesh, writer)
}

/// Export a [`TriMesh`] to glTF 2.0 format (JSON manifest + binary buffer).
///
/// `path` is the stem — `.gltf` and `.bin` extensions are appended automatically.
/// For example, passing `"model"` produces `model.gltf` and `model.bin`.
pub fn export_gltf(mesh: &TriMesh, path: &Path) -> Result<(), MeshError> {
    let stem = path.to_string_lossy();
    let gltf_path_str = format!("{stem}.gltf");
    let bin_path_str = format!("{stem}.bin");
    let bin_filename = std::path::Path::new(&bin_path_str)
        .file_name()
        .unwrap_or_default()
        .to_string_lossy();
    let mut gltf_file = std::fs::File::create(&gltf_path_str)?;
    let mut bin_file = std::fs::File::create(&bin_path_str)?;
    export::gltf::write_gltf(mesh, &mut gltf_file, &mut bin_file, &bin_filename)
}

#[cfg(test)]
mod integration_tests;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tessellation_options_default_values() {
        let opts = TessellationOptions::default();
        assert!((opts.chord_height - 0.005).abs() < f64::EPSILON);
        assert_eq!(opts.min_subdivisions, 4);
        assert_eq!(opts.max_subdivisions, 64);
    }
}
