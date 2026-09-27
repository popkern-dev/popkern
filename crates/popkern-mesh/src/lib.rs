#![warn(missing_docs)]

//! popkern-mesh: Tessellation and mesh export for popkern BRep geometry.
//!
//! This crate provides adaptive tessellation of BRep solids into triangle meshes,
//! and export to common mesh formats (STL, OBJ, glTF).

mod adaptive;
mod angular_refine;
mod boundary_first;
mod boundary_sew;
mod crack_zipper;
mod error;
mod export;
mod face_pass;
mod facet_angle;
mod mesh;
pub mod offset_oracle;
pub mod point_query;
mod tessellate;
mod trim;
mod trim_canonical;
mod trim_clip;
mod trim_sampling;
mod vertex_normal;
mod winding_side;

pub use error::MeshError;
pub use mesh::{FaceMesh, MeshVertex, TriMesh};
pub use tessellate::tessellate_face;

use popkern_core::topology::BRep;
use std::io::Write;
use std::path::Path;

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
    /// Maximum angle (radians) the surface normal / curve tangent may
    /// turn across one facet — the angular deflection bound that
    /// complements the absolute `chord_height`.
    ///
    /// Chord height alone lets small-radius features (fillets, rounds,
    /// thin bosses) facet coarsely: a 0.04-radius round meets a 0.005
    /// chord tolerance with 45° facets, and reflection-line (zebra)
    /// or shaded views then show interpolated-normal wobble the true
    /// surface does not have. Values `>= π` disable the bound (the
    /// facet count is then governed by chord height alone).
    ///
    /// Default: `π` (off — the analysis-gate meshes stay chord-only).
    /// [`TessellationOptions::display`] turns it on.
    pub max_facet_angle: f64,
}

impl Default for TessellationOptions {
    fn default() -> Self {
        Self {
            chord_height: 0.005,
            min_subdivisions: 4,
            max_subdivisions: 64,
            max_facet_angle: std::f64::consts::PI,
        }
    }
}

impl TessellationOptions {
    /// Display-quality preset for viewers and surface analysis
    /// (shaded, zebra, curvature): the default chord height plus a
    /// 12° angular facet bound, so small-radius blends mesh with the
    /// same normal fidelity as large ones.
    pub fn display() -> Self {
        Self {
            max_facet_angle: DISPLAY_MAX_FACET_ANGLE,
            ..Self::default()
        }
    }

    /// Whether the angular facet bound constrains anything (finite,
    /// positive, below π).
    pub(crate) fn angular_bound(&self) -> Option<f64> {
        (self.max_facet_angle > 0.0 && self.max_facet_angle < std::f64::consts::PI)
            .then_some(self.max_facet_angle)
    }
}

/// The angular facet bound of [`TessellationOptions::display`]: 12°.
/// Reflection lines sampled at this density stay visually straight on
/// G2 blends at the playground's 16-stripe zebra frequency while a
/// quarter-round costs 8 facets.
pub const DISPLAY_MAX_FACET_ANGLE: f64 = 12.0 * std::f64::consts::PI / 180.0;

/// Tessellate every face of `brep` into a combined [`TriMesh`].
///
/// The `options` parameter controls chord-height tolerance and subdivision limits.
pub fn tessellate(brep: &BRep, options: &TessellationOptions) -> Result<TriMesh, MeshError> {
    tessellate::tessellate_brep(
        brep,
        options,
        &mut popkern_core::progress::StageTracker::inert(),
    )
}

/// [`tessellate`] with PK-276 interactive control (cooperative
/// cancellation + per-round/stage progress).
///
/// With a cancel token the tessellation runs detached on a worker
/// thread over a private clone of `brep`
/// ([`popkern_core::cancellation::run_detach_on_cancel_with_events`]):
/// the call returns [`MeshError::Cancelled`] within milliseconds of the
/// token firing, and the abandoned worker exits at its next round or
/// stage poll (per-face closures also short-circuit on a fired token).
/// Without a token the pipeline runs inline. A call with a never-fired
/// token and a null reporter produces a mesh identical to
/// [`tessellate`]'s.
pub fn tessellate_with_control(
    brep: &BRep,
    options: &TessellationOptions,
    control: popkern_core::progress::OperationControl<'_>,
) -> Result<TriMesh, MeshError> {
    use popkern_core::cancellation::{DetachOutcome, run_detach_on_cancel_with_events};
    use popkern_core::progress::StageTracker;
    let popkern_core::progress::OperationControl {
        cancel,
        mut reporter,
    } = control;
    match cancel {
        None => tessellate::tessellate_brep(
            brep,
            options,
            &mut StageTracker::new(None, reporter, tessellate::TESSELLATE_STAGE_CEILING),
        ),
        Some(token) => {
            let (worker_brep, worker_options, worker_token) =
                (brep.clone(), options.clone(), token.clone());
            match run_detach_on_cancel_with_events(
                token,
                move |mut sink| {
                    let mut tracker = StageTracker::new(
                        Some(&worker_token),
                        Some(&mut sink),
                        tessellate::TESSELLATE_STAGE_CEILING,
                    );
                    tessellate::tessellate_brep(&worker_brep, &worker_options, &mut tracker)
                },
                |event| {
                    if let Some(r) = reporter.as_mut() {
                        r.on_stage(event);
                    }
                },
            ) {
                DetachOutcome::Completed(result) => result,
                DetachOutcome::Detached => Err(MeshError::Cancelled),
            }
        }
    }
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
