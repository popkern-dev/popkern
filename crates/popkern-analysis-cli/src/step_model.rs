//! STEP file → [`AnalysisModel`].
//!
//! Imports with popkern-io at the pinned commit, walks the product
//! structure from its roots composing each instance's placement, and
//! builds the model with popkern live's display tessellation. Every
//! instance becomes its own placed solid, named after the instance.

use std::path::Path;

use popkern_analysis::{AnalysisModel, Provenance};
use popkern_analysis_brep::{AdapterError, FeatureTolerances, model_from_breps};
use popkern_core::geometry::transform::Transform3;
use popkern_core::topology::BRep;
use popkern_io::{StepDocument, StepImportOptions, import_step};
use popkern_mesh::TessellationOptions;
use sha2::{Digest, Sha256};

/// Deepest product nesting followed (a cycle in a malformed file would
/// otherwise recurse forever).
pub const MAX_ASSEMBLY_DEPTH: usize = 32;

/// Why a STEP file could not become an analysis model.
#[derive(Debug, thiserror::Error)]
pub enum StepModelError {
    /// The file could not be read.
    #[error("cannot read {path}: {detail}")]
    Read {
        /// The path.
        path: String,
        /// The OS message.
        detail: String,
    },
    /// popkern-io refused the file.
    #[error("STEP import failed: {0}")]
    Import(String),
    /// An instance is placed by a transform that is not a similarity
    /// (non-uniform scale or shear), which analytic surfaces cannot
    /// follow.
    #[error("instance `{child}` is placed by a non-rigid transform")]
    NonRigidPlacement {
        /// The instance name.
        child: String,
    },
    /// The product structure nests deeper than [`MAX_ASSEMBLY_DEPTH`]
    /// or refers to a missing product.
    #[error("malformed product structure at `{at}`")]
    MalformedStructure {
        /// Where it was found.
        at: String,
    },
    /// The adapter refused the geometry.
    #[error(transparent)]
    Adapter(#[from] AdapterError),
}

/// One placed body: a BRep in model coordinates and its instance name.
#[derive(Debug, Clone)]
pub struct PlacedBody {
    /// Instance path, e.g. `assembly/left`.
    pub name: String,
    /// The BRep, already placed.
    pub brep: BRep,
}

/// The tessellation the CLI uses: popkern live's display preset, so
/// per-vertex results line up with display buffers at integration.
pub fn tessellation_options() -> TessellationOptions {
    TessellationOptions::display()
}

/// Lowercase hex SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn visit(
    doc: &StepDocument,
    index: usize,
    world: Transform3,
    name: String,
    depth: usize,
    out: &mut Vec<PlacedBody>,
) -> Result<(), StepModelError> {
    let product = doc
        .products
        .get(index)
        .filter(|_| depth <= MAX_ASSEMBLY_DEPTH)
        .ok_or_else(|| StepModelError::MalformedStructure { at: name.clone() })?;
    let identity = world.matrix == Transform3::identity().matrix;
    if !identity && world.similarity_scale().is_none() {
        return Err(StepModelError::NonRigidPlacement { child: name });
    }
    for shape in &product.shapes {
        for brep in &shape.solids {
            let placed = if identity {
                brep.clone()
            } else {
                brep.transform(&world)
                    .map_err(|_| StepModelError::NonRigidPlacement {
                        child: name.clone(),
                    })?
            };
            out.push(PlacedBody {
                name: name.clone(),
                brep: placed,
            });
        }
    }
    for child in &product.children {
        let child_name = child.instance_name.clone().unwrap_or_else(|| {
            doc.products
                .get(child.product_index)
                .map_or_else(|| child.product_index.to_string(), |p| p.name.clone())
        });
        visit(
            doc,
            child.product_index,
            world.compose(&child.transform),
            format!("{name}/{child_name}"),
            depth + 1,
            out,
        )?;
    }
    Ok(())
}

/// The placed bodies of a STEP document, walking from its roots (or
/// every product, unplaced, when it declares none).
pub fn breps_from_document(doc: &StepDocument) -> Result<Vec<PlacedBody>, StepModelError> {
    let mut out = Vec::new();
    if doc.roots.is_empty() {
        for (i, p) in doc.products.iter().enumerate() {
            visit(doc, i, Transform3::identity(), p.name.clone(), 0, &mut out)?;
        }
    } else {
        for &root in &doc.roots {
            let name = doc
                .products
                .get(root)
                .map_or_else(|| root.to_string(), |p| p.name.clone());
            visit(doc, root, Transform3::identity(), name, 0, &mut out)?;
        }
    }
    Ok(out)
}

/// The placed bodies of a STEP file's bytes.
pub fn breps_from_step(bytes: &[u8]) -> Result<Vec<PlacedBody>, StepModelError> {
    let imported = import_step(bytes, &StepImportOptions::default())
        .map_err(|e| StepModelError::Import(e.to_string()))?;
    breps_from_document(&imported.document)
}

/// The analysis model of the STEP file at `path`, with each solid named
/// after its instance.
pub fn model_from_step_file(path: &Path) -> Result<AnalysisModel, StepModelError> {
    let bytes = std::fs::read(path).map_err(|e| StepModelError::Read {
        path: path.display().to_string(),
        detail: e.to_string(),
    })?;
    let bodies = breps_from_step(&bytes)?;
    let provenance = Provenance::Step {
        file_sha256: sha256_hex(&bytes),
    };
    let breps: Vec<BRep> = bodies.iter().map(|b| b.brep.clone()).collect();
    let mut model = model_from_breps(
        &breps,
        &tessellation_options(),
        &FeatureTolerances::default(),
        provenance,
    )?;
    let names = bodies
        .iter()
        .flat_map(|b| std::iter::repeat_n(b.name.clone(), b.brep.solids.len()));
    for (solid, name) in model.solids.iter_mut().zip(names) {
        solid.name = Some(name);
    }
    Ok(model)
}
