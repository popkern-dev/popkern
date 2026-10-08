//! STEP export of one [`Part`] or of several as one document: the path
//! form, the in-memory bytes form (hosts without a filesystem — the
//! popkern live — hand the bytes to a download), and the
//! round-trip gate.
//!
//! Several parts export as ONE product whose shape representation
//! holds every part's solids: a multi-body model (a gearbox's gears, a
//! bearing's races and balls) leaves as separate `MANIFOLD_SOLID_BREP`s
//! side by side, never merged and never reduced to one body. Product
//! structure (a product per part, placed instances) is the assembly
//! layer's job (PK-629).

use popkern_core::topology::BRep;

use crate::error::PartError;
use crate::part::Part;

impl Part {
    /// Export the part as a STEP solid.
    ///
    /// # Errors
    ///
    /// [`PartError::Io`] when export fails.
    pub fn export_step(&self, path: impl AsRef<std::path::Path>) -> Result<(), PartError> {
        export_step_of(std::slice::from_ref(self), path)
    }

    /// The STEP solid as Part 21 bytes, with `name` as the product
    /// name (the file stem when written to disk). This is the in-memory
    /// form [`Part::export_step`] writes: a host without a filesystem
    /// (popkern live) hands these bytes to a download.
    ///
    /// # Errors
    ///
    /// [`PartError::Io`] when the STEP writer declines the body.
    pub fn step_bytes(&self, name: &str) -> Result<Vec<u8>, PartError> {
        step_bytes_of(std::slice::from_ref(self), name)
    }

    /// Round-trip gate: re-import the STEP just written at `path`
    /// and require the full face count back (the dropped-faces
    /// defect class). Returns `self` for chaining.
    ///
    /// # Errors
    ///
    /// [`PartError::Io`] on import failure;
    /// [`PartError::NotApplicable`] describing the mismatch when
    /// faces are missing.
    pub fn assert_step_roundtrip(
        self,
        path: impl AsRef<std::path::Path>,
    ) -> Result<Part, PartError> {
        assert_step_roundtrip_of(std::slice::from_ref(&self), path)?;
        Ok(self)
    }
}

/// Every part's solids as one STEP document's Part 21 bytes, `name`
/// naming the product. One part exports exactly as
/// [`Part::step_bytes`].
///
/// # Errors
///
/// [`PartError::NotApplicable`] for no parts; [`PartError::Io`] when
/// the STEP writer declines a body.
pub fn step_bytes_of(parts: &[Part], name: &str) -> Result<Vec<u8>, PartError> {
    step_bytes_of_in(parts, name, popkern_io::RepresentationContext::default())
}

/// [`step_bytes_of`] under an explicit representation context: the
/// length unit the coordinates are declared in and the geometric
/// uncertainty. A document storing true SI lengths (a mm, cm or inch
/// model) declares metres; the default context declares millimetres,
/// the millimetre-magnitude convention of default-unit documents.
///
/// # Errors
///
/// As [`step_bytes_of`].
pub fn step_bytes_of_in(
    parts: &[Part],
    name: &str,
    context: popkern_io::RepresentationContext,
) -> Result<Vec<u8>, PartError> {
    step_bytes_of_styled(parts, name, context, None)
}

/// [`step_bytes_of_in`] with every solid coloured `colour` (RGB in
/// `[0, 1]`; `None` writes no presentation).
///
/// # Errors
///
/// As [`step_bytes_of_in`].
pub fn step_bytes_of_styled(
    parts: &[Part],
    name: &str,
    context: popkern_io::RepresentationContext,
    colour: Option<[f32; 3]>,
) -> Result<Vec<u8>, PartError> {
    if parts.is_empty() {
        return Err(PartError::NotApplicable(
            "STEP export of no parts: the document has no body".into(),
        ));
    }
    let solids: Vec<BRep> = parts.iter().map(|part| part.brep().clone()).collect();
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.to_string(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids,
                context: context.clone(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata {
                solid_color: colour,
                ..popkern_io::ProductMetadata::default()
            },
        }],
        roots: vec![0],
        // The writer declares units from the document context.
        units: popkern_io::UnitContext { context },
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    popkern_io::export_step(&doc, &popkern_io::StepExportOptions::default()).map_err(|source| {
        PartError::Io {
            path: format!("{name}.step"),
            source,
        }
    })
}

/// Write [`step_bytes_of`] to `path`, the file stem naming the product.
///
/// # Errors
///
/// As [`step_bytes_of`], and [`PartError::Io`] when the file cannot
/// be written.
pub fn export_step_of(parts: &[Part], path: impl AsRef<std::path::Path>) -> Result<(), PartError> {
    export_step_of_in(parts, path, popkern_io::RepresentationContext::default())
}

/// [`export_step_of`] under an explicit representation context (see
/// [`step_bytes_of_in`]).
///
/// # Errors
///
/// As [`export_step_of`].
pub fn export_step_of_in(
    parts: &[Part],
    path: impl AsRef<std::path::Path>,
    context: popkern_io::RepresentationContext,
) -> Result<(), PartError> {
    let path = path.as_ref();
    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "part".into());
    let bytes = step_bytes_of_in(parts, &name, context)?;
    std::fs::write(path, bytes).map_err(|e| PartError::Io {
        path: path.display().to_string(),
        source: popkern_io::error::IoError::Io(e),
    })
}

/// Round-trip gate over a multi-part export: re-import the STEP at
/// `path` and require every face back, and at least one solid per part
/// (a part with several closed shells exports one solid per shell).
///
/// # Errors
///
/// [`PartError::Io`] on import failure; [`PartError::NotApplicable`]
/// describing the mismatch when faces or solids are missing.
pub fn assert_step_roundtrip_of(
    parts: &[Part],
    path: impl AsRef<std::path::Path>,
) -> Result<(), PartError> {
    let path = path.as_ref();
    let imported = popkern_io::import_step_file(path, &popkern_io::StepImportOptions::default())
        .map_err(|source| PartError::Io {
            path: path.display().to_string(),
            source,
        })?;
    let solids_back: Vec<&BRep> = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .collect();
    let faces_back: usize = solids_back.iter().map(|b| b.faces.len()).sum();
    let faces_out: usize = parts.iter().map(Part::face_count).sum();
    if faces_back != faces_out {
        return Err(PartError::NotApplicable(format!(
            "STEP round-trip returned {faces_back} faces, exported {faces_out}"
        )));
    }
    let shells_out: usize = parts.iter().map(|p| p.brep().shells.len().max(1)).sum();
    if solids_back.len() < parts.len() || solids_back.len() > shells_out {
        return Err(PartError::NotApplicable(format!(
            "STEP round-trip returned {} solids for {} part(s) with {shells_out} shell(s)",
            solids_back.len(),
            parts.len()
        )));
    }
    Ok(())
}
