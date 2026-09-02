//! STEP export of a [`Part`]: the path form and the in-memory bytes
//! form (hosts without a filesystem — the browser playground — hand
//! the bytes to a download).

use crate::error::PartError;
use crate::part::Part;

impl Part {
    /// Export the part as a STEP solid.
    ///
    /// # Errors
    ///
    /// [`PartError::Io`] when export fails.
    pub fn export_step(&self, path: impl AsRef<std::path::Path>) -> Result<(), PartError> {
        let path = path.as_ref();
        let name = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "part".into());
        let bytes = self.step_bytes(&name)?;
        std::fs::write(path, bytes).map_err(|e| PartError::Io {
            path: path.display().to_string(),
            source: popkern_io::error::IoError::Io(e),
        })
    }

    /// The STEP solid as Part 21 bytes, with `name` as the product
    /// name (the file stem when written to disk). This is the in-memory
    /// form [`Part::export_step`] writes: a host without a filesystem
    /// (the browser playground) hands these bytes to a download.
    ///
    /// # Errors
    ///
    /// [`PartError::Io`] when the STEP writer declines the body.
    pub fn step_bytes(&self, name: &str) -> Result<Vec<u8>, PartError> {
        let doc = popkern_io::StepDocument {
            products: vec![popkern_io::StepProductNode {
                name: name.to_string(),
                id: name.to_uppercase(),
                children: vec![],
                shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                    kind: popkern_io::ShapeRepKind::AdvancedBrep,
                    solids: vec![self.brep().clone()],
                    context: popkern_io::RepresentationContext::default(),
                    non_brep_items: vec![],
                })],
                metadata: popkern_io::ProductMetadata::default(),
            }],
            roots: vec![0],
            units: popkern_io::UnitContext::default(),
            opaque_bag: popkern_io::OpaqueBag::default(),
        };
        popkern_io::export_step(&doc, &popkern_io::StepExportOptions::default()).map_err(|source| {
            PartError::Io {
                path: format!("{name}.step"),
                source,
            }
        })
    }
}
