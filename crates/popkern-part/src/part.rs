//! [`Part`]: a validated solid plus its tolerance and feature log.
//!
//! Every operation is non-destructive at the kernel level (kernel
//! entry points take `&BRep` and return new bodies); the Part API
//! preserves that: a failed operation leaves its inputs untouched
//! and returns the kernel's diagnostics inside [`PartError`].

use std::collections::HashMap;

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::PartError;
use crate::feature::FeatureLog;
use crate::select::{EdgeQuery, edge_infos};

/// A solid part: a validated closed BRep, the tolerance it was built
/// with, and the log of features that produced it.
#[derive(Debug, Clone)]
pub struct Part {
    brep: BRep,
    tolerance: Tolerance,
    features: FeatureLog,
    /// Blend reports accumulated across the part's history (like the
    /// feature log): one entry per fillet/chamfer operation, oldest
    /// first, surviving boolean and derive operations as records.
    blend_history: Vec<crate::dressup::BlendReport>,
}

impl Part {
    /// Wrap a kernel body, validating it first. Used by every
    /// feature and available for interop with hand-built bodies.
    ///
    /// # Errors
    ///
    /// [`PartError::Validation`] when the body fails kernel
    /// validation.
    pub fn from_brep(
        brep: BRep,
        tolerance: Tolerance,
        origin: impl Into<String>,
    ) -> Result<Part, PartError> {
        Part::from_validated(brep, tolerance, "import", origin.into())
    }

    pub(crate) fn from_validated(
        brep: BRep,
        tolerance: Tolerance,
        feature: &'static str,
        summary: String,
    ) -> Result<Part, PartError> {
        // An EMPTY body passes structural validation vacuously; a Part
        // is a servable product and must carry geometry (PK-347: a
        // match-all openings selector removed every face and served a
        // 0-face Part clean).
        if brep.faces.is_empty() {
            return Err(PartError::NotApplicable(format!(
                "{feature} produced an empty body (0 faces): {summary}"
            )));
        }
        brep.validate().map_err(|source| PartError::Validation {
            feature: summary.clone(),
            source,
        })?;
        // Serving-boundary curve-surface gate (check 9, PK-165):
        // kernel pipelines may carry transient curve lies mid-rewrite
        // (exempted from the fast validate() path), but a Part is a
        // finished serving — every edge curve must lie on its
        // adjacent surfaces here (this is the gate that exposed the
        // PK-166 fake annulus and the PK-168 chamfer corner fold).
        brep.validate_curve_surface()
            .map_err(|source| PartError::Validation {
                feature: summary.clone(),
                source,
            })?;
        let mut features = FeatureLog::default();
        features.push(feature, summary);
        Ok(Part {
            brep,
            tolerance,
            features,
            blend_history: Vec::new(),
        })
    }

    /// Wrap a blended body: validates, merges the feature log, and
    /// appends the blend report to the part's blend history.
    pub(crate) fn rebuild_after_blend(
        &self,
        brep: BRep,
        summary: String,
        report: crate::dressup::BlendReport,
    ) -> Result<Part, PartError> {
        let mut part = Part::from_validated(brep, self.tolerance, "fillet", summary)?;
        let mut log = self.features.clone();
        log.extend_from(&part.features);
        part.features = log;
        part.blend_history = self.blend_history.clone();
        part.blend_history.push(report);
        Ok(part)
    }

    /// The report from the most recent blend that produced this
    /// part, if any (warnings, corner patch count, measured
    /// continuity evidence).
    pub fn last_blend_report(&self) -> Option<&crate::dressup::BlendReport> {
        self.blend_history.last()
    }

    /// Every blend report in this part's history, oldest first —
    /// the PK-188 continuity-evidence record. Reports survive later
    /// operations as records; their handles are provenance for the
    /// body they were measured on, not live references.
    pub fn blend_reports(&self) -> &[crate::dressup::BlendReport] {
        &self.blend_history
    }

    /// Prepend `parent`'s feature history (and blend history) to
    /// this part's (used by every derived-part constructor).
    pub(crate) fn merge_history_from(&mut self, parent: &Part) {
        let mut log = parent.features().clone();
        log.extend_from(&self.features);
        self.features = log;
        let mut blends = parent.blend_history.clone();
        blends.append(&mut self.blend_history);
        self.blend_history = blends;
    }

    /// The underlying kernel body.
    pub fn brep(&self) -> &BRep {
        &self.brep
    }

    /// The tolerance every operation on this part uses.
    pub fn tolerance(&self) -> Tolerance {
        self.tolerance
    }

    /// The features that produced this part, oldest first.
    pub fn features(&self) -> &FeatureLog {
        &self.features
    }

    /// Face count (convenience for gates and tests).
    pub fn face_count(&self) -> usize {
        self.brep.faces.len()
    }

    /// Edge count (convenience for gates and tests).
    pub fn edge_count(&self) -> usize {
        self.brep.edges.len()
    }

    // ───────────────────── combination ─────────────────────

    /// Boolean union with another part.
    ///
    /// Both parts survive unchanged on failure. Mixed tolerances
    /// take the finer of the two (recorded in the feature log).
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn union(&self, other: &Part) -> Result<Part, PartError> {
        self.combine(other, popkern_bool::BoolOperation::Union, "union")
    }

    /// Boolean subtraction (`self − other`).
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn subtract(&self, other: &Part) -> Result<Part, PartError> {
        self.combine(other, popkern_bool::BoolOperation::Subtraction, "subtract")
    }

    /// Boolean intersection.
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn intersect(&self, other: &Part) -> Result<Part, PartError> {
        self.combine(
            other,
            popkern_bool::BoolOperation::Intersection,
            "intersect",
        )
    }

    fn combine(
        &self,
        other: &Part,
        op: popkern_bool::BoolOperation,
        name: &'static str,
    ) -> Result<Part, PartError> {
        let tolerance = finer(self.tolerance, other.tolerance);
        let result =
            popkern_bool::boolean(op, self.brep.clone(), other.brep.clone(), tolerance, None)
                .map_err(|source| PartError::Boolean {
                    operation: name,
                    source,
                })?;
        let mut part = Part::from_validated(
            result,
            tolerance,
            name,
            format!(
                "{name}: {} faces ⊕ {} faces",
                self.face_count(),
                other.face_count()
            ),
        )?;
        let mut log = self.features.clone();
        log.extend_from(&other.features);
        log.extend_from(&part.features);
        part.features = log;
        part.blend_history = self
            .blend_history
            .iter()
            .chain(other.blend_history.iter())
            .cloned()
            .collect();
        Ok(part)
    }

    // ───────────────────── dress-up ─────────────────────

    /// Fillet every edge matched by the selector with a constant
    /// radius. All matched edges blend in ONE kernel operation so
    /// chains and their corner blends resolve together.
    ///
    /// ```no_run
    /// # use popkern_part::prelude::*;
    /// # fn f(part: Part) -> Result<Part, PartError> {
    /// part.fillet(0.08, |e| e.circular().or(|e| e.straight().at_height(0.4)))
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// [`PartError::EmptySelection`] when the query matches nothing
    /// (unless `.allow_empty()`); [`PartError::Blend`] with the
    /// kernel diagnostics on blend failure.
    pub fn fillet(
        &self,
        radius: f64,
        select: impl FnOnce(EdgeQuery) -> EdgeQuery,
    ) -> Result<Part, PartError> {
        let query = select(EdgeQuery::new());
        let infos = edge_infos(&self.brep);
        let mut edges = query.resolve(&infos);
        let selection = query.describe();
        // Seam edges (same face on both sides) can never blend —
        // drop them from the selection rather than failing the
        // whole operation, and record the skip.
        let seams: std::collections::HashSet<_> = infos
            .iter()
            .filter(|i| i.is_seam)
            .map(|i| i.handle)
            .collect();
        let skipped_seams = edges.iter().filter(|e| seams.contains(e)).count();
        edges.retain(|e| !seams.contains(e));
        if edges.is_empty() {
            if query.empty_permitted() {
                let mut part = self.clone();
                part.features
                    .push("fillet", format!("no-op (empty selection): {selection}"));
                return Ok(part);
            }
            return Err(PartError::EmptySelection {
                selection,
                entity: "edges",
                available: infos.len(),
            });
        }
        let specs: Vec<FilletSpec> = edges
            .iter()
            .map(|&edge| FilletSpec {
                edge,
                radius: FilletRadius::Constant(radius),
            })
            .collect();
        let options = FilletOptions::default();
        let pre_faces: std::collections::BTreeSet<_> = self.brep.faces.keys().collect();
        let result =
            fillet_edges(&self.brep, &specs, &options, self.tolerance, None).map_err(|source| {
                PartError::Blend {
                    selection: selection.clone(),
                    source,
                }
            })?;
        let seam_note = if skipped_seams > 0 {
            format!(" (skipped {skipped_seams} seam edges)")
        } else {
            String::new()
        };
        let evidence =
            crate::dressup::measure_blend_evidence(&result.brep, &pre_faces, options.continuity);
        let report = crate::dressup::BlendReport {
            edges_blended: edges.len(),
            warnings: result
                .edge_warnings
                .iter()
                .map(|(eh, w)| format!("{eh:?}: {w:?}"))
                .collect(),
            corner_patches: result.corner_patches.len(),
            requested_continuity: options.continuity,
            continuity: Some(evidence),
        };
        self.rebuild_after_blend(
            result.brep,
            format!(
                "radius {radius} on {} edges: {selection}{seam_note}",
                edges.len()
            ),
            report,
        )
    }

    // ───────────────────── interrogation & gates ─────────────────────

    /// Tessellate the part.
    ///
    /// # Errors
    ///
    /// [`PartError::Mesh`] when tessellation fails.
    pub fn mesh(
        &self,
        options: &popkern_mesh::TessellationOptions,
    ) -> Result<popkern_mesh::TriMesh, PartError> {
        popkern_mesh::tessellate(&self.brep, options).map_err(|source| PartError::Mesh { source })
    }

    /// Count of open (single-use) edges in the welded default-quality
    /// mesh — 0 means watertight.
    ///
    /// # Errors
    ///
    /// [`PartError::Mesh`] when tessellation fails.
    pub fn open_edge_count(&self) -> Result<usize, PartError> {
        let mesh = self.mesh(&popkern_mesh::TessellationOptions::default())?;
        Ok(welded_open_edges(&mesh))
    }

    /// Watertightness gate: error unless the tessellated part has
    /// zero open mesh edges. Returns `self` for chaining.
    ///
    /// # Errors
    ///
    /// [`PartError::NotWatertight`] with the open-edge count.
    pub fn assert_watertight(self) -> Result<Part, PartError> {
        let open_edges = self.open_edge_count()?;
        if open_edges != 0 {
            return Err(PartError::NotWatertight { open_edges });
        }
        Ok(self)
    }

    /// Export the tessellated part as Wavefront OBJ.
    ///
    /// # Errors
    ///
    /// [`PartError::Mesh`] on tessellation failure; [`PartError::Io`]
    /// on write failure.
    pub fn export_obj(&self, path: impl AsRef<std::path::Path>) -> Result<(), PartError> {
        let path = path.as_ref();
        let mesh = self.mesh(&popkern_mesh::TessellationOptions::default())?;
        let file = std::fs::File::create(path).map_err(|e| {
            PartError::NotApplicable(format!("cannot create '{}': {e}", path.display()))
        })?;
        let mut writer = std::io::BufWriter::new(file);
        popkern_mesh::export_obj(&mesh, &mut writer).map_err(|source| PartError::Mesh { source })
    }

    /// Export the tessellated part as STL (binary by default).
    ///
    /// # Errors
    ///
    /// [`PartError::Mesh`] on tessellation or write failure.
    pub fn export_stl(&self, path: impl AsRef<std::path::Path>) -> Result<(), PartError> {
        let path = path.as_ref();
        let mesh = self.mesh(&popkern_mesh::TessellationOptions::default())?;
        popkern_mesh::export_stl(&mesh, path, true).map_err(|source| PartError::Mesh { source })
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
        let path = path.as_ref();
        let imported =
            popkern_io::import_step_file(path, &popkern_io::StepImportOptions::default()).map_err(
                |source| PartError::Io {
                    path: path.display().to_string(),
                    source,
                },
            )?;
        let faces_back: usize = imported
            .document
            .products
            .iter()
            .flat_map(|p| p.shapes.iter())
            .flat_map(|s| s.solids.iter())
            .map(|b| b.faces.len())
            .sum();
        if faces_back != self.face_count() {
            return Err(PartError::NotApplicable(format!(
                "STEP round-trip returned {faces_back} faces, exported {}",
                self.face_count()
            )));
        }
        Ok(self)
    }
}

/// The finer (smaller linear value) of two tolerances.
fn finer(a: Tolerance, b: Tolerance) -> Tolerance {
    if a.linear() <= b.linear() { a } else { b }
}

/// Welded open-edge census (the watertightness measure every example
/// gate uses): weld mesh vertices at 1e-7, count edges used by
/// exactly one triangle.
fn welded_open_edges(mesh: &popkern_mesh::TriMesh) -> usize {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: popkern_core::geometry::Point3| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}
