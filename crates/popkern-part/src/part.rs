//! [`Part`]: a validated solid plus its tolerance and feature log.
//!
//! Every operation is non-destructive at the kernel level (kernel
//! entry points take `&BRep` and return new bodies); the Part API
//! preserves that: a failed operation leaves its inputs untouched
//! and returns the kernel's diagnostics inside [`PartError`].

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

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
    last_blend: Option<crate::dressup::BlendReport>,
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
        brep.validate().map_err(|source| PartError::Validation {
            feature: summary.clone(),
            source,
        })?;
        let mut features = FeatureLog::default();
        features.push(feature, summary);
        Ok(Part {
            brep,
            tolerance,
            features,
            last_blend: None,
        })
    }

    /// Wrap a blended body: validates, merges the feature log, and
    /// stores the blend report.
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
        part.last_blend = Some(report);
        Ok(part)
    }

    /// The report from the most recent blend that produced this
    /// part, if any (warnings, corner patch count).
    pub fn last_blend_report(&self) -> Option<&crate::dressup::BlendReport> {
        self.last_blend.as_ref()
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

    // ───────────────────── polyhedra ─────────────────────

    /// Build a part as the convex hull of `points` (coplanar hull
    /// triangles merged into polygons). Intended for canonical
    /// solids — the hull is O(n⁴) plane enumeration.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] for degenerate input;
    /// validation errors if wiring fails.
    pub fn polyhedron(
        points: &[popkern_core::geometry::Point3],
        tolerance: Tolerance,
    ) -> Result<Part, PartError> {
        let faces = crate::poly::convex_hull_faces(points);
        if faces.len() < 4 {
            return Err(PartError::NotApplicable(format!(
                "convex hull of {} points produced only {} faces",
                points.len(),
                faces.len()
            )));
        }
        let specs: Vec<crate::poly::PolyFace> = faces
            .into_iter()
            .map(crate::poly::PolyFace::simple)
            .collect();
        let brep = crate::poly::build_polyhedron_brep(points, &specs, tolerance)?;
        Part::from_validated(
            brep,
            tolerance,
            "polyhedron",
            format!("convex hull of {} points", points.len()),
        )
    }

    /// Uniform vertex truncation on a polyhedral part: every vertex
    /// becomes a flat corner facet cut `fraction` (0 < f < 0.5)
    /// along each incident edge.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the part is not polyhedral
    /// (non-planar faces or curved edges) or the fraction is out of
    /// range.
    pub fn truncate_vertices(&self, fraction: f64) -> Result<Part, PartError> {
        if !(fraction > 0.0 && fraction < 0.5) {
            return Err(PartError::NotApplicable(format!(
                "truncation fraction must be in (0, 0.5), got {fraction}"
            )));
        }
        let (points, faces) = crate::poly::extract_polyhedron(&self.brep)?;
        if faces.iter().any(|f| !f.inners.is_empty()) {
            return Err(PartError::NotApplicable(
                "truncate_vertices requires faces without inner loops (truncate before engraving)"
                    .into(),
            ));
        }
        let outers: Vec<Vec<usize>> = faces.iter().map(|f| f.outer.clone()).collect();
        let (tp, tf) = crate::poly::truncate_vertices(&points, &outers, fraction);
        let specs: Vec<crate::poly::PolyFace> =
            tf.into_iter().map(crate::poly::PolyFace::simple).collect();
        let brep = crate::poly::build_polyhedron_brep(&tp, &specs, self.tolerance)?;
        let mut part = Part::from_validated(
            brep,
            self.tolerance,
            "truncate_vertices",
            format!("fraction {fraction}"),
        )?;
        let mut log = self.features.clone();
        log.extend_from(&part.features);
        part.features = log;
        Ok(part)
    }

    /// Edit the part as a polyhedron representation: the closure
    /// receives `(points, faces)` (see [`crate::poly`]) and may
    /// engrave pockets or otherwise extend them; the result is
    /// rebuilt and validated. The escape hatch for per-face work
    /// like dice numbering.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the part is not polyhedral;
    /// whatever the closure returns; validation errors on rebuild.
    pub fn edit_polyhedron(
        &self,
        summary: impl Into<String>,
        edit: impl FnOnce(
            &mut Vec<popkern_core::geometry::Point3>,
            &mut Vec<crate::poly::PolyFace>,
        ) -> Result<(), PartError>,
    ) -> Result<Part, PartError> {
        let (mut points, mut faces) = crate::poly::extract_polyhedron(&self.brep)?;
        edit(&mut points, &mut faces)?;
        let brep = crate::poly::build_polyhedron_brep(&points, &faces, self.tolerance)?;
        let mut part =
            Part::from_validated(brep, self.tolerance, "edit_polyhedron", summary.into())?;
        let mut log = self.features.clone();
        log.extend_from(&part.features);
        part.features = log;
        Ok(part)
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
        let result = fillet_edges(
            &self.brep,
            &specs,
            &FilletOptions::default(),
            self.tolerance,
            None,
        )
        .map_err(|source| PartError::Blend {
            selection: selection.clone(),
            source,
        })?;
        let seam_note = if skipped_seams > 0 {
            format!(" (skipped {skipped_seams} seam edges)")
        } else {
            String::new()
        };
        let report = crate::dressup::BlendReport {
            edges_blended: edges.len(),
            warnings: result
                .edge_warnings
                .iter()
                .map(|(eh, w)| format!("{eh:?}: {w:?}"))
                .collect(),
            corner_patches: result.corner_patches.len(),
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

    // ───────────────────── export ─────────────────────

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
        let doc = popkern_io::StepDocument {
            products: vec![popkern_io::StepProductNode {
                name: name.clone(),
                id: name.to_uppercase(),
                children: vec![],
                shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                    kind: popkern_io::ShapeRepKind::AdvancedBrep,
                    solids: vec![self.brep.clone()],
                    context: popkern_io::RepresentationContext::default(),
                    non_brep_items: vec![],
                })],
                metadata: popkern_io::ProductMetadata::default(),
            }],
            roots: vec![0],
            units: popkern_io::UnitContext::default(),
            opaque_bag: popkern_io::OpaqueBag::default(),
        };
        popkern_io::export_step_file(&doc, path, &popkern_io::StepExportOptions::default()).map_err(
            |source| PartError::Io {
                path: path.display().to_string(),
                source,
            },
        )
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
