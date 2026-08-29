//! Shaping operations beyond blends: chamfer, shell, draft, and
//! rigid transforms.

use popkern_blend::{ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};
use popkern_core::geometry::transform::Transform3;

use crate::error::PartError;
use crate::part::Part;
use crate::select::{EdgeQuery, FaceQuery, edge_infos, face_infos};

impl Part {
    /// Equal-distance (45°) chamfer on every selected edge, in one
    /// kernel operation.
    ///
    /// # Errors
    ///
    /// As [`Part::fillet`]: empty selection unless permitted; kernel
    /// blend diagnostics on failure.
    pub fn chamfer(
        &self,
        distance: f64,
        select: impl FnOnce(EdgeQuery) -> EdgeQuery,
    ) -> Result<Part, PartError> {
        let query = select(EdgeQuery::new());
        let infos = edge_infos(self.brep());
        let mut edges = query.resolve(&infos);
        let selection = query.describe();
        let seams: std::collections::HashSet<_> = infos
            .iter()
            .filter(|i| i.is_seam)
            .map(|i| i.handle)
            .collect();
        edges.retain(|e| !seams.contains(e));
        if edges.is_empty() {
            if query.empty_permitted() {
                return Ok(self.clone());
            }
            return Err(PartError::EmptySelection {
                selection,
                entity: "edges",
                available: infos.len(),
            });
        }
        let specs: Vec<ChamferSpec> = edges
            .iter()
            .map(|&edge| ChamferSpec {
                edge,
                kind: ChamferKind::Equal(distance),
            })
            .collect();
        // Recorded-survivors semantics: the kernel is transactional
        // (ANY failed leg fails the whole op — PK-140 closed the
        // silent partial-apply hole). Product dress-up keeps the
        // best-effort UX: on PartialFailure, retry once WITHOUT the
        // failing edges and surface every skipped leg as a report
        // warning — recorded, never silent (measured: a
        // slot-interrupted chamfer shipped for months with 2 of 4
        // legs missing and no diagnostic).
        let mut skipped: Vec<String> = Vec::new();
        let attempt = chamfer_edges(
            self.brep(),
            &specs,
            &ChamferOptions::default(),
            self.tolerance(),
            None,
        );
        let result = match attempt {
            Ok(r) => r,
            Err(popkern_blend::BlendError::PartialFailure { ref failures, .. })
                if failures.len() < specs.len() =>
            {
                let failed: std::collections::HashSet<_> =
                    failures.iter().map(|(eh, _)| *eh).collect();
                for (eh, failure) in failures {
                    let msg = format!("{failure:?}");
                    skipped.push(format!(
                        "chamfer leg skipped {eh:?}: {}",
                        &msg[..msg.len().min(160)]
                    ));
                }
                let surviving: Vec<ChamferSpec> = specs
                    .iter()
                    .filter(|s| !failed.contains(&s.edge))
                    .cloned()
                    .collect();
                chamfer_edges(
                    self.brep(),
                    &surviving,
                    &ChamferOptions::default(),
                    self.tolerance(),
                    None,
                )
                .map_err(|source| PartError::Blend {
                    selection: selection.clone(),
                    source,
                })?
            }
            Err(source) => {
                return Err(PartError::Blend {
                    selection: selection.clone(),
                    source,
                });
            }
        };
        let report = crate::dressup::BlendReport {
            edges_blended: edges.len() - skipped.len(),
            warnings: skipped
                .into_iter()
                .chain(
                    result
                        .edge_warnings
                        .iter()
                        .map(|(eh, w)| format!("{eh:?}: {w:?}")),
                )
                .collect(),
            corner_patches: result.corner_patches.len(),
            // A chamfer's boundary is a designed crease: nothing
            // smooth was requested, so there is no continuity claim
            // to measure evidence for.
            requested_continuity: crate::dressup::Continuity::G0,
            continuity: None,
        };
        self.rebuild_after_blend(
            result.brep,
            format!("chamfer {distance} on {} edges: {selection}", edges.len()),
            report,
        )
    }

    /// Hollow the part to a wall `thickness`, removing the selected
    /// faces as openings (an empty-permitted query with no matches
    /// produces a closed hollow solid). Material is removed inward.
    ///
    /// # Errors
    ///
    /// [`PartError::EmptySelection`] when the opening query matches
    /// nothing and did not opt into `.allow_empty()`;
    /// [`PartError::NotApplicable`] carrying the kernel diagnostics
    /// on shell failure.
    pub fn shell(
        &self,
        thickness: f64,
        openings: impl FnOnce(FaceQuery) -> FaceQuery,
    ) -> Result<Part, PartError> {
        self.shell_with_options(
            thickness,
            popkern_shell::ShellDirection::Inward,
            &popkern_shell::ShellOptions::default(),
            openings,
        )
    }

    /// Draft the selected faces by `angle_radians` for mold
    /// withdrawal along `pull` (neutral plane through the part's
    /// lowest point, perpendicular to the pull direction).
    ///
    /// # Errors
    ///
    /// [`PartError::EmptySelection`] on an empty face selection;
    /// [`PartError::NotApplicable`] carrying kernel diagnostics.
    pub fn draft(
        &self,
        angle_radians: f64,
        pull: crate::select::Dir,
        select: impl FnOnce(FaceQuery) -> FaceQuery,
    ) -> Result<Part, PartError> {
        let query = select(FaceQuery::new());
        let infos = face_infos(self.brep());
        let faces = query.resolve(&infos);
        let selection = query.describe();
        if faces.is_empty() {
            if query.empty_permitted() {
                return Ok(self.clone());
            }
            return Err(PartError::EmptySelection {
                selection,
                entity: "faces",
                available: infos.len(),
            });
        }
        let pull_direction = pull.vector();
        // Neutral plane: perpendicular to pull, through the lowest
        // vertex along the pull axis (the parting side stays put).
        let zero = popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
        let min_t = self
            .brep()
            .vertices
            .iter()
            .map(|(_, v)| (v.position - zero).dot(&pull_direction))
            .fold(f64::INFINITY, f64::min);
        let origin = zero + pull_direction * min_t;
        let specs: Vec<popkern_draft::DraftSpec> = faces
            .iter()
            .map(|&face| popkern_draft::DraftSpec {
                face,
                angle: angle_radians,
                pull_direction,
                neutral: popkern_draft::NeutralReference::Plane {
                    origin,
                    normal: pull_direction,
                },
                angle_variation: popkern_draft::AngleVariation::Constant,
            })
            .collect();
        let result = popkern_draft::draft_faces(
            self.brep(),
            &specs,
            &popkern_draft::DraftOptions::default(),
            self.tolerance(),
            None,
        )
        .map_err(|e| PartError::NotApplicable(format!("draft failed: {e}")))?;
        let mut part = Part::from_validated(
            result.brep,
            self.tolerance(),
            "draft",
            format!(
                "{angle_radians} rad along {pull:?} on {} faces: {selection}",
                faces.len()
            ),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Blend two faces that share no edge with a constant-radius
    /// rolling ball, trimming both walls back to the contact curves
    /// and seaming the strip in (the PK-142 face–face blend,
    /// `FaceBlendAttach::TrimBoth`).
    ///
    /// Each selector must resolve to exactly one face — a face–face
    /// blend is a pairwise operation, and an ambiguous selection is
    /// an error rather than a guess.
    ///
    /// # Errors
    ///
    /// [`PartError::EmptySelection`] when a selector matches nothing;
    /// [`PartError::NotApplicable`] when a selector is ambiguous or
    /// the kernel declines the blend (structured diagnostics
    /// preserved in the message).
    pub fn face_blend(
        &self,
        radius: f64,
        select_a: impl FnOnce(FaceQuery) -> FaceQuery,
        select_b: impl FnOnce(FaceQuery) -> FaceQuery,
    ) -> Result<Part, PartError> {
        let infos = face_infos(self.brep());
        let resolve_one = |label: &'static str,
                           query: FaceQuery|
         -> Result<popkern_core::topology::FaceHandle, PartError> {
            let faces = query.resolve(&infos);
            let selection = query.describe();
            match faces.as_slice() {
                [] => Err(PartError::EmptySelection {
                    selection,
                    entity: "faces",
                    available: infos.len(),
                }),
                [face] => Ok(*face),
                many => Err(PartError::NotApplicable(format!(
                    "face_blend {label} selector must resolve to exactly one face; \
                     `{selection}` matched {}",
                    many.len()
                ))),
            }
        };
        let face_a = resolve_one("first", select_a(FaceQuery::new()))?;
        let face_b = resolve_one("second", select_b(FaceQuery::new()))?;
        let spec = popkern_blend::FaceBlendSpec {
            face_a,
            face_b,
            radius,
            radius_end: None,
            attach: popkern_blend::FaceBlendAttach::TrimBoth,
        };
        let result = popkern_blend::face_blend(self.brep(), &spec, self.tolerance())
            .map_err(|e| PartError::NotApplicable(format!("face_blend failed: {e}")))?;
        let attached = result.attached.ok_or_else(|| {
            PartError::NotApplicable(
                "face_blend attach declined: kernel returned the strip sheet only".into(),
            )
        })?;
        let mut part = Part::from_validated(
            attached,
            self.tolerance(),
            "face_blend",
            format!("radius {radius} between {face_a:?} and {face_b:?}"),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Rigid translation.
    ///
    /// # Errors
    ///
    /// Kernel transform or validation errors.
    pub fn translate(&self, dx: f64, dy: f64, dz: f64) -> Result<Part, PartError> {
        self.transformed(
            &Transform3::translation(dx, dy, dz),
            format!("translate ({dx}, {dy}, {dz})"),
        )
    }

    /// Rigid rotation about the global X axis (radians).
    ///
    /// # Errors
    ///
    /// Kernel transform or validation errors.
    pub fn rotate_x(&self, angle_radians: f64) -> Result<Part, PartError> {
        self.transformed(
            &Transform3::rotation_x(angle_radians),
            format!("rotate_x {angle_radians}"),
        )
    }

    /// Rigid rotation about the global Y axis (radians).
    ///
    /// # Errors
    ///
    /// Kernel transform or validation errors.
    pub fn rotate_y(&self, angle_radians: f64) -> Result<Part, PartError> {
        self.transformed(
            &Transform3::rotation_y(angle_radians),
            format!("rotate_y {angle_radians}"),
        )
    }

    /// Rigid rotation about the global Z axis (radians).
    ///
    /// # Errors
    ///
    /// Kernel transform or validation errors.
    pub fn rotate_z(&self, angle_radians: f64) -> Result<Part, PartError> {
        self.transformed(
            &Transform3::rotation_z(angle_radians),
            format!("rotate_z {angle_radians}"),
        )
    }

    /// Apply an arbitrary rigid transform.
    ///
    /// # Errors
    ///
    /// Kernel transform or validation errors.
    pub fn transformed(&self, t: &Transform3, summary: String) -> Result<Part, PartError> {
        let brep = self.brep().transform(t).map_err(|source| PartError::Core {
            feature: summary.clone(),
            source,
        })?;
        let mut part = Part::from_validated(brep, self.tolerance(), "transform", summary)?;
        part.merge_history_from(self);
        Ok(part)
    }
}
