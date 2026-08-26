//! Offset-family shaping methods (PK-328): options pass-through for
//! shell/thicken/offset plus the pre-flight `max_safe_offset` query.
//! Split from `shaping.rs` (500-line ratchet).

use crate::error::PartError;
use crate::part::Part;
use crate::select::{FaceQuery, face_infos};

impl Part {
    /// [`Part::shell`] with an explicit direction and the full
    /// [`popkern_shell::ShellOptions`] knob set (join style, trim
    /// policy, partial-face handling — PK-328 options pass-through).
    ///
    /// # Errors
    ///
    /// As [`Part::shell`].
    pub fn shell_with_options(
        &self,
        thickness: f64,
        direction: popkern_shell::ShellDirection,
        options: &popkern_shell::ShellOptions,
        openings: impl FnOnce(FaceQuery) -> FaceQuery,
    ) -> Result<Part, PartError> {
        let query = openings(FaceQuery::new());
        let infos = face_infos(self.brep());
        let faces = query.resolve(&infos);
        let selection = query.describe();
        if faces.is_empty() && !query.empty_permitted() {
            return Err(PartError::EmptySelection {
                selection,
                entity: "faces",
                available: infos.len(),
            });
        }
        let result = popkern_shell::shell(
            self.brep(),
            thickness,
            direction,
            &faces,
            &[],
            options,
            self.tolerance(),
            None,
        )
        .map_err(|e| PartError::NotApplicable(format!("shell failed: {e}")))?;
        let mut part = Part::from_validated(
            result.brep,
            self.tolerance(),
            "shell",
            format!(
                "thickness {thickness}, {direction:?}, {:?}/{:?}, {} openings: {selection}",
                options.join_style,
                options.trim_policy,
                faces.len()
            ),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Thicken the body symmetrically by ±`thickness`/2 (PK-328).
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] carrying the kernel diagnostics.
    pub fn thicken(&self, thickness: f64) -> Result<Part, PartError> {
        self.thicken_with_options(thickness, &popkern_shell::ShellOptions::default())
    }

    /// [`Part::thicken`] with the full options knob set.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] carrying the kernel diagnostics.
    pub fn thicken_with_options(
        &self,
        thickness: f64,
        options: &popkern_shell::ShellOptions,
    ) -> Result<Part, PartError> {
        let result =
            popkern_shell::thicken::thicken(self.brep(), thickness, options, self.tolerance())
                .map_err(|e| PartError::NotApplicable(format!("thicken failed: {e}")))?;
        let mut part = Part::from_validated(
            result.brep,
            self.tolerance(),
            "thicken",
            format!(
                "thickness {thickness}, {:?}/{:?}",
                options.join_style, options.trim_policy
            ),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Offset every face of the solid: positive `distance` grows the
    /// body, negative shrinks it (PK-328).
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] carrying the kernel diagnostics.
    pub fn offset(&self, distance: f64) -> Result<Part, PartError> {
        self.offset_with_options(distance, &popkern_shell::ShellOptions::default())
    }

    /// [`Part::offset`] with the full options knob set.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] carrying the kernel diagnostics.
    pub fn offset_with_options(
        &self,
        distance: f64,
        options: &popkern_shell::ShellOptions,
    ) -> Result<Part, PartError> {
        let result =
            popkern_shell::offset::offset_solid(self.brep(), distance, options, self.tolerance())
                .map_err(|e| PartError::NotApplicable(format!("offset failed: {e}")))?;
        let mut part = Part::from_validated(
            result.brep,
            self.tolerance(),
            "offset",
            format!(
                "distance {distance}, {:?}/{:?}",
                options.join_style, options.trim_policy
            ),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// A cheap approximate preview of the offset skin (exact
    /// offset-map samples, approximate join strips, no trim) for
    /// interactive display while a slider drags (PK-328). Never a
    /// model — display only.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when classification fails.
    pub fn shell_preview(
        &self,
        thickness: f64,
        direction: popkern_shell::ShellDirection,
    ) -> Result<popkern_shell::PreviewMesh, PartError> {
        popkern_shell::shell_preview(self.brep(), thickness, direction, self.tolerance())
            .map_err(|e| PartError::NotApplicable(format!("shell preview failed: {e}")))
    }

    /// The largest safe offset magnitude for `direction` — the slider
    /// upper bound a UI reads before running shell/offset/thicken
    /// (PK-328, interactive).
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the material side cannot be
    /// resolved.
    pub fn max_safe_offset(
        &self,
        direction: popkern_shell::ShellDirection,
    ) -> Result<popkern_shell::SafeOffsetReport, PartError> {
        popkern_shell::max_safe_offset(self.brep(), direction, self.tolerance())
            .map_err(|e| PartError::NotApplicable(format!("max_safe_offset failed: {e}")))
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    fn box_part() -> Part {
        Profile::on(Plane::XY)
            .rect(1.0, 1.0)
            .extrude(0.5)
            .expect("box extrudes")
    }

    #[test]
    fn shell_with_options_outward_grows_join_faces() {
        let part = box_part();
        let rounded = part
            .shell_with_options(
                0.1,
                popkern_shell::ShellDirection::Outward,
                &popkern_shell::ShellOptions::default(),
                // No openings: an explicitly-empty selection (a bare
                // query matches every face).
                |f| f.filter(|_| false).allow_empty(),
            )
            .expect("outward shell serves");
        // Rolling-ball joins: 6+6 skins + 12 strips + 8 corners.
        assert!(
            rounded.brep().faces.len() > 12,
            "arc joins must appear, got {}",
            rounded.brep().faces.len()
        );
        assert!(rounded.brep().validate_all().is_empty());
    }

    #[test]
    fn thicken_and_offset_serve_and_log() {
        let part = box_part();
        let grown = part.offset(0.1).expect("offset serves");
        assert!(
            grown
                .features()
                .entries()
                .iter()
                .any(|f| f.name == "offset")
        );
        let thick = part.thicken(0.1).expect("thicken serves");
        assert!(
            thick
                .features()
                .entries()
                .iter()
                .any(|f| f.name == "thicken")
        );
    }

    #[test]
    fn max_safe_offset_reports_slider_bounds() {
        let part = box_part();
        let inward = part
            .max_safe_offset(popkern_shell::ShellDirection::Inward)
            .expect("inward bound");
        assert!(
            (inward.max_distance - 0.25 * 0.99).abs() < 1e-6,
            "box half-thickness binds: {inward}"
        );
        let outward = part
            .max_safe_offset(popkern_shell::ShellDirection::Outward)
            .expect("outward bound");
        assert!(outward.max_distance.is_infinite(), "{outward}");
    }
}
