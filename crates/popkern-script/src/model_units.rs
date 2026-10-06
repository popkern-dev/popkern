//! Document length units and the STEP export that declares them
//! (split from `model.rs` under the 500-line file budget; PK-462).

use popkern_history::unit::Unit;

use crate::error::ScriptError;
use crate::model::Model;

impl Model {
    /// Set the length unit new literals are recorded in.
    ///
    /// The build tolerance is a length like any other literal, so it
    /// follows the unit (PK-462): a document's tolerance means the
    /// same fraction of its geometry whatever unit it is authored in.
    /// A mm model with the default tolerance builds at 1e-6 mm (1e-9
    /// m) — at 1e-6 m an M6 thread ran with a tolerance 1000x coarser
    /// relative to its features than the same part authored in mm
    /// magnitudes, and every weld gate landed at feature scale.
    pub fn unit(mut self, unit: Unit) -> Model {
        let factor = |u: Unit| match u {
            Unit::Meters | Unit::Millimeters | Unit::Centimeters | Unit::Inches => {
                Some(u.to_si(1.0))
            }
            _ => None,
        };
        if let (Some(from), Some(to)) = (factor(self.unit), factor(unit)) {
            let linear = self.tolerance.linear() * to / from;
            self.tolerance = self.tolerance.with_linear(linear);
            self.tree.set_build_tolerance(Some(self.tolerance));
        }
        self.unit = unit;
        self
    }

    /// The length unit literals are recorded in (default metres).
    pub fn length_unit(&self) -> Unit {
        self.unit
    }

    /// Export the whole model as STEP: every alive body as its own
    /// solid in one product (a gearbox's six gears, not the active
    /// one). A single-body model writes exactly what
    /// [`Part::export_step`] writes for that body.
    pub fn to_step(&mut self, path: impl AsRef<std::path::Path>) -> Result<(), ScriptError> {
        let parts = self.to_parts_all()?;
        let context = self.step_context();
        Ok(popkern_part::export_step_of_in(&parts, path, context)?)
    }

    /// The STEP representation context this document's lengths are
    /// declared in (PK-462). A mm, cm or inch document stores true SI
    /// lengths and declares metres, so a consumer reads an M6 as 6 mm;
    /// declaring millimetres made it 0.006 mm. A metre (default-unit)
    /// document keeps the millimetre declaration of its
    /// millimetre-magnitude convention. The uncertainty is the
    /// document's own tolerance.
    fn step_context(&self) -> popkern_io::RepresentationContext {
        let mut context = popkern_io::RepresentationContext::default();
        if matches!(
            self.unit,
            Unit::Millimeters | Unit::Centimeters | Unit::Inches
        ) {
            context.length_unit = popkern_io::LengthUnit::Metre;
            context.uncertainty = self.tolerance.linear();
        }
        context
    }

    /// The whole model as STEP Part 21 bytes with `name` as the
    /// product name — what [`Model::to_step`] writes for a path whose
    /// file stem is `name`. Hosts without a filesystem (popkern
    /// live) hand these to a download.
    pub fn to_step_bytes(&mut self, name: &str) -> Result<Vec<u8>, ScriptError> {
        let parts = self.to_parts_all()?;
        let context = self.step_context();
        Ok(popkern_part::step_bytes_of_in(&parts, name, context)?)
    }
}
