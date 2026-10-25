//! Conversions shared by the structural and thermal setup files: units,
//! required keys, face tables, tolerances and options.

use std::collections::BTreeMap;

use popkern_analysis::{AnalysisOptions, Expect, FaceDef, FaceRef, Tolerance};

use super::convert::{face_def, selector_file};
use super::{ExpectFile, FaceEntry, OptionsFile, SetupFileError, ToleranceFile};

/// `analysis` must be absent or equal `expected`.
pub(crate) fn check_analysis(
    found: Option<&str>,
    expected: &'static str,
) -> Result<(), SetupFileError> {
    match found {
        None => Ok(()),
        Some(a) if a == expected => Ok(()),
        Some(a) => Err(SetupFileError::WrongAnalysis {
            expected,
            found: a.to_string(),
        }),
    }
}

/// `units` must be present and `"mm"`.
pub(crate) fn check_units(units: Option<&str>) -> Result<(), SetupFileError> {
    match units {
        None => Err(SetupFileError::MissingField {
            field: "units".to_string(),
        }),
        Some("mm") => Ok(()),
        Some(other) => Err(SetupFileError::UnsupportedUnits(other.to_string())),
    }
}

/// A required key's value.
pub(crate) fn required<T>(field: &str, value: Option<T>) -> Result<T, SetupFileError> {
    value.ok_or_else(|| SetupFileError::MissingField {
        field: field.to_string(),
    })
}

/// The `[faces]` table as face definitions, in name order.
pub(crate) fn face_defs(
    faces: &BTreeMap<String, FaceEntry>,
) -> Result<Vec<FaceDef>, SetupFileError> {
    faces
        .iter()
        .map(|(name, entry)| face_def(name, entry))
        .collect()
}

/// Face definitions as a `[faces]` table.
pub(crate) fn face_entries(faces: &[FaceDef]) -> BTreeMap<String, FaceEntry> {
    faces
        .iter()
        .map(|def| {
            let mut entry = FaceEntry::default();
            match &def.reference {
                FaceRef::Name(n) => entry.name = Some(n.clone()),
                FaceRef::Id(id) => entry.id = Some(id.0),
                FaceRef::Select(s) => {
                    let f = selector_file(s);
                    entry.plane = f.plane;
                    entry.cylinder = f.cylinder;
                    entry.sphere = f.sphere;
                    entry.nearest_to = f.nearest_to;
                    entry.all_of = f.all_of;
                    entry.any_of = f.any_of;
                }
            }
            entry.expect = match def.expect {
                Expect::Exactly(1) => None,
                Expect::Exactly(n) => Some(ExpectFile::Exactly(n)),
                Expect::AtLeast(n) => Some(ExpectFile::AtLeast { at_least: n }),
            };
            (def.name.clone(), entry)
        })
        .collect()
}

/// The `[tolerance]` table, or the default.
pub(crate) fn tolerance(t: Option<ToleranceFile>) -> Tolerance {
    t.map_or_else(Tolerance::default, |t| Tolerance {
        angle_deg: t.angle_deg,
        distance: t.distance,
    })
}

/// A tolerance as a `[tolerance]` table.
pub(crate) fn tolerance_file(t: &Tolerance) -> ToleranceFile {
    ToleranceFile {
        angle_deg: t.angle_deg,
        distance: t.distance,
    }
}

/// Every option written out, so the file states what the run used.
pub(crate) fn options_file(o: &AnalysisOptions) -> OptionsFile {
    OptionsFile {
        target_cell_size: o.target_cell_size,
        max_active_cells: Some(o.max_active_cells),
        subcell_depth: Some(o.subcell_depth),
        aggregation_threshold: Some(o.aggregation_threshold),
        support_penalty_factor: Some(o.support_penalty_factor),
        nitsche_factor: Some(o.nitsche_factor),
        cg_rel_tol: Some(o.cg_rel_tol),
        cg_max_iters: Some(o.cg_max_iters),
        coarsest_dofs: Some(o.coarsest_dofs),
        grid_shift: Some(o.grid_shift),
    }
}
