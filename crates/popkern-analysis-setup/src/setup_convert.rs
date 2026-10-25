//! Conversion between the setup file shape and the solver's setup
//! types, with typed errors for everything the file can get wrong.

use popkern_analysis::setup::{GenerativeSpec, InfillSpec, TemperatureSource, ThermalSpec};
use popkern_analysis::structural::generative::print::BuildDirection;
use popkern_analysis::{
    AnalysisOptions, Axis, Expect, FaceDef, FaceId, FaceRef, LoadDirection, LoadSpec, Selector,
    StructuralSetup, SupportKind, SupportSpec,
};

pub use super::error::SetupFileError;
use super::shared::{
    check_analysis, check_units, face_defs, face_entries, options_file, required, tolerance,
    tolerance_file,
};
use super::{
    CylinderFile, DirectionFile, ExpectFile, FaceEntry, GenerativeFile, InfillFile, LoadFile,
    ModalFile, OptionsFile, PlaneFile, SelectorFile, SetupFile, SphereFile, SupportFile,
    SupportKindFile, ThermalLoadFile,
};

fn selector(face: &str, s: &SelectorFile) -> Result<Selector, SetupFileError> {
    let found = usize::from(s.plane.is_some())
        + usize::from(s.cylinder.is_some())
        + usize::from(s.sphere.is_some())
        + usize::from(s.nearest_to.is_some())
        + usize::from(s.all_of.is_some())
        + usize::from(s.any_of.is_some());
    if found != 1 {
        return Err(SetupFileError::FaceReferenceCount {
            face: face.to_string(),
            found,
        });
    }
    let list = |parts: &[SelectorFile], tag: &str| -> Result<Vec<Selector>, SetupFileError> {
        let nested = format!("{face}/{tag}");
        parts.iter().map(|p| selector(&nested, p)).collect()
    };
    Ok(if let Some(p) = &s.plane {
        Selector::Plane {
            normal: p.normal,
            through: p.through,
            angle_tol: p.angle_tol,
            dist_tol: p.dist_tol,
        }
    } else if let Some(c) = &s.cylinder {
        Selector::Cylinder {
            axis: c.axis,
            radius: c.radius,
            through: c.through,
            angle_tol: c.angle_tol,
            dist_tol: c.dist_tol,
        }
    } else if let Some(sp) = &s.sphere {
        Selector::Sphere {
            centre: sp.centre,
            radius: sp.radius,
            dist_tol: sp.dist_tol,
        }
    } else if let Some(point) = s.nearest_to {
        Selector::NearestTo(point)
    } else if let Some(parts) = &s.all_of {
        Selector::AllOf(list(parts, "all_of")?)
    } else {
        Selector::AnyOf(list(s.any_of.as_deref().unwrap_or(&[]), "any_of")?)
    })
}

pub(crate) fn selector_file(s: &Selector) -> SelectorFile {
    let mut f = SelectorFile::default();
    match s {
        Selector::Plane {
            normal,
            through,
            angle_tol,
            dist_tol,
        } => {
            f.plane = Some(PlaneFile {
                normal: *normal,
                through: *through,
                angle_tol: *angle_tol,
                dist_tol: *dist_tol,
            })
        }
        Selector::Cylinder {
            axis,
            radius,
            through,
            angle_tol,
            dist_tol,
        } => {
            f.cylinder = Some(CylinderFile {
                axis: *axis,
                radius: *radius,
                through: *through,
                angle_tol: *angle_tol,
                dist_tol: *dist_tol,
            })
        }
        Selector::Sphere {
            centre,
            radius,
            dist_tol,
        } => {
            f.sphere = Some(SphereFile {
                centre: *centre,
                radius: *radius,
                dist_tol: *dist_tol,
            })
        }
        Selector::NearestTo(p) => f.nearest_to = Some(*p),
        Selector::AllOf(parts) => f.all_of = Some(parts.iter().map(selector_file).collect()),
        Selector::AnyOf(parts) => f.any_of = Some(parts.iter().map(selector_file).collect()),
    }
    f
}

pub(crate) fn face_def(name: &str, e: &FaceEntry) -> Result<FaceDef, SetupFileError> {
    let as_selector = SelectorFile {
        plane: e.plane,
        cylinder: e.cylinder,
        sphere: e.sphere,
        nearest_to: e.nearest_to,
        all_of: e.all_of.clone(),
        any_of: e.any_of.clone(),
    };
    let selectors = [
        e.plane.is_some(),
        e.cylinder.is_some(),
        e.sphere.is_some(),
        e.nearest_to.is_some(),
        e.all_of.is_some(),
        e.any_of.is_some(),
    ];
    let selector_keys = selectors.iter().filter(|&&b| b).count();
    let found = selector_keys + usize::from(e.name.is_some()) + usize::from(e.id.is_some());
    if found != 1 {
        return Err(SetupFileError::FaceReferenceCount {
            face: name.to_string(),
            found,
        });
    }
    let reference = match (&e.name, e.id) {
        (Some(n), _) => FaceRef::Name(n.clone()),
        (None, Some(id)) => FaceRef::Id(FaceId(id)),
        (None, None) => FaceRef::Select(selector(name, &as_selector)?),
    };
    let expect = match e.expect {
        None => Expect::default(),
        Some(ExpectFile::Exactly(n)) => Expect::Exactly(n),
        Some(ExpectFile::AtLeast { at_least }) => Expect::AtLeast(at_least),
    };
    Ok(FaceDef {
        name: name.to_string(),
        reference,
        expect,
    })
}

const AXES: [(&str, Axis); 6] = [
    ("+x", Axis::PlusX),
    ("-x", Axis::MinusX),
    ("+y", Axis::PlusY),
    ("-y", Axis::MinusY),
    ("+z", Axis::PlusZ),
    ("-z", Axis::MinusZ),
];

fn direction(load: usize, d: &DirectionFile) -> Result<LoadDirection, SetupFileError> {
    match d {
        DirectionFile::Vector(v) => Ok(LoadDirection::Vector(*v)),
        DirectionFile::Keyword(k) if k == "face_normal" => Ok(LoadDirection::FaceNormal),
        DirectionFile::Keyword(k) => AXES
            .iter()
            .find(|(name, _)| name == k)
            .map(|&(_, axis)| LoadDirection::Axis(axis))
            .ok_or_else(|| SetupFileError::UnknownDirection {
                load,
                value: k.clone(),
            }),
    }
}

/// The signed axes a part may grow along, by name.
const BUILD_AXES: [(&str, BuildDirection); 6] = [
    (
        "+x",
        BuildDirection {
            axis: 0,
            negative: false,
        },
    ),
    (
        "-x",
        BuildDirection {
            axis: 0,
            negative: true,
        },
    ),
    (
        "+y",
        BuildDirection {
            axis: 1,
            negative: false,
        },
    ),
    (
        "-y",
        BuildDirection {
            axis: 1,
            negative: true,
        },
    ),
    (
        "+z",
        BuildDirection {
            axis: 2,
            negative: false,
        },
    ),
    (
        "-z",
        BuildDirection {
            axis: 2,
            negative: true,
        },
    ),
];

fn build_direction(name: &str) -> Result<BuildDirection, SetupFileError> {
    BUILD_AXES
        .iter()
        .find(|(n, _)| *n == name)
        .map(|&(_, d)| d)
        .ok_or_else(|| SetupFileError::UnknownBuildDirection(name.to_string()))
}

fn build_direction_name(d: BuildDirection) -> String {
    BUILD_AXES
        .iter()
        .find(|(_, b)| *b == d)
        .map_or("+z", |(n, _)| n)
        .to_string()
}

fn direction_file(d: &LoadDirection) -> DirectionFile {
    match d {
        LoadDirection::Vector(v) => DirectionFile::Vector(*v),
        LoadDirection::FaceNormal => DirectionFile::Keyword("face_normal".to_string()),
        LoadDirection::Axis(axis) => {
            let name = AXES
                .iter()
                .find(|(_, a)| a == axis)
                .map_or("+z", |(n, _)| n);
            DirectionFile::Keyword(name.to_string())
        }
    }
}

pub(crate) fn options(o: Option<&OptionsFile>) -> AnalysisOptions {
    let d = AnalysisOptions::default();
    let Some(o) = o else { return d };
    AnalysisOptions {
        target_cell_size: o.target_cell_size.or(d.target_cell_size),
        max_active_cells: o.max_active_cells.unwrap_or(d.max_active_cells),
        subcell_depth: o.subcell_depth.unwrap_or(d.subcell_depth),
        aggregation_threshold: o.aggregation_threshold.unwrap_or(d.aggregation_threshold),
        support_penalty_factor: o.support_penalty_factor.unwrap_or(d.support_penalty_factor),
        nitsche_factor: o.nitsche_factor.unwrap_or(d.nitsche_factor),
        cg_rel_tol: o.cg_rel_tol.unwrap_or(d.cg_rel_tol),
        cg_max_iters: o.cg_max_iters.unwrap_or(d.cg_max_iters),
        coarsest_dofs: o.coarsest_dofs.unwrap_or(d.coarsest_dofs),
        grid_shift: o.grid_shift.unwrap_or(d.grid_shift),
    }
}

impl SetupFile {
    /// Converts to solver setup types. Validation of values is left to
    /// `StructuralSetup::validate`.
    pub fn into_setup(self) -> Result<StructuralSetup, SetupFileError> {
        check_analysis(self.analysis.as_deref(), "structural")?;
        check_units(self.units.as_deref())?;
        let material = required("material", self.material)?;
        let faces = face_defs(&self.faces)?;
        let supports = self
            .support
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let kind = match (s.kind, s.normal, s.direction) {
                    (SupportKindFile::Fixed, None, None) => SupportKind::Fixed,
                    (SupportKindFile::Fixed, _, _) => {
                        return Err(SetupFileError::NormalOnFixed(i));
                    }
                    (SupportKindFile::Roller, Some(normal), None) => SupportKind::Roller { normal },
                    (SupportKindFile::Roller, None, Some(direction)) => {
                        SupportKind::DirectionalRoller { direction }
                    }
                    (SupportKindFile::Roller, None, None) => SupportKind::SurfaceRoller,
                    (SupportKindFile::Roller, Some(_), Some(_)) => {
                        return Err(SetupFileError::RollerNormalAndDirection(i));
                    }
                };
                Ok(SupportSpec {
                    faces: s.faces.clone(),
                    kind,
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let loads = self
            .load
            .iter()
            .enumerate()
            .map(|(i, l)| {
                Ok(LoadSpec {
                    faces: l.faces.clone(),
                    total_force_n: l.force_n,
                    direction: direction(i, &l.direction)?,
                })
            })
            .collect::<Result<Vec<_>, SetupFileError>>()?;
        let thermal = match &self.thermal {
            None => None,
            Some(t) => Some(ThermalSpec {
                reference_c: t.reference_c,
                temperature: match (t.uniform_c, &t.setup) {
                    (Some(c), None) => TemperatureSource::Uniform(c),
                    (None, Some(path)) => TemperatureSource::ThermalSetup(path.clone()),
                    _ => return Err(SetupFileError::ThermalSourceCount),
                },
            }),
        };
        Ok(StructuralSetup {
            material,
            build_direction: self.build_direction,
            tolerance: tolerance(self.tolerance),
            faces,
            supports,
            loads,
            thermal,
            modal: self.modal.map(|m| m.modes),
            infill: self.infill.as_ref().map(|i| InfillSpec {
                pattern: i.pattern.clone(),
                density: i.density,
                walls_mm: i.walls_mm,
                line_width_mm: i.line_width_mm,
            }),
            generative: match &self.generative {
                Some(g) => Some(GenerativeSpec {
                    volume_fraction: g.volume_fraction,
                    filter_mm: g.filter_mm,
                    keep: g.keep.clone(),
                    keep_mm: g.keep_mm,
                    max_iterations: g.max_iterations,
                    crisp: g.crisp,
                    robust: g.robust,
                    build_direction: g
                        .build_direction
                        .as_deref()
                        .map(build_direction)
                        .transpose()?,
                }),
                None => None,
            },
            options: options(self.options.as_ref()),
        })
    }

    /// The file shape of a setup. Every option and the tolerance are
    /// written out, so the file states exactly what the run used.
    pub fn from_setup(setup: &StructuralSetup) -> Self {
        SetupFile {
            analysis: None,
            units: Some("mm".to_string()),
            material: Some(setup.material.clone()),
            build_direction: setup.build_direction,
            tolerance: Some(tolerance_file(&setup.tolerance)),
            faces: face_entries(&setup.faces),
            support: setup
                .supports
                .iter()
                .map(|s| match s.kind {
                    SupportKind::Fixed => SupportFile {
                        faces: s.faces.clone(),
                        kind: SupportKindFile::Fixed,
                        normal: None,
                        direction: None,
                    },
                    SupportKind::Roller { normal } => SupportFile {
                        faces: s.faces.clone(),
                        kind: SupportKindFile::Roller,
                        normal: Some(normal),
                        direction: None,
                    },
                    SupportKind::SurfaceRoller => SupportFile {
                        faces: s.faces.clone(),
                        kind: SupportKindFile::Roller,
                        normal: None,
                        direction: None,
                    },
                    SupportKind::DirectionalRoller { direction } => SupportFile {
                        faces: s.faces.clone(),
                        kind: SupportKindFile::Roller,
                        normal: None,
                        direction: Some(direction),
                    },
                })
                .collect(),
            load: setup
                .loads
                .iter()
                .map(|l| LoadFile {
                    faces: l.faces.clone(),
                    force_n: l.total_force_n,
                    direction: direction_file(&l.direction),
                })
                .collect(),
            thermal: setup.thermal.as_ref().map(|t| ThermalLoadFile {
                reference_c: t.reference_c,
                uniform_c: match t.temperature {
                    TemperatureSource::Uniform(c) => Some(c),
                    TemperatureSource::ThermalSetup(_) => None,
                },
                setup: match &t.temperature {
                    TemperatureSource::ThermalSetup(p) => Some(p.clone()),
                    TemperatureSource::Uniform(_) => None,
                },
            }),
            modal: setup.modal.map(|modes| ModalFile { modes }),
            infill: setup.infill.as_ref().map(|i| InfillFile {
                pattern: i.pattern.clone(),
                density: i.density,
                walls_mm: i.walls_mm,
                line_width_mm: i.line_width_mm,
            }),
            generative: setup.generative.as_ref().map(|g| GenerativeFile {
                volume_fraction: g.volume_fraction,
                filter_mm: g.filter_mm,
                keep: g.keep.clone(),
                keep_mm: g.keep_mm,
                max_iterations: g.max_iterations,
                crisp: g.crisp,
                robust: g.robust,
                build_direction: g.build_direction.map(build_direction_name),
            }),
            options: Some(options_file(&setup.options)),
        }
    }
}
