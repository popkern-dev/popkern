//! Parameter types for [`Feature::Generative`](crate::feature::Feature::Generative)
//! (PK-956, docs/architecture/generative-features-plan.md): the body is
//! the design space; supports, loads and kept material name its faces;
//! replay replaces it with the stiffest design at the volume fraction.

use popkern_core::naming::EntityName;
use serde::{Deserialize, Serialize};

use crate::parameter::ParameterId;

/// What a support holds (mirrors the analysis setup's support kinds a
/// script offers).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum GenerativeSupportKindDesc {
    /// Every displacement component.
    Fixed,
    /// Only the displacement along each face's own normal (radial on a
    /// bore).
    SurfaceRoller,
}

/// A support on design space faces.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenerativeSupportDesc {
    /// The faces, by stable name.
    pub faces: Vec<EntityName>,
    /// What it holds.
    pub kind: GenerativeSupportKindDesc,
}

/// A total force spread over design space faces.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenerativeLoadDesc {
    /// The faces, by stable name.
    pub faces: Vec<EntityName>,
    /// Total force, N.
    pub force: ParameterId,
    /// Direction it acts in (non-zero; normalised by the solver).
    pub direction: [f64; 3],
}

/// The build direction a design prints along: a grid axis and sense.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct BuildDirectionDesc {
    /// 0 x, 1 y, 2 z.
    pub axis: u8,
    /// Whether the part grows towards −axis.
    pub negative: bool,
}

impl BuildDirectionDesc {
    /// The skin spelling (`"+x"` … `"-z"`).
    pub fn spelling(self) -> String {
        let axis = ["x", "y", "z"]
            .get(self.axis as usize)
            .copied()
            .unwrap_or("?");
        format!("{}{axis}", if self.negative { "-" } else { "+" })
    }
}

/// A generative design's recipe.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenerativeDesc {
    /// Material preset name.
    pub material: String,
    /// Supports (at least one).
    pub supports: Vec<GenerativeSupportDesc>,
    /// Loads (at least one).
    pub loads: Vec<GenerativeLoadDesc>,
    /// Faces whose nearby material is kept, one list per selection the
    /// recipe was given.
    #[serde(default)]
    pub keep: Vec<Vec<EntityName>>,
    /// Material kept, as a fraction of the design space (a ratio).
    pub volume_fraction: ParameterId,
    /// Density filter radius (a length): members are about twice it.
    pub filter: ParameterId,
    /// How deep kept material reaches (a length; `None`: the filter
    /// radius).
    #[serde(default)]
    pub keep_depth: Option<ParameterId>,
    /// Analysis cell size (a length; `None`: a third of the filter).
    #[serde(default)]
    pub cell: Option<ParameterId>,
    /// The build direction the design prints along without supports.
    #[serde(default)]
    pub build_direction: Option<BuildDirectionDesc>,
    /// Most optimisation updates.
    pub max_iterations: u32,
    /// Thicken the thinnest members (the robust formulation).
    #[serde(default)]
    pub robust: bool,
    /// Most active analysis cells (`None`: the solver's default).
    #[serde(default)]
    pub max_cells: Option<u32>,
    /// Millimetres per stored length unit (lengths cross into the
    /// solver in millimetres; 1000 for a document storing metres).
    pub mm_per_unit: f64,
}

impl GenerativeDesc {
    /// Every parameter with its role, in a fixed order.
    pub fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        let mut v = vec![
            ("volume_fraction", self.volume_fraction),
            ("filter", self.filter),
        ];
        if let Some(id) = self.keep_depth {
            v.push(("keep_depth", id));
        }
        if let Some(id) = self.cell {
            v.push(("cell", id));
        }
        for load in &self.loads {
            v.push(("force", load.force));
        }
        v
    }

    /// Every face name the recipe references.
    pub fn entity_names(&self) -> Vec<&EntityName> {
        self.supports
            .iter()
            .flat_map(|s| s.faces.iter())
            .chain(self.loads.iter().flat_map(|l| l.faces.iter()))
            .chain(self.keep.iter().flatten())
            .collect()
    }

    /// The same, mutably (reference rebinding).
    pub fn entity_names_mut(&mut self) -> Vec<&mut EntityName> {
        self.supports
            .iter_mut()
            .flat_map(|s| s.faces.iter_mut())
            .chain(self.loads.iter_mut().flat_map(|l| l.faces.iter_mut()))
            .chain(self.keep.iter_mut().flatten())
            .collect()
    }
}

#[cfg(test)]
#[path = "feature_generative_tests.rs"]
mod tests;
