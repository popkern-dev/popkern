//! Sizing feasibility score diagnostics for corner/apex planning.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ApexFeasibilityScoreIngredient {
    ApexRangeMargin,
    DihedralMargin,
    CurvatureRadiusMargin,
    FaceAreaMargin,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct ApexFeasibilityScoreFact {
    pub apex_range_margin_nm: Option<i64>,
    pub dihedral_margin_nm: Option<i64>,
    pub curvature_radius_margin_nm: Option<i64>,
    pub face_area_margin_nm: Option<i64>,
    pub missing: ApexFeasibilityScoreMissing,
}

impl ApexFeasibilityScoreFact {
    pub fn missing_ingredients(&self) -> Vec<ApexFeasibilityScoreIngredient> {
        let mut missing = Vec::new();
        if self.missing.apex_range_margin {
            missing.push(ApexFeasibilityScoreIngredient::ApexRangeMargin);
        }
        if self.missing.dihedral_margin {
            missing.push(ApexFeasibilityScoreIngredient::DihedralMargin);
        }
        if self.missing.curvature_radius_margin {
            missing.push(ApexFeasibilityScoreIngredient::CurvatureRadiusMargin);
        }
        if self.missing.face_area_margin {
            missing.push(ApexFeasibilityScoreIngredient::FaceAreaMargin);
        }
        missing
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct ApexFeasibilityScoreMissing {
    pub apex_range_margin: bool,
    pub dihedral_margin: bool,
    pub curvature_radius_margin: bool,
    pub face_area_margin: bool,
}
