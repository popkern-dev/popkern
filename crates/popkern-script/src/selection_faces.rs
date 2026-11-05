//! Face selection builders, spelling and compilation (split from
//! `selection` under the file budget). A face selection is a
//! disjunction of conjunctions: `or` joins whole selections,
//! `except` (and `!`, `faces().except(self)`) removes what another
//! selection matches, so "every face but the bottom" is
//! `faces().except(faces().facing_neg("z"))`.

use popkern_core::naming::EntityName;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_history::unit::Unit;
use popkern_part::{Dir, FaceInfo, FaceQuery, SurfaceKind};

use super::{EntityRef, FacePredicate, FaceSelector, faces_named};
use crate::error::ScriptError;

impl FaceSelector {
    fn push(mut self, p: FacePredicate) -> Self {
        if let Some(last) = self.branches.last_mut() {
            last.push(p);
        }
        self
    }
    /// Keep faces of this surface kind.
    pub fn kind(self, k: SurfaceKind) -> Self {
        self.push(FacePredicate::Kind(k))
    }
    /// Keep faces whose outward normal is along `dir`.
    pub fn facing(self, dir: Dir) -> Self {
        self.push(FacePredicate::Facing(dir))
    }
    /// Keep faces whose outward normal is against `dir`.
    pub fn facing_neg(self, dir: Dir) -> Self {
        self.push(FacePredicate::FacingNeg(dir))
    }
    /// Keep faces whose centroid is at height `z`.
    pub fn at_height(self, z: f64) -> Self {
        self.push(FacePredicate::AtHeight(z))
    }
    /// Exactly the labelled entity.
    pub fn labeled(self, label: &str) -> Self {
        self.push(FacePredicate::Labeled(label.to_string()))
    }
    /// Exactly this entity.
    pub fn named(self, name: EntityName) -> Self {
        self.push(FacePredicate::Named(name))
    }
    /// Exactly the entity a spelling names.
    pub fn spelled(self, text: &str) -> Self {
        self.push(FacePredicate::Spelled(text.to_string()))
    }
    /// Union with another whole selection.
    pub fn or(mut self, other: FaceSelector) -> Self {
        self.branches.extend(other.branches);
        self
    }
    /// Drop the faces `other` matches (from the last branch).
    pub fn except(self, other: FaceSelector) -> Self {
        self.push(FacePredicate::Except(Box::new(other)))
    }
    /// An empty result is a recorded no-op instead of an error.
    pub fn allow_empty(mut self) -> Self {
        self.allow_empty = true;
        self
    }
    /// Canonical spelling (spec §6.1).
    pub fn describe(&self) -> String {
        let parts: Vec<String> = self
            .branches
            .iter()
            .map(|b| {
                if b.is_empty() {
                    "faces()".to_string()
                } else {
                    format!(
                        "faces().{}",
                        b.iter()
                            .map(FacePredicate::spell)
                            .collect::<Vec<_>>()
                            .join(".")
                    )
                }
            })
            .collect();
        let mut text = parts
            .first()
            .cloned()
            .unwrap_or_else(|| "faces()".to_string());
        for part in parts.iter().skip(1) {
            text.push_str(&format!(".or({part})"));
        }
        if self.allow_empty {
            text.push_str(".allow_empty()");
        }
        text
    }
    /// Labels this selection cites (inner `except` selections included).
    pub fn labels_used(&self) -> Vec<String> {
        self.branches
            .iter()
            .flatten()
            .flat_map(|p| match p {
                FacePredicate::Labeled(l) => vec![l.clone()],
                FacePredicate::Except(inner) => inner.labels_used(),
                _ => Vec::new(),
            })
            .collect()
    }
    pub(crate) fn allow_empty_flag(&self) -> bool {
        self.allow_empty
    }
    /// The faces of `brep` this selection matches, in face order
    /// (empty allowed: the caller judges emptiness).
    fn matching(
        &self,
        brep: &BRep,
        infos: &[FaceInfo],
        resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
        unit: Unit,
    ) -> Result<std::collections::HashSet<FaceHandle>, ScriptError> {
        let mut out = std::collections::HashSet::new();
        for branch in &self.branches {
            out.extend(branch_query(branch, brep, infos, resolve, unit)?.resolve(infos));
        }
        Ok(out)
    }
    /// Compile to the `popkern-part` query used for resolution against
    /// `brep` (names follow lineage, as for edges).
    /// `unit` as for [`super::EdgeSelector::to_query`].
    pub(crate) fn to_query(
        &self,
        brep: &BRep,
        resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
        unit: Unit,
    ) -> Result<FaceQuery, ScriptError> {
        let infos = popkern_part::face_infos(brep);
        let mut q = match self.branches.as_slice() {
            [only] => branch_query(only, brep, &infos, resolve, unit)?,
            _ => {
                let set = self.matching(brep, &infos, resolve, unit)?;
                FaceQuery::new().filter(move |f: &FaceInfo| set.contains(&f.handle))
            }
        };
        if self.allow_empty {
            q = q.allow_empty();
        }
        Ok(q)
    }
}

/// Every face a selection does not match: `faces().except(self)`.
impl std::ops::Not for FaceSelector {
    type Output = FaceSelector;
    fn not(self) -> FaceSelector {
        let allow_empty = self.allow_empty;
        let mut inner = self;
        inner.allow_empty = false;
        let mut out = FaceSelector::default().except(inner);
        out.allow_empty = allow_empty;
        out
    }
}

/// One conjunction as a `popkern-part` query.
fn branch_query(
    branch: &[FacePredicate],
    brep: &BRep,
    infos: &[FaceInfo],
    resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
    unit: Unit,
) -> Result<FaceQuery, ScriptError> {
    let mut q = FaceQuery::new();
    for p in branch {
        q = match p {
            FacePredicate::Kind(k) => q.kind(*k),
            FacePredicate::Facing(d) => q.facing(*d),
            FacePredicate::FacingNeg(d) => q.facing_neg(*d),
            FacePredicate::AtHeight(z) => q.at_height(unit.to_si(*z)),
            FacePredicate::Labeled(l) => {
                let set = faces_named(brep, &resolve(EntityRef::Label(l))?);
                q.filter(move |f: &FaceInfo| set.contains(&f.handle))
            }
            FacePredicate::Named(n) => {
                let set = faces_named(brep, n);
                q.filter(move |f: &FaceInfo| set.contains(&f.handle))
            }
            FacePredicate::Spelled(t) => {
                let set = faces_named(brep, &resolve(EntityRef::Spelled(t))?);
                q.filter(move |f: &FaceInfo| set.contains(&f.handle))
            }
            FacePredicate::Except(inner) => {
                let set = inner.matching(brep, infos, resolve, unit)?;
                q.filter(move |f: &FaceInfo| !set.contains(&f.handle))
            }
        };
    }
    Ok(q)
}
