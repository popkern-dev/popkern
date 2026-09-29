//! Selections as data (spec §6): a query is a value with a canonical
//! spelling, evaluated once at record time against the live body, and
//! the document stores the resolved names — never the query.

use popkern_core::naming::EntityName;
use popkern_part::{CurveKind, Dir, EdgeInfo, EdgeQuery, FaceInfo, FaceQuery, SurfaceKind};

use crate::error::ScriptError;

/// One edge predicate of the surface grammar.
#[derive(Debug, Clone, PartialEq)]
pub enum EdgePredicate {
    /// Circles and arcs.
    Circular,
    /// Straight lines.
    Straight,
    /// Conic (non-circular) curves.
    Conic,
    /// Freeform NURBS curves.
    Nurbs,
    /// Adjacent faces have these surface kinds (either order).
    Between(SurfaceKind, SurfaceKind),
    /// Midpoint height.
    AtHeight(f64),
    /// Chord direction parallel to an axis.
    Direction(Dir),
    /// An end vertex at this point (the edges meeting at one corner).
    AtVertex([f64; 3]),
    /// Chord length within a closed range.
    ChordLengthBetween(f64, f64),
    /// Not a seam edge.
    NotSeam,
    /// Exactly the entity behind an entity label.
    Labeled(String),
    /// Exactly this entity.
    Named(EntityName),
    /// Exactly the entity a spelling names (`rounds/Generated/0`,
    /// see `Model::spell_entity`); resolved at record time.
    Spelled(String),
}

/// A reference to an entity by text, resolved by the model at record
/// time: an entity label, or a spelled name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntityRef<'a> {
    /// An entity label (`labeled("mount")`).
    Label(&'a str),
    /// A spelled name (`named("rounds/Generated/0")`).
    Spelled(&'a str),
}

/// One face predicate of the surface grammar.
#[derive(Debug, Clone, PartialEq)]
pub enum FacePredicate {
    /// Surface kind.
    Kind(SurfaceKind),
    /// Outward normal along an axis.
    Facing(Dir),
    /// Outward normal against an axis.
    FacingNeg(Dir),
    /// Centroid height.
    AtHeight(f64),
    /// Exactly the entity behind an entity label.
    Labeled(String),
    /// Exactly this entity.
    Named(EntityName),
    /// Exactly the entity a spelling names.
    Spelled(String),
}

/// An edge selection: a disjunction of conjunctions.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct EdgeSelector {
    branches: Vec<Vec<EdgePredicate>>,
    allow_empty: bool,
}

/// A face selection: one conjunction.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FaceSelector {
    predicates: Vec<FacePredicate>,
    allow_empty: bool,
}

/// Start an edge selection.
pub fn edges() -> EdgeSelector {
    EdgeSelector {
        branches: vec![Vec::new()],
        allow_empty: false,
    }
}

/// Start a face selection.
pub fn faces() -> FaceSelector {
    FaceSelector::default()
}

fn kind_name(k: SurfaceKind) -> &'static str {
    match k {
        SurfaceKind::Plane => "plane",
        SurfaceKind::Cylinder => "cylinder",
        SurfaceKind::Cone => "cone",
        SurfaceKind::Sphere => "sphere",
        SurfaceKind::Torus => "torus",
        SurfaceKind::Nurbs => "nurbs",
        SurfaceKind::Gregory => "gregory",
    }
}

fn dir_name(d: Dir) -> &'static str {
    match d {
        Dir::X => "x",
        Dir::Y => "y",
        Dir::Z => "z",
    }
}

impl EdgePredicate {
    fn spell(&self) -> String {
        match self {
            EdgePredicate::Circular => "circular()".into(),
            EdgePredicate::Straight => "straight()".into(),
            EdgePredicate::Conic => "conic()".into(),
            EdgePredicate::Nurbs => "nurbs()".into(),
            EdgePredicate::Between(a, b) => {
                format!("between({}, {})", kind_name(*a), kind_name(*b))
            }
            EdgePredicate::AtHeight(z) => format!("at_height({z})"),
            EdgePredicate::Direction(d) => format!("direction({})", dir_name(*d)),
            EdgePredicate::AtVertex([x, y, z]) => format!("at_vertex([{x}, {y}, {z}])"),
            EdgePredicate::ChordLengthBetween(a, b) => format!("chord_length_between({a}, {b})"),
            EdgePredicate::NotSeam => "not_seam()".into(),
            EdgePredicate::Labeled(l) => format!("labeled({l:?})"),
            EdgePredicate::Named(n) => format!("named({n:?})"),
            EdgePredicate::Spelled(t) => format!("named({t:?})"),
        }
    }
}

impl FacePredicate {
    fn spell(&self) -> String {
        match self {
            FacePredicate::Kind(k) => format!("kind({})", kind_name(*k)),
            FacePredicate::Facing(d) => format!("facing({})", dir_name(*d)),
            FacePredicate::FacingNeg(d) => format!("facing_neg({})", dir_name(*d)),
            FacePredicate::AtHeight(z) => format!("at_height({z})"),
            FacePredicate::Labeled(l) => format!("labeled({l:?})"),
            FacePredicate::Named(n) => format!("named({n:?})"),
            FacePredicate::Spelled(t) => format!("named({t:?})"),
        }
    }
}

impl EdgeSelector {
    fn push(mut self, p: EdgePredicate) -> Self {
        if let Some(last) = self.branches.last_mut() {
            last.push(p);
        }
        self
    }
    /// Keep circular edges.
    pub fn circular(self) -> Self {
        self.push(EdgePredicate::Circular)
    }
    /// Keep straight edges.
    pub fn straight(self) -> Self {
        self.push(EdgePredicate::Straight)
    }
    /// Keep conic edges.
    pub fn conic(self) -> Self {
        self.push(EdgePredicate::Conic)
    }
    /// Keep NURBS edges.
    pub fn nurbs(self) -> Self {
        self.push(EdgePredicate::Nurbs)
    }
    /// Keep edges between these two surface kinds.
    pub fn between(self, a: SurfaceKind, b: SurfaceKind) -> Self {
        self.push(EdgePredicate::Between(a, b))
    }
    /// Keep edges whose midpoint is at height `z`.
    pub fn at_height(self, z: f64) -> Self {
        self.push(EdgePredicate::AtHeight(z))
    }
    /// Keep edges parallel to `dir`.
    pub fn direction(self, dir: Dir) -> Self {
        self.push(EdgePredicate::Direction(dir))
    }
    /// Keep edges with an end vertex at `point` (within 1e-6).
    pub fn at_vertex(self, point: [f64; 3]) -> Self {
        self.push(EdgePredicate::AtVertex(point))
    }
    /// Keep edges whose chord length is within `[min, max]`.
    pub fn chord_length_between(self, min: f64, max: f64) -> Self {
        self.push(EdgePredicate::ChordLengthBetween(min, max))
    }
    /// Drop seam edges explicitly (they are always dropped at record).
    pub fn not_seam(self) -> Self {
        self.push(EdgePredicate::NotSeam)
    }
    /// Exactly the labelled entity.
    pub fn labeled(self, label: &str) -> Self {
        self.push(EdgePredicate::Labeled(label.to_string()))
    }
    /// Exactly this entity.
    pub fn named(self, name: EntityName) -> Self {
        self.push(EdgePredicate::Named(name))
    }
    /// Exactly the entity a spelling names (`rounds/Generated/0`).
    pub fn spelled(self, text: &str) -> Self {
        self.push(EdgePredicate::Spelled(text.to_string()))
    }
    /// Union with another whole selection.
    pub fn or(mut self, other: EdgeSelector) -> Self {
        self.branches.extend(other.branches);
        self
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
                    "edges()".to_string()
                } else {
                    format!(
                        "edges().{}",
                        b.iter()
                            .map(EdgePredicate::spell)
                            .collect::<Vec<_>>()
                            .join(".")
                    )
                }
            })
            .collect();
        let mut text = parts
            .first()
            .cloned()
            .unwrap_or_else(|| "edges()".to_string());
        for part in parts.iter().skip(1) {
            text.push_str(&format!(".or({part})"));
        }
        if self.allow_empty {
            text.push_str(".allow_empty()");
        }
        text
    }
    /// Labels this selection cites.
    pub fn labels_used(&self) -> Vec<String> {
        self.branches
            .iter()
            .flatten()
            .filter_map(|p| match p {
                EdgePredicate::Labeled(l) => Some(l.clone()),
                _ => None,
            })
            .collect()
    }
    pub(crate) fn allow_empty_flag(&self) -> bool {
        self.allow_empty
    }
    /// Compile to the `popkern-part` query used for resolution.
    /// `resolve` maps an entity label or spelling to its name.
    pub(crate) fn to_query(
        &self,
        resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
    ) -> Result<EdgeQuery, ScriptError> {
        let mut query: Option<EdgeQuery> = None;
        for branch in &self.branches {
            let mut q = EdgeQuery::new();
            for p in branch {
                q = match p {
                    EdgePredicate::Circular => q.circular(),
                    EdgePredicate::Straight => q.straight(),
                    EdgePredicate::Conic => {
                        q.filter(|e: &EdgeInfo| e.curve == Some(CurveKind::Conic))
                    }
                    EdgePredicate::Nurbs => {
                        q.filter(|e: &EdgeInfo| e.curve == Some(CurveKind::Nurbs))
                    }
                    EdgePredicate::Between(a, b) => q.between(*a, *b),
                    EdgePredicate::AtHeight(z) => q.at_height(*z),
                    EdgePredicate::Direction(d) => q.direction(*d),
                    EdgePredicate::AtVertex([x, y, z]) => {
                        q.at_vertex(popkern_core::geometry::Point3::new(*x, *y, *z))
                    }
                    EdgePredicate::ChordLengthBetween(a, b) => q.chord_length_between(*a, *b),
                    EdgePredicate::NotSeam => q.not_seam(),
                    EdgePredicate::Labeled(l) => {
                        let name = resolve(EntityRef::Label(l))?;
                        q.filter(move |e: &EdgeInfo| e.name.as_ref() == Some(&name))
                    }
                    EdgePredicate::Named(n) => {
                        let name = n.clone();
                        q.filter(move |e: &EdgeInfo| e.name.as_ref() == Some(&name))
                    }
                    EdgePredicate::Spelled(t) => {
                        let name = resolve(EntityRef::Spelled(t))?;
                        q.filter(move |e: &EdgeInfo| e.name.as_ref() == Some(&name))
                    }
                };
            }
            query = Some(match query {
                None => q,
                Some(acc) => acc.or(move |_| q),
            });
        }
        let mut q = query.unwrap_or_default();
        if self.allow_empty {
            q = q.allow_empty();
        }
        Ok(q)
    }
}

impl FaceSelector {
    fn push(mut self, p: FacePredicate) -> Self {
        self.predicates.push(p);
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
    /// An empty result is a recorded no-op instead of an error.
    pub fn allow_empty(mut self) -> Self {
        self.allow_empty = true;
        self
    }
    /// Canonical spelling (spec §6.1).
    pub fn describe(&self) -> String {
        let mut text = if self.predicates.is_empty() {
            "faces()".to_string()
        } else {
            format!(
                "faces().{}",
                self.predicates
                    .iter()
                    .map(FacePredicate::spell)
                    .collect::<Vec<_>>()
                    .join(".")
            )
        };
        if self.allow_empty {
            text.push_str(".allow_empty()");
        }
        text
    }
    /// Labels this selection cites.
    pub fn labels_used(&self) -> Vec<String> {
        self.predicates
            .iter()
            .filter_map(|p| match p {
                FacePredicate::Labeled(l) => Some(l.clone()),
                _ => None,
            })
            .collect()
    }
    pub(crate) fn allow_empty_flag(&self) -> bool {
        self.allow_empty
    }
    pub(crate) fn to_query(
        &self,
        resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
    ) -> Result<FaceQuery, ScriptError> {
        let mut q = FaceQuery::new();
        for p in &self.predicates {
            q = match p {
                FacePredicate::Kind(k) => q.kind(*k),
                FacePredicate::Facing(d) => q.facing(*d),
                FacePredicate::FacingNeg(d) => q.facing_neg(*d),
                FacePredicate::AtHeight(z) => q.at_height(*z),
                FacePredicate::Labeled(l) => {
                    let name = resolve(EntityRef::Label(l))?;
                    q.filter(move |f: &FaceInfo| f.name.as_ref() == Some(&name))
                }
                FacePredicate::Named(n) => {
                    let name = n.clone();
                    q.filter(move |f: &FaceInfo| f.name.as_ref() == Some(&name))
                }
                FacePredicate::Spelled(t) => {
                    let name = resolve(EntityRef::Spelled(t))?;
                    q.filter(move |f: &FaceInfo| f.name.as_ref() == Some(&name))
                }
            };
        }
        if self.allow_empty {
            q = q.allow_empty();
        }
        Ok(q)
    }
}

#[path = "selection_resolve.rs"]
mod selection_resolve;
pub use selection_resolve::SelectionSnapshot;
pub(crate) use selection_resolve::{exactly_one, resolve_edges, resolve_faces};
