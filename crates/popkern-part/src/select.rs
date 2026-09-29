//! Edge selection: lazy predicate queries over a part's topology.
//!
//! Selectors are the load-bearing ergonomic piece of the API — they
//! replace the hand-rolled "walk every edge, check both faces'
//! surface kinds, test the midpoint height" loops that dominate the
//! hand-written examples. Queries are built inside a closure and
//! resolved when the operation runs:
//!
//! ```no_run
//! # use popkern_part::prelude::*;
//! # fn f(part: Part) -> Result<Part, PartError> {
//! part.fillet(0.08, |e| {
//!     e.circular()
//!         .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
//!         .or(|e| e.straight().at_height(0.4))
//! })
//! # }
//! ```
//!
//! An EMPTY resolved selection is an error by default (silent no-op
//! dress-ups are how scripts rot); opt out with
//! [`EdgeQuery::allow_empty`].

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::EdgeHandle;

pub use crate::select_surface_kind::SurfaceKind;

/// Cardinal directions for direction predicates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dir {
    /// Global +X (sign-insensitive in predicates).
    X,
    /// Global +Y.
    Y,
    /// Global +Z.
    Z,
}

impl Dir {
    /// The unit vector for this direction.
    pub fn vector(self) -> Vec3 {
        self.vec()
    }

    fn vec(self) -> Vec3 {
        match self {
            Dir::X => Vec3::new(1.0, 0.0, 0.0),
            Dir::Y => Vec3::new(0.0, 1.0, 0.0),
            Dir::Z => Vec3::new(0.0, 0.0, 1.0),
        }
    }
}

/// Everything a predicate can see about one edge. Also the payload
/// of the [`EdgeQuery::filter`] escape hatch.
#[derive(Debug, Clone)]
pub struct EdgeInfo {
    /// The edge.
    pub handle: EdgeHandle,
    /// Curve family of the edge geometry (None when the edge has no
    /// curve attached).
    pub curve: Option<CurveKind>,
    /// Curve midpoint in 3D.
    pub midpoint: Point3,
    /// Unit chord direction (endpoint-to-endpoint; None for closed
    /// edges, whose chord is degenerate).
    pub direction: Option<Vec3>,
    /// Surface kinds of the two adjacent faces.
    pub faces: (Option<SurfaceKind>, Option<SurfaceKind>),
    /// Whether the edge is a SEAM (both half-edges on the same face,
    /// e.g. a cylinder's closure line). Seams cannot be filleted.
    pub is_seam: bool,
    /// Endpoint-to-endpoint chord length (0 for closed edges).
    pub chord_length: f64,
    /// The edge's stable name, when the body carries provenance names.
    pub name: Option<popkern_core::naming::EntityName>,
    /// Positions of the edge's two end vertices (equal for a closed
    /// edge; `None` when the topology has no vertices for it).
    pub endpoints: Option<(Point3, Point3)>,
}

/// The analytic family of an edge's curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    /// Straight line.
    Line,
    /// Circle or circular arc.
    Circular,
    /// Other conics.
    Conic,
    /// Freeform NURBS curve.
    Nurbs,
}

type EdgePredicate = Box<dyn Fn(&EdgeInfo) -> bool>;

/// A lazy edge query: predicates AND together; whole queries union
/// via [`EdgeQuery::or`].
pub struct EdgeQuery {
    branches: Vec<Vec<EdgePredicate>>,
    description: Vec<String>,
    allow_empty: bool,
}

impl Default for EdgeQuery {
    fn default() -> Self {
        EdgeQuery::new()
    }
}

impl EdgeQuery {
    /// An empty query: matches every entity until predicates narrow it.
    pub fn new() -> EdgeQuery {
        EdgeQuery {
            branches: vec![Vec::new()],
            description: vec![String::new()],
            allow_empty: false,
        }
    }

    fn push(mut self, label: &str, pred: impl Fn(&EdgeInfo) -> bool + 'static) -> EdgeQuery {
        if let (Some(branch), Some(desc)) = (self.branches.last_mut(), self.description.last_mut())
        {
            branch.push(Box::new(pred));
            if !desc.is_empty() {
                desc.push_str(" ∧ ");
            }
            desc.push_str(label);
        }
        self
    }

    /// Keep circular edges (circles and arcs).
    pub fn circular(self) -> EdgeQuery {
        self.push("circular", |e| e.curve == Some(CurveKind::Circular))
    }

    /// Keep straight-line edges.
    pub fn straight(self) -> EdgeQuery {
        self.push("straight", |e| e.curve == Some(CurveKind::Line))
    }

    /// Keep edges whose two adjacent faces have these surface kinds
    /// (order-insensitive).
    pub fn between(self, a: SurfaceKind, b: SurfaceKind) -> EdgeQuery {
        self.push(&format!("between({a:?}, {b:?})"), move |e| {
            matches!(e.faces, (Some(x), Some(y)) if (x == a && y == b) || (x == b && y == a))
        })
    }

    /// Keep edges whose midpoint height is `z` (within `1e-6`).
    pub fn at_height(self, z: f64) -> EdgeQuery {
        self.push(&format!("at_height({z})"), move |e| {
            (e.midpoint.z - z).abs() < 1e-6
        })
    }

    /// Keep edges with an end vertex at `point` (within 1e-6, the
    /// `at_height` tolerance): the edges meeting at one corner.
    pub fn at_vertex(self, point: Point3) -> EdgeQuery {
        self.push(
            &format!("at_vertex({}, {}, {})", point.x, point.y, point.z),
            move |e| {
                e.endpoints.is_some_and(|(a, b)| {
                    a.distance_to(&point) < 1e-6 || b.distance_to(&point) < 1e-6
                })
            },
        )
    }

    /// Keep edges whose chord direction is parallel to `dir`
    /// (sign-insensitive, within ~0.6°).
    pub fn direction(self, dir: Dir) -> EdgeQuery {
        let v = dir.vec();
        self.push(&format!("direction({dir:?})"), move |e| {
            e.direction
                .map(|d| d.dot(&v).abs() > 1.0 - 1e-4)
                .unwrap_or(false)
        })
    }

    /// Keep edges whose chord length lies in `[min, max]`.
    pub fn chord_length_between(self, min: f64, max: f64) -> EdgeQuery {
        self.push(&format!("chord∈[{min},{max}]"), move |e| {
            e.chord_length >= min && e.chord_length <= max
        })
    }

    /// Exclude seam edges (both sides on one face — a cylinder's
    /// closure line). [`Part::fillet`](crate::Part::fillet) already
    /// skips seams automatically; this predicate makes the exclusion
    /// explicit in queries used for other purposes.
    pub fn not_seam(self) -> EdgeQuery {
        self.push("¬seam", |e| !e.is_seam)
    }

    /// Escape hatch: keep edges satisfying an arbitrary predicate
    /// over [`EdgeInfo`].
    pub fn filter(self, pred: impl Fn(&EdgeInfo) -> bool + 'static) -> EdgeQuery {
        self.push("filter(..)", pred)
    }

    /// Union with another query branch: edges matching EITHER side
    /// are selected.
    pub fn or(mut self, build: impl FnOnce(EdgeQuery) -> EdgeQuery) -> EdgeQuery {
        let other = build(EdgeQuery::new());
        self.branches.extend(other.branches);
        self.description.extend(other.description);
        self.allow_empty |= other.allow_empty;
        self
    }

    /// Permit the query to resolve to nothing (the operation becomes
    /// a recorded no-op instead of an error).
    pub fn allow_empty(mut self) -> EdgeQuery {
        self.allow_empty = true;
        self
    }

    /// Whether an empty resolution is permitted.
    /// Whether `.allow_empty()` was requested.
    pub fn empty_permitted(&self) -> bool {
        self.allow_empty
    }

    /// Human-readable description for feature logs and errors.
    /// The query's canonical spelling (diagnostics; the scripting surface stores it).
    pub fn describe(&self) -> String {
        self.description
            .iter()
            .map(|d| {
                if d.is_empty() {
                    "any edge".to_string()
                } else {
                    d.clone()
                }
            })
            .collect::<Vec<_>>()
            .join(" ∨ ")
    }

    /// Resolve against a body: every edge matching any branch.
    /// Resolve against precomputed edge facts; seams are NOT filtered here.
    pub fn resolve(&self, infos: &[EdgeInfo]) -> Vec<EdgeHandle> {
        infos
            .iter()
            .filter(|info| {
                self.branches
                    .iter()
                    .any(|branch| branch.iter().all(|p| p(info)))
            })
            .map(|info| info.handle)
            .collect()
    }
}

/// Everything a face predicate can see. Also the payload of the
/// [`FaceQuery::filter`] escape hatch.
#[derive(Debug, Clone)]
pub struct FaceInfo {
    /// The face.
    pub handle: popkern_core::topology::FaceHandle,
    /// Surface family.
    pub kind: Option<SurfaceKind>,
    /// Outer-loop vertex centroid.
    pub centroid: Point3,
    /// Plane normal for planar faces (chart normal, orientation
    /// applied); None otherwise.
    pub normal: Option<Vec3>,
    /// The face's stable name, when the body carries provenance names.
    pub name: Option<popkern_core::naming::EntityName>,
}

type FacePredicate = Box<dyn Fn(&FaceInfo) -> bool>;

/// A lazy face query (predicates AND together), mirroring
/// [`EdgeQuery`] for face-consuming operations (shell openings,
/// draft).
pub struct FaceQuery {
    preds: Vec<FacePredicate>,
    description: String,
    allow_empty: bool,
}

impl Default for FaceQuery {
    fn default() -> Self {
        FaceQuery::new()
    }
}

impl FaceQuery {
    /// An empty query: matches every entity until predicates narrow it.
    pub fn new() -> FaceQuery {
        FaceQuery {
            preds: Vec::new(),
            description: String::new(),
            allow_empty: false,
        }
    }

    fn push(mut self, label: &str, pred: impl Fn(&FaceInfo) -> bool + 'static) -> FaceQuery {
        self.preds.push(Box::new(pred));
        if !self.description.is_empty() {
            self.description.push_str(" ∧ ");
        }
        self.description.push_str(label);
        self
    }

    /// Keep faces of this surface family.
    pub fn kind(self, k: SurfaceKind) -> FaceQuery {
        self.push(&format!("kind({k:?})"), move |f| f.kind == Some(k))
    }

    /// Keep planar faces whose outward normal points along `dir`
    /// (signed, within ~0.8°).
    pub fn facing(self, dir: Dir) -> FaceQuery {
        let v = dir.vec();
        self.push(&format!("facing(+{dir:?})"), move |f| {
            f.normal.map(|n| n.dot(&v) > 1.0 - 1e-4).unwrap_or(false)
        })
    }

    /// Keep planar faces whose outward normal points along `-dir`.
    pub fn facing_neg(self, dir: Dir) -> FaceQuery {
        let v = dir.vec();
        self.push(&format!("facing(-{dir:?})"), move |f| {
            f.normal.map(|n| n.dot(&v) < -(1.0 - 1e-4)).unwrap_or(false)
        })
    }

    /// Keep faces whose centroid height is `z` (within 1e-6).
    pub fn at_height(self, z: f64) -> FaceQuery {
        self.push(&format!("at_height({z})"), move |f| {
            (f.centroid.z - z).abs() < 1e-6
        })
    }

    /// Escape hatch: arbitrary predicate over [`FaceInfo`].
    pub fn filter(self, pred: impl Fn(&FaceInfo) -> bool + 'static) -> FaceQuery {
        self.push("filter(..)", pred)
    }

    /// Permit an empty resolution.
    pub fn allow_empty(mut self) -> FaceQuery {
        self.allow_empty = true;
        self
    }

    /// Whether `.allow_empty()` was requested.
    pub fn empty_permitted(&self) -> bool {
        self.allow_empty
    }

    /// The query's canonical spelling (diagnostics; the scripting surface stores it).
    pub fn describe(&self) -> String {
        if self.description.is_empty() {
            "any face".into()
        } else {
            self.description.clone()
        }
    }

    /// Resolve against a body.
    /// Resolve against precomputed face facts.
    pub fn resolve(&self, infos: &[FaceInfo]) -> Vec<popkern_core::topology::FaceHandle> {
        infos
            .iter()
            .filter(|info| self.preds.iter().all(|p| p(info)))
            .map(|info| info.handle)
            .collect()
    }
}

#[path = "select_infos.rs"]
mod select_infos;
pub use select_infos::{edge_infos, face_infos};
