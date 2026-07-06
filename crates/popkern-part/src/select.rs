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

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle};

/// The analytic family of a face's surface, for adjacency selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    /// Analytic plane.
    Plane,
    /// Analytic cylinder.
    Cylinder,
    /// Analytic sphere.
    Sphere,
    /// Analytic cone.
    Cone,
    /// Analytic torus.
    Torus,
    /// Freeform NURBS surface.
    Nurbs,
    /// Gregory blend patch.
    Gregory,
}

impl SurfaceKind {
    fn of(geometry: &SurfaceGeometry) -> SurfaceKind {
        match geometry {
            SurfaceGeometry::Plane { .. } => SurfaceKind::Plane,
            SurfaceGeometry::Cylinder { .. } => SurfaceKind::Cylinder,
            SurfaceGeometry::Sphere { .. } => SurfaceKind::Sphere,
            SurfaceGeometry::Cone { .. } => SurfaceKind::Cone,
            SurfaceGeometry::Torus { .. } => SurfaceKind::Torus,
            SurfaceGeometry::Nurbs(_) => SurfaceKind::Nurbs,
            SurfaceGeometry::Gregory(_) => SurfaceKind::Gregory,
        }
    }
}

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

/// Collect the queryable facts about every edge of a body.
pub(crate) fn edge_infos(brep: &BRep) -> Vec<EdgeInfo> {
    let mut out = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let (curve, midpoint, direction, chord_length) = match edge.curve.as_ref() {
            Some(c) => {
                let ceval: &dyn CurveEval = &c.geometry;
                let mid = ceval.point_at(0.5 * (c.t_start + c.t_end));
                let p0 = ceval.point_at(c.t_start);
                let p1 = ceval.point_at(c.t_end);
                let chord = p1 - p0;
                let len = chord.length();
                let dir = if len > 1e-9 {
                    Some(chord * (1.0 / len))
                } else {
                    None
                };
                let kind = match &c.geometry {
                    CurveGeometry::Line { .. } => CurveKind::Line,
                    CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. } => CurveKind::Circular,
                    CurveGeometry::Ellipse { .. }
                    | CurveGeometry::Parabola { .. }
                    | CurveGeometry::Hyperbola { .. } => CurveKind::Conic,
                    CurveGeometry::Nurbs(_) => CurveKind::Nurbs,
                };
                (Some(kind), mid, dir, len)
            }
            None => {
                let he = edge.half_edge;
                let p = brep
                    .half_edges
                    .get(he)
                    .and_then(|h| brep.vertices.get(h.origin))
                    .map(|v| v.position)
                    .unwrap_or(Point3::new(0.0, 0.0, 0.0));
                (None, p, None, 0.0)
            }
        };
        let he = edge.half_edge;
        let (faces, is_seam) = match brep.half_edges.get(he) {
            Some(h) => {
                let f1 = brep.loops.get(h.loop_).map(|l| l.face);
                let f2 = brep
                    .half_edges
                    .get(h.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| l.face);
                let kind_of = |f: Option<popkern_core::topology::FaceHandle>| {
                    f.and_then(|fh| brep.faces.get(fh))
                        .and_then(|face| face.surface.as_ref())
                        .map(|s| SurfaceKind::of(&s.geometry))
                };
                let seam = matches!((f1, f2), (Some(a), Some(b)) if a == b);
                ((kind_of(f1), kind_of(f2)), seam)
            }
            None => ((None, None), false),
        };
        out.push(EdgeInfo {
            handle: eh,
            curve,
            midpoint,
            direction,
            faces,
            is_seam,
            chord_length,
        });
    }
    out
}

type EdgePredicate = Box<dyn Fn(&EdgeInfo) -> bool>;

/// A lazy edge query: predicates AND together; whole queries union
/// via [`EdgeQuery::or`].
pub struct EdgeQuery {
    branches: Vec<Vec<EdgePredicate>>,
    description: Vec<String>,
    allow_empty: bool,
}

impl EdgeQuery {
    pub(crate) fn new() -> EdgeQuery {
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
    pub(crate) fn empty_permitted(&self) -> bool {
        self.allow_empty
    }

    /// Human-readable description for feature logs and errors.
    pub(crate) fn describe(&self) -> String {
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
    pub(crate) fn resolve(&self, infos: &[EdgeInfo]) -> Vec<EdgeHandle> {
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
