//! Where the pieces of a body touch, found from their geometry
//! (PK-1117): the points and curves `BRep::glue_touching` makes
//! topology of.
//!
//! The pieces are the face sets joined through shared edges. Between
//! two pieces:
//!
//! - a point contact where a vertex of one lies within the linear
//!   tolerance of the other: at a vertex, on an edge, or on a face, or
//!   where an edge of one touches a face of the other (`touch_points`);
//! - a curve contact where an edge of one lies on a face of the other
//!   (every sample on the surface), clipped to the runs inside the
//!   face's trim, each run's ends point contacts too; a run is kept only
//!   where it crosses a face's interior (its midpoint inside), since a
//!   run along a face's boundary needs no imprint, only its ends, and a
//!   curve grazing a boundary edge stays within the tolerance of it
//!   over a short run that is no contact curve;
//! - a tangent ruling between surfaces that touch along a line (plane
//!   and cylinder, two cylinders with parallel axes, outside or inside
//!   each other), clipped to the runs inside both faces, and a tangent
//!   point between a sphere and a plane or another sphere.
//!
//! A run's end found by sampling is moved onto the face boundary it
//! crosses. Not found: edges crossing at a point interior to both, and
//! tangency between other surface pairs (a cone on a plane, a torus);
//! those pieces stay apart where they touch.

use std::collections::BTreeMap;

use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use crate::validation::FaceUvContainment;

/// The parametric band of a sharp (even-odd) trim reading.
const SHARP: f64 = 1e-12;

/// What [`BRep::touch_plan`] found.
#[derive(Debug, Clone, Default)]
pub(super) struct TouchPlan {
    /// Point contacts, one per place.
    pub points: Vec<Point3>,
    /// Curve contacts, each running between two of `points`.
    pub curves: Vec<CurveOnEdge>,
}

/// An axis-aligned box.
#[derive(Debug, Clone, Copy)]
pub(super) struct Box3 {
    lo: Point3,
    hi: Point3,
}

impl Box3 {
    fn empty() -> Self {
        let (a, b) = (f64::INFINITY, f64::NEG_INFINITY);
        Self {
            lo: Point3::new(a, a, a),
            hi: Point3::new(b, b, b),
        }
    }

    fn grow(&mut self, p: Point3) {
        self.lo = Point3::new(self.lo.x.min(p.x), self.lo.y.min(p.y), self.lo.z.min(p.z));
        self.hi = Point3::new(self.hi.x.max(p.x), self.hi.y.max(p.y), self.hi.z.max(p.z));
    }

    fn union(mut self, o: Self) -> Self {
        self.grow(o.lo);
        self.grow(o.hi);
        self
    }

    /// Grown by `d` plus a twentieth of its diagonal on every side.
    fn padded(self, d: f64) -> Self {
        let pad = d + 0.05 * self.lo.distance_to(&self.hi);
        let v = Vec3::new(pad, pad, pad);
        Self {
            lo: self.lo - v,
            hi: self.hi + v,
        }
    }

    pub(super) fn meets(&self, o: &Self) -> bool {
        self.lo.x <= o.hi.x
            && o.lo.x <= self.hi.x
            && self.lo.y <= o.hi.y
            && o.lo.y <= self.hi.y
            && self.lo.z <= o.hi.z
            && o.lo.z <= self.hi.z
    }

    fn holds(&self, p: Point3) -> bool {
        self.meets(&Self { lo: p, hi: p })
    }

    fn corners(&self) -> [Point3; 2] {
        [self.lo, self.hi]
    }
}

/// `curve` evaluated `k` of `n` steps from its start to its end.
fn at(curve: &CurveOnEdge, k: usize, n: usize) -> f64 {
    curve.t_start + (curve.t_end - curve.t_start) * k as f64 / n as f64
}

impl BRep {
    /// Each face's piece: faces joined through a shared edge are one.
    pub(super) fn face_pieces(&self) -> BTreeMap<FaceHandle, usize> {
        let faces: Vec<FaceHandle> = self.faces.keys().collect();
        let index: BTreeMap<FaceHandle, usize> =
            faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn root(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        for e in self.edges.keys() {
            let owners: Vec<usize> = self
                .radial_fins(e)
                .into_iter()
                .filter_map(|h| self.fin_face(h))
                .filter_map(|f| index.get(&f).copied())
                .collect();
            for w in owners.windows(2) {
                let (a, b) = (root(&mut parent, w[0]), root(&mut parent, w[1]));
                parent[a] = b;
            }
        }
        let mut label: BTreeMap<usize, usize> = BTreeMap::new();
        faces
            .iter()
            .enumerate()
            .map(|(i, &f)| {
                let r = root(&mut parent, i);
                let next = label.len();
                (f, *label.entry(r).or_insert(next))
            })
            .collect()
    }

    /// A box holding `face`: its edges sampled (a sphere's whole ball).
    pub(super) fn face_box(&self, face: FaceHandle) -> Box3 {
        let mut b = Box3::empty();
        for e in self.face_edges(face) {
            if let Some(c) = self.edges[e].curve.as_ref() {
                for k in 0..=8 {
                    b.grow(c.geometry.point_at(at(c, k, 8)));
                }
            }
        }
        for v in self.face_vertices(face) {
            b.grow(self.vertices[v].position);
        }
        if let Some(SurfaceGeometry::Sphere { center, radius, .. }) =
            self.faces[face].surface.as_ref().map(|s| &s.geometry)
        {
            let r = Vec3::new(*radius, *radius, *radius);
            b.grow(*center - r);
            b.grow(*center + r);
        }
        b.padded(self.tolerance.linear())
    }

    /// Every fin of `face`'s loops.
    pub(super) fn face_fins(&self, face: FaceHandle) -> Vec<crate::topology::HalfEdgeHandle> {
        let Some(f) = self.faces.get(face) else {
            return Vec::new();
        };
        std::iter::once(f.outer_loop)
            .chain(f.inner_loops.iter().copied())
            .flat_map(|l| self.loop_cycle(l))
            .collect()
    }

    /// The edges of `face`'s boundary loops.
    pub(super) fn face_edges(&self, face: FaceHandle) -> Vec<EdgeHandle> {
        let mut out: Vec<EdgeHandle> = self
            .face_fins(face)
            .into_iter()
            .map(|h| self.half_edges[h].edge)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// The vertices of `face`: its loops' and its vertex loops'.
    pub(super) fn face_vertices(&self, face: FaceHandle) -> Vec<VertexHandle> {
        let lone = self
            .face_vertex_loops(face)
            .into_iter()
            .filter_map(|l| self.loops[l].vertex);
        let mut out: Vec<VertexHandle> = self
            .face_fins(face)
            .into_iter()
            .map(|h| self.half_edges[h].origin)
            .chain(lone)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Whether `p` lies on `face` within `tol`, inside or on its trim.
    /// The trim classifier's boundary band is parametric, ten
    /// tolerances wide (wider than the tolerance on a plane, other on a
    /// curved chart), so containment is read sharply (even-odd) and
    /// nearness to the boundary measured: inside the trim, or within
    /// `tol` of one of the face's vertices or edges.
    pub(super) fn face_holds(&self, face: FaceHandle, p: &Point3, tol: f64) -> bool {
        self.face_point_trim_containment_within(face, p, SHARP)
            .is_some_and(|c| {
                c.distance <= tol
                    && (c.containment != FaceUvContainment::Outside
                        || self.face_boundary_gap(face, p) <= tol)
            })
    }

    /// The distance from `p` to `face`'s nearest vertex or edge.
    pub(super) fn face_boundary_gap(&self, face: FaceHandle, p: &Point3) -> f64 {
        let vertices = self
            .face_vertices(face)
            .into_iter()
            .map(|v| self.vertices[v].position.distance_to(p));
        let edges = self
            .face_edges(face)
            .into_iter()
            .map(|e| self.edge_gap(e, p));
        vertices.chain(edges).fold(f64::INFINITY, f64::min)
    }

    /// The distance from `p` to `edge`'s curve within its run: the
    /// curve's closest point when it lies in the run, else the nearest
    /// of samples refined about the best.
    pub(super) fn edge_gap(&self, edge: EdgeHandle, p: &Point3) -> f64 {
        let Some(c) = self.edges.get(edge).and_then(|e| e.curve.as_ref()) else {
            return f64::INFINITY;
        };
        let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        let slack = 1e-9 * (hi - lo).abs().max(1.0);
        if let Some(u) = c.geometry.closest_point(p)
            && u >= lo - slack
            && u <= hi + slack
        {
            return p.distance_to(&c.geometry.point_at(u.clamp(lo, hi)));
        }
        let gap = |u: f64| p.distance_to(&c.geometry.point_at(u));
        let n = 32;
        let best = (0..=n)
            .map(|k| lo + (hi - lo) * k as f64 / n as f64)
            .min_by(|a, b| gap(*a).total_cmp(&gap(*b)))
            .unwrap_or(lo);
        let h = (hi - lo) / n as f64;
        gap(runs::golden_min(
            gap,
            (best - h).max(lo),
            (best + h).min(hi),
        ))
    }

    /// Whether `p` lies on `face`'s interior within `tol`.
    pub(super) fn face_interior_holds(&self, face: FaceHandle, p: &Point3, tol: f64) -> bool {
        self.face_point_trim_containment(face, p)
            .is_some_and(|c| c.distance <= tol && c.containment == FaceUvContainment::Interior)
    }

    /// Whether `p` lies within `tol` of `edge`'s interior, its ends
    /// farther.
    pub(super) fn edge_interior_holds(&self, edge: EdgeHandle, p: &Point3, tol: f64) -> bool {
        let Some((a, b)) = self.edge_vertices(edge) else {
            return false;
        };
        let Some(c) = self.edges.get(edge).and_then(|d| d.curve.as_ref()) else {
            return false;
        };
        let far = |v: VertexHandle| self.vertices[v].position.distance_to(p) > tol;
        let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
        far(a)
            && far(b)
            && c.geometry
                .closest_point(p)
                .is_some_and(|t| t > lo && t < hi && c.geometry.point_at(t).distance_to(p) <= tol)
    }

    /// Whether the piece of faces `faces` holds `p` at a vertex, on an
    /// edge or on a face, within `tol`.
    fn piece_holds(&self, faces: &[FaceHandle], p: &Point3, tol: f64) -> bool {
        faces.iter().any(|&f| {
            self.face_vertices(f)
                .into_iter()
                .any(|v| self.vertices[v].position.distance_to(p) <= tol)
                || self
                    .face_edges(f)
                    .into_iter()
                    .any(|e| self.edge_interior_holds(e, p, tol))
                || self.face_interior_holds(f, p, tol)
        })
    }

    /// Every point and curve where two pieces touch.
    pub(super) fn touch_plan(&self) -> TouchPlan {
        let tol = self.tolerance.linear();
        let pieces = self.face_pieces();
        let count = pieces.values().max().map_or(0, |m| m + 1);
        let mut faces_of: Vec<Vec<FaceHandle>> = vec![Vec::new(); count];
        for (&f, &k) in &pieces {
            faces_of[k].push(f);
        }
        let boxes: BTreeMap<FaceHandle, Box3> =
            pieces.keys().map(|&f| (f, self.face_box(f))).collect();
        let piece_box: Vec<Box3> = faces_of
            .iter()
            .map(|fs| fs.iter().fold(Box3::empty(), |b, f| b.union(boxes[f])))
            .collect();
        let mut plan = TouchPlan::default();
        for i in 0..count {
            for j in 0..count {
                if i == j || !piece_box[i].meets(&piece_box[j]) {
                    continue;
                }
                let mut verts: Vec<VertexHandle> = faces_of[i]
                    .iter()
                    .flat_map(|&f| self.face_vertices(f))
                    .collect();
                verts.sort_unstable();
                verts.dedup();
                for v in verts {
                    let p = self.vertices[v].position;
                    if piece_box[j].holds(p) && self.piece_holds(&faces_of[j], &p, tol) {
                        plan.add_point(p, tol);
                    }
                }
                let mut edges: Vec<EdgeHandle> = faces_of[i]
                    .iter()
                    .flat_map(|&f| self.face_edges(f))
                    .collect();
                edges.sort_unstable();
                edges.dedup();
                for e in edges {
                    for &g in &faces_of[j] {
                        if boxes[&g].meets(&piece_box[i]) {
                            for run in self.edge_runs_on_face(e, g) {
                                plan.add_curve(run, tol);
                            }
                        }
                    }
                }
                if i > j {
                    continue;
                }
                for &f in &faces_of[i] {
                    for &g in &faces_of[j] {
                        if boxes[&f].meets(&boxes[&g]) {
                            self.tangent_contacts(f, g, &boxes, &mut plan);
                            // An edge of one resting on the other: a
                            // ball on a box's edge (PK-1155).
                            self.edge_touch_points(f, g, false, &mut plan);
                            self.edge_touch_points(g, f, false, &mut plan);
                        }
                    }
                }
            }
        }
        // A piece touching itself: two of its faces sharing no vertex
        // tangent along a ruling or at a point (a shell's wall closed
        // exactly, its hole's offset cylinder lying on the floor's
        // offset plane), or an edge of one touching the other
        // (`touch_points`).
        for faces in &faces_of {
            for (k, &f) in faces.iter().enumerate() {
                let fv = self.face_vertices(f);
                for &g in &faces[k + 1..] {
                    if boxes[&f].meets(&boxes[&g])
                        && !self.face_vertices(g).iter().any(|v| fv.contains(v))
                    {
                        self.tangent_contacts(f, g, &boxes, &mut plan);
                        // An edge of one resting on the other (PK-1155).
                        self.edge_touch_points(f, g, true, &mut plan);
                        self.edge_touch_points(g, f, true, &mut plan);
                    }
                }
            }
        }
        plan
    }
}

impl TouchPlan {
    pub(super) fn add_point(&mut self, p: Point3, tol: f64) {
        if !self.points.iter().any(|q| q.distance_to(&p) <= tol) {
            self.points.push(p);
        }
    }

    pub(super) fn add_curve(&mut self, c: CurveOnEdge, tol: f64) {
        let ends = |c: &CurveOnEdge| {
            (
                c.geometry.point_at(c.t_start),
                c.geometry.point_at(c.t_end),
                c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
            )
        };
        let (a, b, m) = ends(&c);
        if a.distance_to(&b) <= tol {
            return;
        }
        let same = self.curves.iter().any(|d| {
            let (p, q, n) = ends(d);
            n.distance_to(&m) <= tol
                && ((p.distance_to(&a) <= tol && q.distance_to(&b) <= tol)
                    || (p.distance_to(&b) <= tol && q.distance_to(&a) <= tol))
        });
        self.add_point(a, tol);
        self.add_point(b, tol);
        if !same {
            self.curves.push(c);
        }
    }
}

#[path = "touch_points.rs"]
mod points;
#[path = "touch_runs.rs"]
mod runs;
