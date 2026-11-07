//! glue_faces / unglue_face: two coincident faces facing each other
//! become one face between the regions behind them (a box standing on
//! a box), and an internal face splits back into two.
//!
//! The glue is one step over the face pair, its edges and the space
//! between: gluing only the boundary edges of two touching faces would
//! close the zero-thickness gap between them off as a region of its own,
//! whose sense no geometry decides. The faces must share their vertices
//! ([`BRep::glue_vertices`]); their edges may be one already or still
//! two coincident edges, which merge. The kept face takes the far side
//! of the glued-away one, so each of its sides now faces a region
//! behind one of the faces, and the shells of the region between lose
//! those two sides (a shell left in pieces separates; a bounded region
//! left in pieces declines).

use std::collections::{BTreeMap, BTreeSet};

use super::EulerError;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;
use crate::topology::{
    BRep, EdgeHandle, Face, FaceHandle, FaceSide, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
    RadialLinks, Region, RegionKind, Shell, ShellHandle,
};

/// What [`BRep::glue_faces`] changed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GlueFaceResult {
    /// The face kept, now between the two regions.
    pub face: FaceHandle,
    /// Coincident edges of the glued-away face that merged into the
    /// kept face's.
    pub merged_edges: Vec<EdgeHandle>,
    /// Shells that went: joined into another of their region, or left
    /// holding nothing.
    pub removed_shells: Vec<ShellHandle>,
}

impl BRep {
    /// Glue face `gone` onto coincident face `keep` (see the module
    /// docs). Counts: faces −1, loops −(gone's loops), fins −(gone's
    /// fins), edges −(edges merged), shells −(shells joined or emptied).
    ///
    /// # Errors
    ///
    /// [`EulerError::FacesDoNotMatch`] when the boundaries do not pair
    /// fin for fin, running opposite ways along coincident edges;
    /// [`EulerError::NotCoincident`] when the surfaces part inside;
    /// [`EulerError::FacesNotFacing`] when no side of each lies in one
    /// region; [`EulerError::RadialOrder`] / [`EulerError::RegionConflict`]
    /// from merging an edge's rings; [`EulerError::UnglueSplitsRegion`]
    /// when the region between would fall into pieces. The body is then
    /// unchanged.
    pub fn glue_faces(
        &mut self,
        keep: FaceHandle,
        gone: FaceHandle,
    ) -> Result<GlueFaceResult, EulerError> {
        if !self.faces.contains_key(keep) || !self.faces.contains_key(gone) || keep == gone {
            return Err(EulerError::FacesDoNotMatch { keep, gone });
        }
        let pairs = self.pair_face_fins(keep, gone)?;
        let tolerance = self.tolerance.linear();
        let bulge = self.surface_gap(keep, gone);
        if bulge > tolerance {
            return Err(EulerError::NotCoincident {
                distance: bulge,
                tolerance,
            });
        }
        let (k, g) = (self.faces[keep].clone(), self.faces[gone].clone());
        let region = |s: ShellHandle| self.shells.get(s).map(|s| s.region);
        // The sides facing each other lie in one region; the face keeps
        // its other side and takes the glued-away face's other one.
        let fronts = region(k.front).is_some() && region(k.front) == region(g.front);
        let backs = region(k.back).is_some() && region(k.back) == region(g.back);
        let between = match (fronts, backs) {
            (true, _) => FaceSide::Front,
            (false, true) => FaceSide::Back,
            _ => return Err(EulerError::FacesNotFacing { keep, gone }),
        };
        let mut work = self.clone();
        let mut merged_edges = Vec::new();
        let mut edges = BTreeSet::new();
        for &(kf, gf) in &pairs {
            let (ek, eg) = (work.half_edges[kf].edge, work.half_edges[gf].edge);
            edges.insert(ek);
            let rest: Vec<_> = work
                .radial_fins(eg)
                .into_iter()
                .filter(|&f| f != gf)
                .collect();
            work.keep_first_fin_off(eg, gf);
            work.remove_half_edge_reanchored(gf);
            if ek != eg {
                for &f in &rest {
                    work.half_edges[f].edge = ek;
                }
                let mut ring = work.radial_fins(ek);
                ring.extend(rest);
                let order = work
                    .radial_order_with_sense(ek, ring, &self.tolerance)
                    .map_err(|e| EulerError::RadialOrder(Box::new(e)))?;
                let ring: Vec<_> = order.iter().map(|&(f, _)| f).collect();
                work.link_ring(&ring);
                work.edges.remove(eg);
                merged_edges.push(eg);
            }
        }
        for l in std::iter::once(g.outer_loop).chain(g.inner_loops.iter().copied()) {
            work.loops.remove(l);
        }
        work.faces.remove(gone);
        work.face_use_pcurves.retain(|key, _| key.face != gone);
        let old: BTreeSet<ShellHandle> = match between {
            FaceSide::Front => {
                work.faces[keep].front = g.back;
                [k.front, g.front].into()
            }
            FaceSide::Back => {
                work.faces[keep].back = g.front;
                [k.back, g.back].into()
            }
        };
        let mut joins = Vec::new();
        let mut adopt = Vec::new();
        for &e in &edges {
            let order = work
                .radial_order_with_sense(e, work.radial_fins(e), &self.tolerance)
                .map_err(|err| EulerError::RadialOrder(Box::new(err)))?;
            let wedge = work.wedge_joins(e, &order)?;
            joins.extend(wedge.joins);
            adopt.extend(wedge.adopt);
        }
        work.adopt_shells(&adopt);
        let mut removed_shells = work.join_shells(joins);
        for shell in old {
            if !work.shells.contains_key(shell) {
                continue;
            }
            let region = work.shells[shell].region;
            if !work.shell_has_faces(shell)
                && work.shells[shell].wire_edges.is_empty()
                && work.shells[shell].acorn.is_none()
            {
                work.shells.remove(shell);
                removed_shells.push(shell);
                if region != work.outside && work.region_shells(region).is_empty() {
                    work.regions.remove(region);
                } else if work.regions.get(region).is_some_and(|r| r.shell == shell) {
                    return Err(EulerError::UnglueSplitsRegion { region });
                }
                continue;
            }
            if !work.separate_shell_pieces(shell).is_empty() && region != work.outside {
                return Err(EulerError::UnglueSplitsRegion { region });
            }
        }
        *self = work;
        Ok(GlueFaceResult {
            face: keep,
            merged_edges,
            removed_shells,
        })
    }

    /// Split internal face `face` (a solid region behind each side, in
    /// different shells) into
    /// itself and a reversed twin on the same edges, with a new void
    /// region of no thickness between them, the inverse of
    /// [`BRep::glue_faces`] up to that region (gluing the two again
    /// restores the face). Returns the twin. Counts: faces +1, loops +
    /// (its loops), fins + (its fins), shells +1, regions +1.
    ///
    /// # Errors
    ///
    /// [`EulerError::NotAnInternalFace`] when both sides lie in one
    /// shell or a fin's sense along its edge is undecidable.
    pub fn unglue_face(&mut self, face: FaceHandle) -> Result<FaceHandle, EulerError> {
        let data = self
            .faces
            .get(face)
            .cloned()
            .ok_or(EulerError::NotAnInternalFace { face })?;
        let solid = |s: ShellHandle| {
            self.shell_region(s).and_then(|r| self.region_kind(r)) == Some(RegionKind::Solid)
        };
        if data.back == data.front || !solid(data.back) || !solid(data.front) {
            return Err(EulerError::NotAnInternalFace { face });
        }
        let loops: Vec<LoopHandle> = std::iter::once(data.outer_loop)
            .chain(data.inner_loops.iter().copied())
            .collect();
        let fins: Vec<Vec<HalfEdgeHandle>> = loops.iter().map(|&l| self.loop_cycle(l)).collect();
        let senses: BTreeMap<HalfEdgeHandle, bool> = fins
            .iter()
            .flatten()
            .map(|&f| self.fin_is_forward(f).map(|s| (f, s)))
            .collect::<Option<_>>()
            .ok_or(EulerError::NotAnInternalFace { face })?;
        let mut work = self.clone();
        let region = work.regions.insert(Region {
            kind: RegionKind::Void,
            shell: ShellHandle::default(),
        });
        let gap = work.shells.insert(Shell::new(face, region));
        work.regions[region].shell = gap;
        let mut surface = data.surface.clone();
        if let Some(s) = surface.as_mut() {
            s.orientation = !s.orientation;
        }
        let twin = work.faces.insert(Face {
            outer_loop: LoopHandle::default(),
            inner_loops: Vec::new(),
            back: data.front,
            front: gap,
            surface,
            name: None,
            tolerance: data.tolerance,
            role: data.role,
        });
        work.faces[face].front = gap;
        if work.shells[data.front].face == face {
            work.shells[data.front].face = twin;
        }
        for (i, cycle) in fins.iter().enumerate() {
            let new_loop = work.loops.insert(Loop {
                half_edge: HalfEdgeHandle::default(),
                face: twin,
                vertex: None,
            });
            if i == 0 {
                work.faces[twin].outer_loop = new_loop;
            } else {
                work.faces[twin].inner_loops.push(new_loop);
            }
            // The twin runs each fin backwards, the loop reversed.
            let mut made = Vec::with_capacity(cycle.len());
            for &k in cycle {
                let end = work.half_edges[work.half_edges[k].next].origin;
                let edge = work.half_edges[k].edge;
                let pcurve = work.half_edges[k].pcurve.as_ref().map(|p| p.reversed());
                let g = work.half_edges.insert(HalfEdge {
                    origin: end,
                    ring: RadialLinks::alone(),
                    next: HalfEdgeHandle::default(),
                    prev: HalfEdgeHandle::default(),
                    edge,
                    loop_: new_loop,
                    pcurve,
                });
                // The twin sits in the wedge on the face's front: after
                // a forward fin, before a reversed one.
                let mut ring = work.radial_fins(edge);
                let at = ring.iter().position(|&f| f == k).unwrap_or(0);
                let place = if senses[&k] { at + 1 } else { at };
                ring.insert(place, g);
                work.link_ring(&ring);
                made.push(g);
            }
            let n = made.len();
            for (j, &g) in made.iter().enumerate() {
                // made[j] runs k_j backwards: it follows made[j + 1].
                work.half_edges[g].next = made[(j + n - 1) % n];
                work.half_edges[g].prev = made[(j + 1) % n];
            }
            work.loops[new_loop].half_edge = made[0];
        }
        *self = work;
        Ok(twin)
    }

    /// Pair each fin of `gone` with the fin of `keep` running the other
    /// way between the same vertices along a coincident edge.
    fn pair_face_fins(
        &self,
        keep: FaceHandle,
        gone: FaceHandle,
    ) -> Result<Vec<(HalfEdgeHandle, HalfEdgeHandle)>, EulerError> {
        let fins_of = |f: FaceHandle| -> Vec<HalfEdgeHandle> {
            let data = &self.faces[f];
            std::iter::once(data.outer_loop)
                .chain(data.inner_loops.iter().copied())
                .flat_map(|l| self.loop_cycle(l))
                .collect()
        };
        let (kf, gf) = (fins_of(keep), fins_of(gone));
        let mismatch = EulerError::FacesDoNotMatch { keep, gone };
        if kf.len() != gf.len() || kf.is_empty() {
            return Err(mismatch);
        }
        let tolerance = self.tolerance.linear();
        let mut used = BTreeSet::new();
        let mut pairs = Vec::with_capacity(gf.len());
        for &g in &gf {
            let (go, ge) = (
                self.half_edges[g].origin,
                self.fin_end(g).ok_or(mismatch.clone())?,
            );
            let best = kf
                .iter()
                .copied()
                .filter(|k| !used.contains(k))
                .filter(|&k| self.half_edges[k].origin == ge && self.fin_end(k) == Some(go))
                .map(|k| {
                    let (ek, eg) = (self.half_edges[k].edge, self.half_edges[g].edge);
                    let d = if ek == eg {
                        0.0
                    } else {
                        self.curve_distance(eg, ek)
                    };
                    (k, d)
                })
                .filter(|&(_, d)| d <= tolerance)
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .ok_or(mismatch.clone())?;
            used.insert(best.0);
            pairs.push((best.0, g));
        }
        Ok(pairs)
    }

    /// Before fin `fin` leaves `edge`: when it is the edge's first fin
    /// (the one running the stored curve forward), make another fin
    /// leaving the same vertex the first, or with none, the first of the
    /// rest with the curve's parameter window reversed so it still runs
    /// the curve forward.
    fn keep_first_fin_off(&mut self, edge: EdgeHandle, fin: HalfEdgeHandle) {
        if self.edges.get(edge).map(|e| e.half_edge) != Some(fin) {
            return;
        }
        let origin = self.half_edges[fin].origin;
        let rest: Vec<_> = self
            .radial_fins(edge)
            .into_iter()
            .filter(|&f| f != fin)
            .collect();
        if let Some(&f) = rest.iter().find(|&&f| self.half_edges[f].origin == origin) {
            self.edges[edge].half_edge = f;
        } else if let Some(&f) = rest.first() {
            self.edges[edge].half_edge = f;
            if let Some(c) = self.edges[edge].curve.as_mut() {
                std::mem::swap(&mut c.t_start, &mut c.t_end);
            }
        }
    }

    /// How far `gone`'s surface lies from `keep`'s inside the faces: at
    /// the centroid of `keep`'s outer-loop vertices, projected to its
    /// surface (zero where either face has no surface).
    fn surface_gap(&self, keep: FaceHandle, gone: FaceHandle) -> f64 {
        let (Some(sk), Some(sg)) = (
            self.faces[keep].surface.as_ref(),
            self.faces[gone].surface.as_ref(),
        ) else {
            return 0.0;
        };
        let pts: Vec<Point3> = self
            .loop_cycle(self.faces[keep].outer_loop)
            .iter()
            .map(|&f| self.vertices[self.half_edges[f].origin].position)
            .collect();
        if pts.is_empty() {
            return 0.0;
        }
        let n = pts.len() as f64;
        let c = pts.iter().fold(Point3::new(0.0, 0.0, 0.0), |a, p| {
            Point3::new(a.x + p.x / n, a.y + p.y / n, a.z + p.z / n)
        });
        let Some((u, v)) = sk.geometry.closest_point(&c) else {
            return 0.0;
        };
        let on_keep = sk.geometry.point_at(u, v);
        sg.geometry
            .closest_point(&on_keep)
            .map_or(f64::INFINITY, |(a, b)| {
                sg.geometry.point_at(a, b).distance_to(&on_keep)
            })
    }
}
