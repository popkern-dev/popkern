//! Radial order of the fins around every complex edge (PK-973).
//!
//! At a station on the edge, each fin's face leaves the edge in the
//! direction `d = folded_normal × fin_tangent` (the face interior lies
//! left of the fin seen from the folded normal). Projected onto the
//! plane normal to the edge tangent `t`, the directions sort by angle
//! about `t`; the wedge between two consecutive fins is one region's
//! slice of space near the edge. Each fin also records which of its
//! face's sides faces the NEXT wedge (counter-clockwise about `t`):
//! the Ahead side does iff the face's outward normal has a positive
//! component along `t × d`.
//!
//! Ties (faces leaving in the same direction: tangent contact, a
//! cylinder resting on a plane) go down the tie ladder
//! (`radial_ladder`): cross-section curvature, then an offset probe;
//! still-tied fins are coincident faces the join did not merge, a
//! typed `Ambiguous` finding with the measured keys, never a guess.
//!
//! The order is computed at several stations (three, five on a curved
//! edge) and must agree at all of them: a change means faces cross
//! inside the edge, i.e. an intersection curve was missed
//! (`Inconsistent`, naming the disagreeing stations).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use super::ids::{EdgeId, FinId};
use super::model::CellComplex;
use super::radial_ladder::{TieEvidence, TieLadder, cross_section_curvature, runs};

/// Why radial ordering declined.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RadialDecline {
    /// Two or more fins leave the edge in the same direction and no
    /// rung of the tie ladder separates them.
    Ambiguous {
        /// The edge.
        edge: EdgeId,
        /// Fins that stayed tied.
        fins: Vec<FinId>,
        /// The station (fraction of the edge from its start).
        station: f64,
        /// What the ladder measured for every fin of the tie.
        evidence: Vec<TieEvidence>,
    },
    /// The cyclic order differs between stations along the edge: faces
    /// cross inside it, an intersection curve the intersection stage
    /// missed.
    Inconsistent {
        /// The edge.
        edge: EdgeId,
        /// The first station's order and the first that disagrees.
        stations: Vec<(f64, Vec<FinId>)>,
    },
    /// A face's normal or the edge tangent could not be evaluated.
    Unevaluable {
        /// The edge.
        edge: EdgeId,
        /// The fin whose face failed, if the edge itself did not.
        fin: Option<FinId>,
    },
}

/// Per-fin result of the radial stage.
#[derive(Debug, Clone, Default)]
pub(crate) struct RadialSides {
    /// For each fin: whether its face's Ahead side faces the next
    /// (counter-clockwise) wedge about the edge tangent.
    pub ahead_faces_next: Vec<bool>,
}

/// Geometry the radial stage needs for a fin's face.
pub(crate) struct FaceGeometry<'a> {
    /// The operand BReps.
    pub operands: &'a [&'a BRep],
    /// The intersection graph the complex lifts, when there is one:
    /// its half-edges' pcurves seed each fin's projection onto its face.
    pub graph: Option<&'a crate::merge::graph::IntersectionGraph>,
}

/// One station: the point, the unit tangent along start → end.
pub(super) fn station(c: &CellComplex, edge: EdgeId, s: f64) -> Option<(Point3, Vec3)> {
    let e = c.edges.get(edge.index())?;
    let ps = c.vertices.get(e.start.index())?.position;
    let pe = c.vertices.get(e.end.index())?.position;
    match &e.curve {
        Some(curve) => {
            let (t0, t1) = (curve.t_start, curve.t_end);
            let a = curve.geometry.point_at(t0);
            let b = curve.geometry.point_at(t1);
            // Which parameter end is the edge's start (closed edges:
            // both; keep the curve's own direction).
            let reversed = e.start != e.end && a.distance_to(&pe) < a.distance_to(&ps);
            let s = if reversed { 1.0 - s } else { s };
            let t = t0 + (t1 - t0) * s;
            let p = curve.geometry.point_at(t);
            let mut d = curve.geometry.derivative_at(t, 1);
            if (t1 - t0) < 0.0 {
                d = d * -1.0;
            }
            if reversed {
                d = d * -1.0;
            }
            let _ = b; // the far end is implied by `reversed`
            Some((p, d.normalized()?))
        }
        None => {
            let d = (pe - ps).normalized()?;
            Some((ps + (pe - ps) * s, d))
        }
    }
}

/// A fin's face at a station.
struct FinFrame<'a> {
    /// Unit direction from the edge into the face, normal to `t`.
    d: Vec3,
    /// The face's surface (primary owner).
    surf: &'a SurfaceOnFace,
    /// The station's parameters on it.
    uv: (f64, f64),
}

/// The fin's face at the point nearest `p`.
fn fin_frame<'a>(
    c: &CellComplex,
    geo: &FaceGeometry<'a>,
    fin: FinId,
    p: Point3,
    t: Vec3,
) -> Option<FinFrame<'a>> {
    let f = c.fins.get(fin.index())?;
    let face = c.faces.get(f.face.index())?;
    let owner = face.owners.first()?;
    let surf = geo
        .operands
        .get(owner.operand)?
        .faces
        .get(owner.face)?
        .surface
        .as_ref()?;
    let reach = c.edges.get(f.edge.index())?.tolerance.max(1e-9) * 1e3;
    let uv = seeded_uv(geo, f.source_half_edge, surf, p, reach)
        .or_else(|| surf.geometry.closest_point(&p))?;
    let n = surf.folded_normal(uv.0, uv.1)?;
    let tf = if f.forward { t } else { t * -1.0 };
    let d = n.cross(&tf);
    let d = (d - t * d.dot(&t)).normalized()?;
    Some(FinFrame { d, surf, uv })
}

/// Pcurve stations a fin's projection seed is chosen from.
const SEED_STATIONS: usize = 16;

/// The station's parameters on `surf` by single-start projection from
/// the fin's graph pcurve station whose image lies nearest `p`, kept
/// when its foot lies within `reach`. A multi-start projection onto a
/// multi-turn NURBS wall costs in proportion to its turns, once per fin
/// per station (PK-1039: 0.58 s / 2.12 s of a 40 / 80 thread cut, now
/// 0.05 s / 0.08 s). The pcurve's parameter need not track the edge
/// curve's (record uv fits drift: PK-1031), hence the nearest station
/// rather than the one at the station's fraction; and either half of
/// the graph edge may carry this face's pcurve (the other carries its
/// partner face's, whose uv means nothing here: seeds 1e6 away).
fn seeded_uv(
    geo: &FaceGeometry<'_>,
    half_edge: Option<usize>,
    surf: &SurfaceOnFace,
    p: Point3,
    reach: f64,
) -> Option<(f64, f64)> {
    let graph = geo.graph?;
    let he = graph.half_edges.get(half_edge?)?;
    let off = |q: (f64, f64)| surf.geometry.point_at(q.0, q.1).distance_to(&p);
    let (u0, v0) = [Some(he), graph.half_edges.get(he.twin)]
        .into_iter()
        .flatten()
        .filter_map(|h| h.pcurve.as_ref())
        .filter(|pc| {
            let (a, b) = pc.parameter_range();
            a.is_finite() && b.is_finite()
        })
        .flat_map(|pc| {
            let (a, b) = pc.parameter_range();
            (0..=SEED_STATIONS)
                .map(move |k| pc.point_at(a + (b - a) * k as f64 / SEED_STATIONS as f64))
        })
        .map(|q| (q.x, q.y))
        .min_by(|x, y| off(*x).total_cmp(&off(*y)))?;
    popkern_core::geometry::uv_projection::project_point_to_uv_near(&surf.geometry, &p, u0, v0)
        .filter(|&q| off(q) <= reach)
}

/// Whether `fin`'s face's Ahead side faces the next wedge about its
/// edge. With `d = n × tf` normal to `t`, `t × d` is `+n` when the fin
/// runs along `t` and `−n` against it, so the outward normal points
/// into the next wedge exactly when the face's folds outward and the
/// fin runs forward, or neither: topology alone, no evaluation.
fn ahead_faces_next(c: &CellComplex, fin: FinId) -> bool {
    let f = &c.fins[fin.index()];
    c.faces[f.face.index()].folded_is_outward == f.forward
}

/// Angle of `d` about `t` from `reference`, in [0, 2π).
fn angle(reference: Vec3, d: Vec3, t: Vec3) -> f64 {
    let a = reference.cross(&d).dot(&t).patan2(reference.dot(&d));
    if a < 0.0 {
        a + std::f64::consts::TAU
    } else {
        a
    }
}

/// The scales one edge is ordered at.
#[derive(Debug, Clone, Copy)]
struct EdgeScale {
    /// Angular tolerance of a first-order tie.
    tie_angle: f64,
    /// Offset probe distance.
    probe: f64,
    /// The edge's feature length: curvatures closer than
    /// `tie_angle·1e-3 / feature` turn the faces apart by less than a
    /// thousandth of the angular tolerance across it.
    feature: f64,
    /// The rungs in use.
    ladder: TieLadder,
    /// The shortest probe the ladder takes (a few linear tolerances).
    floor: f64,
}

/// One fin's sort keys at a station.
#[derive(Debug, Clone, Copy)]
struct Key {
    /// Index into the station's frames.
    index: usize,
    fin: FinId,
    angle: f64,
    ahead_next: bool,
    curvature: Option<f64>,
    probe: Option<f64>,
}

impl Key {
    fn evidence(&self) -> TieEvidence {
        TieEvidence {
            fin: self.fin,
            angle: self.angle,
            curvature: self.curvature,
            probe: self.probe,
        }
    }
}

/// Order a first-order tie down the ladder (see `radial_ladder`).
/// The signed angle, in (-π, π], of `f`'s face a distance `dist` along
/// its first-order ray from `p`, measured from that ray about `t`.
fn probe_angle(f: &FinFrame<'_>, (p, t): (Point3, Vec3), dist: f64) -> Option<f64> {
    let target = p + f.d * dist;
    let q = popkern_core::geometry::uv_projection::project_point_to_uv_near(
        &f.surf.geometry,
        &target,
        f.uv.0,
        f.uv.1,
    )
    .or_else(|| f.surf.geometry.closest_point(&target))
    .map(|(u, v)| f.surf.geometry.point_at(u, v))?;
    let w = q - p;
    let w = w - t * w.dot(&t);
    Some(
        w.normalized()
            .map_or(0.0, |w| f.d.cross(&w).dot(&t).patan2(f.d.dot(&w))),
    )
}

fn break_tie(
    c: &CellComplex,
    mut cluster: Vec<Key>,
    frames: &[FinFrame<'_>],
    (p, t): (Point3, Vec3),
    (edge, station): (EdgeId, f64),
    scale: &EdgeScale,
) -> Result<Vec<Key>, RadialDecline> {
    let unevaluable = |fin| RadialDecline::Unevaluable {
        edge,
        fin: Some(fin),
    };
    // The faces' own order first: probed within the tied faces (half
    // the narrowest one's reach, at most the edge's probe), the probe
    // sees first and second order together. Curvature alone orders
    // faces as they lie beyond a crossing: a cylinder crossing a plane
    // at 5.9e-4 rad turns into the 1e-5 lens between them at first
    // order and away from it at second, and curvature put it on the far
    // side, joining the lens's two sides. For faces truly tangent the
    // probe reads the curvature term and agrees with that rung, which
    // stays for ties the probe cannot resolve.
    let reach = cluster
        .iter()
        .map(|k| super::radial_stations::face_clearance(c, k.fin, p))
        .fold(f64::INFINITY, f64::min);
    let near = (0.5 * reach).max(scale.floor).min(scale.probe);
    if scale.ladder.offset_probe {
        let mut narrow = cluster.clone();
        for k in &mut narrow {
            k.probe = Some(
                probe_angle(&frames[k.index], (p, t), near).ok_or_else(|| unevaluable(k.fin))?,
            );
        }
        narrow.sort_by(|a, b| a.probe.unwrap_or(0.0).total_cmp(&b.probe.unwrap_or(0.0)));
        let resolved = narrow.windows(2).all(|w| {
            w[1].probe.unwrap_or(0.0) - w[0].probe.unwrap_or(0.0) >= scale.tie_angle * 1e-3
        });
        if resolved {
            return Ok(narrow);
        }
    }
    let mut groups = if scale.ladder.curvature {
        for k in &mut cluster {
            let f = &frames[k.index];
            let kappa = cross_section_curvature(&f.surf.geometry, f.uv, f.d, t)
                .ok_or_else(|| unevaluable(k.fin))?;
            k.curvature = Some(kappa);
        }
        cluster.sort_by(|a, b| {
            a.curvature
                .unwrap_or(0.0)
                .total_cmp(&b.curvature.unwrap_or(0.0))
        });
        runs(
            cluster,
            |k| k.curvature.unwrap_or(0.0),
            scale.tie_angle * 1e-3 / scale.feature,
        )
    } else {
        vec![cluster]
    };
    let mut tied = false;
    for g in groups.iter_mut().filter(|g| g.len() > 1) {
        if !scale.ladder.offset_probe {
            tied = true;
            continue;
        }
        for k in g.iter_mut() {
            k.probe = Some(
                probe_angle(&frames[k.index], (p, t), scale.probe)
                    .ok_or_else(|| unevaluable(k.fin))?,
            );
        }
        g.sort_by(|a, b| a.probe.unwrap_or(0.0).total_cmp(&b.probe.unwrap_or(0.0)));
        tied |= g.windows(2).any(|w| {
            w[1].probe.unwrap_or(0.0) - w[0].probe.unwrap_or(0.0) < scale.tie_angle * 1e-3
        });
    }
    let ordered: Vec<Key> = groups.into_iter().flatten().collect();
    if tied {
        return Err(RadialDecline::Ambiguous {
            edge,
            fins: ordered.iter().map(|k| k.fin).collect(),
            station,
            evidence: ordered.iter().map(Key::evidence).collect(),
        });
    }
    Ok(ordered)
}

/// Order the fins of one edge at station `s`; returns the order and,
/// per fin in that order, whether its Ahead side faces the next wedge.
fn order_at(
    c: &CellComplex,
    geo: &FaceGeometry<'_>,
    edge: EdgeId,
    s: f64,
    scale: &EdgeScale,
) -> Result<Vec<(FinId, bool)>, RadialDecline> {
    let fins = &c.edges[edge.index()].radial;
    let (p, t) = station(c, edge, s).ok_or(RadialDecline::Unevaluable { edge, fin: None })?;
    let mut frames = Vec::with_capacity(fins.len());
    for &fin in fins {
        frames.push(
            fin_frame(c, geo, fin, p, t).ok_or(RadialDecline::Unevaluable {
                edge,
                fin: Some(fin),
            })?,
        );
    }
    let reference = frames[0].d;
    let mut keyed: Vec<Key> = fins
        .iter()
        .zip(&frames)
        .enumerate()
        .map(|(index, (&fin, f))| {
            let mut a = angle(reference, f.d, t);
            // A direction just below 2π is tied with the reference at
            // 0: wrap it below zero so the tie forms one cluster.
            if a > std::f64::consts::TAU - scale.tie_angle {
                a -= std::f64::consts::TAU;
            }
            Key {
                index,
                fin,
                angle: a,
                ahead_next: ahead_faces_next(c, fin),
                curvature: None,
                probe: None,
            }
        })
        .collect();
    keyed.sort_by(|a, b| a.angle.total_cmp(&b.angle));
    let mut out = Vec::with_capacity(keyed.len());
    for cluster in runs(keyed, |k| k.angle, scale.tie_angle) {
        let cluster = if cluster.len() > 1 {
            break_tie(c, cluster, &frames, (p, t), (edge, s), scale)?
        } else {
            cluster
        };
        out.extend(cluster.iter().map(|k| (k.fin, k.ahead_next)));
    }
    Ok(out)
}

/// Order the fins of every edge with three or more fins, and record
/// every fin's side-facing (two-fin edges need it too). `tie_angle` is
/// the angular tolerance for first-order ties; the tie-ladder probe
/// distance is relative to each edge's length, floored at a hundred
/// linear tolerances — never an absolute constant.
pub(crate) fn order_radially(
    c: &mut CellComplex,
    geo: &FaceGeometry<'_>,
    tie_angle: f64,
    linear: f64,
) -> Result<RadialSides, RadialDecline> {
    order_radially_with(c, geo, tie_angle, linear, TieLadder::default())
}

/// `order_radially` with the given ladder rungs.
pub(crate) fn order_radially_with(
    c: &mut CellComplex,
    geo: &FaceGeometry<'_>,
    tie_angle: f64,
    linear: f64,
    ladder: TieLadder,
) -> Result<RadialSides, RadialDecline> {
    let mut sides = RadialSides {
        ahead_faces_next: vec![false; c.fins.len()],
    };
    for ei in 0..c.edges.len() {
        let edge = EdgeId::from_index(ei);
        let m = c.edges[ei].radial.len();
        if m == 0 {
            continue;
        }
        // One or two fins: any cyclic order is the order, and the sides
        // need no geometry. Evaluating every fin's face cost a
        // projection onto it per edge (PK-1039: on a multi-turn thread
        // cut, 0.58 s / 2.12 s at 205 / 365 edges, the complex's whole
        // cost).
        if m < 3 {
            for &fin in &c.edges[ei].radial {
                sides.ahead_faces_next[fin.index()] = ahead_faces_next(c, fin);
            }
            c.edges[ei].radial_ordered = true;
            continue;
        }
        let length = super::radial_stations::edge_length(c, edge);
        let scale = EdgeScale {
            tie_angle,
            probe: (length * 1.0e-2).max(linear * 100.0),
            feature: length.max(linear * 100.0),
            ladder,
            floor: linear * 10.0,
        };
        let stations = super::radial_stations::stations(c, geo, edge, length, linear);
        let first = order_at(c, geo, edge, stations[0], &scale)?;
        for &s in &stations[1..] {
            let other = order_at(c, geo, edge, s, &scale)?;
            if !super::radial_stations::same_cycle(&first, &other) {
                let ids = |o: &[(FinId, bool)]| o.iter().map(|x| x.0).collect();
                return Err(RadialDecline::Inconsistent {
                    edge,
                    stations: vec![(stations[0], ids(&first)), (s, ids(&other))],
                });
            }
        }
        c.edges[ei].radial = first.iter().map(|x| x.0).collect();
        c.edges[ei].radial_ordered = true;
        for (fin, ahead_next) in first {
            sides.ahead_faces_next[fin.index()] = ahead_next;
        }
    }
    Ok(sides)
}

#[cfg(test)]
#[path = "radial_tests.rs"]
mod tests;
