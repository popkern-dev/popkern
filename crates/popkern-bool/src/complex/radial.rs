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
//! cylinder resting on a plane) are broken by an offset probe: each
//! tied face is followed a short distance `h` into itself and the
//! probe points are compared in the normal plane — the second-order
//! (curvature) separation of tangent faces. Still-tied fins are
//! coincident faces the lift did not merge; that is a typed
//! `Ambiguous` finding, never a guess. The order is computed at three
//! stations and must agree at all of them: a change means faces cross
//! inside the edge, i.e. an intersection curve was missed
//! (`Inconsistent`).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use super::ids::{EdgeId, FinId};
use super::model::CellComplex;

/// Why radial ordering declined.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RadialDecline {
    /// Two or more fins leave the edge in the same direction and the
    /// offset probe cannot separate them.
    Ambiguous {
        /// The edge.
        edge: EdgeId,
        /// Fins that stayed tied.
        fins: Vec<FinId>,
    },
    /// The cyclic order differs between stations along the edge.
    Inconsistent {
        /// The edge.
        edge: EdgeId,
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
    pub operands: [&'a BRep; 2],
}

/// One station: the point, the unit tangent along start → end.
fn station(c: &CellComplex, edge: EdgeId, s: f64) -> Option<(Point3, Vec3)> {
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

/// The fin's face surface (primary owner) at the point nearest `p`:
/// (folded normal, outward multiplier).
fn face_frame(
    c: &CellComplex,
    geo: &FaceGeometry<'_>,
    fin: FinId,
    p: Point3,
) -> Option<(Vec3, f64, popkern_core::geometry::surface::SurfaceOnFace)> {
    let f = c.fins.get(fin.index())?;
    let face = c.faces.get(f.face.index())?;
    let owner = face.owners.first()?;
    let surf = geo
        .operands
        .get(owner.operand)?
        .faces
        .get(owner.face)?
        .surface
        .clone()?;
    let (u, v) = surf.geometry.closest_point(&p)?;
    let n = surf.folded_normal(u, v)?;
    let out = if face.folded_is_outward { 1.0 } else { -1.0 };
    Some((n, out, surf))
}

/// Unit direction from the edge into the fin's face, in the plane
/// normal to `t`, and the face's outward normal there.
fn fin_direction(
    c: &CellComplex,
    geo: &FaceGeometry<'_>,
    fin: FinId,
    p: Point3,
    t: Vec3,
) -> Option<(Vec3, Vec3, popkern_core::geometry::surface::SurfaceOnFace)> {
    let f = c.fins.get(fin.index())?;
    let (n, out, surf) = face_frame(c, geo, fin, p)?;
    let tf = if f.forward { t } else { t * -1.0 };
    let d = n.cross(&tf);
    let d = (d - t * d.dot(&t)).normalized()?;
    Some((d, n * out, surf))
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

/// Order the fins of one edge at station `s`; returns the order and,
/// per fin in that order, whether its Ahead side faces the next wedge.
fn order_at(
    c: &CellComplex,
    geo: &FaceGeometry<'_>,
    edge: EdgeId,
    s: f64,
    tie_angle: f64,
    probe: f64,
) -> Result<Vec<(FinId, bool)>, RadialDecline> {
    let fins = c.edges[edge.index()].radial.clone();
    let (p, t) = station(c, edge, s).ok_or(RadialDecline::Unevaluable { edge, fin: None })?;
    let mut entries = Vec::with_capacity(fins.len());
    for &fin in &fins {
        let (d, n_out, surf) =
            fin_direction(c, geo, fin, p, t).ok_or(RadialDecline::Unevaluable {
                edge,
                fin: Some(fin),
            })?;
        entries.push((fin, d, n_out, surf));
    }
    let reference = entries[0].1;
    let mut keyed: Vec<(f64, f64, FinId, bool)> = entries
        .iter()
        .map(|(fin, d, n_out, _)| {
            let ahead_next = n_out.dot(&t.cross(d)) > 0.0;
            (angle(reference, *d, t), 0.0, *fin, ahead_next)
        })
        .collect();
    // Tie ladder: second key = probe-point angle for fins whose
    // first-order directions agree.
    let n = keyed.len();
    for i in 0..n {
        let tied = (0..n).any(|j| {
            j != i && {
                let da = (keyed[i].0 - keyed[j].0).abs();
                da.min(std::f64::consts::TAU - da) < tie_angle
            }
        });
        if tied {
            let (_, d, _, surf) = &entries[i];
            let q0 = p + *d * probe;
            let q = surf
                .geometry
                .closest_point(&q0)
                .map(|(u, v)| surf.geometry.point_at(u, v))
                .ok_or(RadialDecline::Unevaluable {
                    edge,
                    fin: Some(entries[i].0),
                })?;
            let w = q - p;
            let w = w - t * w.dot(&t);
            keyed[i].1 = match w.normalized() {
                Some(w) => {
                    let a = angle(reference, w, t);
                    // Signed small deviation from the first-order ray.
                    let dev = a - keyed[i].0;
                    if dev > std::f64::consts::PI {
                        dev - std::f64::consts::TAU
                    } else if dev < -std::f64::consts::PI {
                        dev + std::f64::consts::TAU
                    } else {
                        dev
                    }
                }
                None => 0.0,
            };
        }
    }
    // A direction just below 2π is tied with the reference at 0: wrap
    // it below zero so the tie forms one cluster.
    for k in &mut keyed {
        if k.0 > std::f64::consts::TAU - tie_angle {
            k.0 -= std::f64::consts::TAU;
        }
    }
    // Cluster first-order ties, then order inside each cluster by the
    // probe deviation.
    keyed.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut out: Vec<(FinId, bool)> = Vec::with_capacity(n);
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && keyed[j].0 - keyed[i].0 < tie_angle {
            j += 1;
        }
        let mut cluster = keyed[i..j].to_vec();
        cluster.sort_by(|a, b| a.1.total_cmp(&b.1));
        for w in cluster.windows(2) {
            if (w[1].1 - w[0].1).abs() < tie_angle * 1e-3 {
                return Err(RadialDecline::Ambiguous {
                    edge,
                    fins: cluster.iter().map(|k| k.2).collect(),
                });
            }
        }
        out.extend(cluster.iter().map(|k| (k.2, k.3)));
        i = j;
    }
    Ok(out)
}

/// Length of an edge, sampled (closed edges have no chord).
fn edge_length(c: &CellComplex, edge: EdgeId) -> f64 {
    let pts: Vec<Point3> = (0..=8)
        .filter_map(|k| station(c, edge, f64::from(k) / 8.0).map(|(p, _)| p))
        .collect();
    pts.windows(2).map(|w| w[0].distance_to(&w[1])).sum()
}

/// Whether two cyclic orders are rotations of each other.
fn same_cycle(a: &[(FinId, bool)], b: &[(FinId, bool)]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let Some(start) = b.iter().position(|x| x.0 == a[0].0) else {
        return false;
    };
    (0..a.len()).all(|k| a[k] == b[(start + k) % b.len()])
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
    let mut sides = RadialSides {
        ahead_faces_next: vec![false; c.fins.len()],
    };
    for ei in 0..c.edges.len() {
        let edge = EdgeId::from_index(ei);
        let m = c.edges[ei].radial.len();
        if m == 0 {
            continue;
        }
        let probe = (edge_length(c, edge) * 1.0e-2).max(linear * 100.0);
        let stations = if m >= 3 {
            &[0.5, 0.25, 0.75][..]
        } else {
            &[0.5][..]
        };
        let first = order_at(c, geo, edge, stations[0], tie_angle, probe)?;
        for &s in &stations[1..] {
            let other = order_at(c, geo, edge, s, tie_angle, probe)?;
            if !same_cycle(&first, &other) {
                return Err(RadialDecline::Inconsistent { edge });
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
