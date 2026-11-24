//! Radial order of the fins about an edge.
//!
//! At a station on the edge, each fin's face leaves the edge in the
//! direction `d = folded_normal × fin_tangent` (the face interior lies
//! left of the fin seen from the folded normal). Projected onto the
//! plane normal to the edge tangent `t`, the directions sort by angle
//! about `t` (right-handed): the wedge between two consecutive fins is
//! one region's slice of space near the edge.
//!
//! Ties (faces leaving in the same direction: tangent contact, a
//! cylinder resting on a plane) go down the tie ladder ([`ladder`]):
//! cross-section curvature, then an offset probe; still-tied fins are
//! coincident faces, a typed [`RadialOrderDecline::Ambiguous`] with the
//! measured keys, never a guess.
//!
//! The order is taken at several stations (three, five on a curved
//! edge, one on a straight edge between planes) and must agree at all
//! of them: a change means faces cross inside the edge, an intersection
//! the edge does not record ([`RadialOrderDecline::Inconsistent`]).
//!
//! The engine is generic over where the fins come from
//! ([`RadialFins`]): a body's edge ([`BRep::radial_order`]) or the
//! boolean cell complex's.
//!
//! [`BRep::radial_order`]: crate::topology::BRep::radial_order

#[cfg(test)]
pub(crate) mod book_test_helpers;
mod brep_order;
pub mod ladder;

pub use brep_order::RadialOrderError;
pub use ladder::{TieEvidence, TieLadder, cross_section_curvature, runs, tangent_gap};

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// One fin's face at a station.
#[derive(Debug, Clone, Copy)]
pub struct RadialFrame<'a> {
    /// Unit direction from the edge into the face, normal to `t`.
    pub d: Vec3,
    /// The face's surface.
    pub surface: &'a SurfaceGeometry,
    /// The station's parameters on it.
    pub uv: (f64, f64),
}

/// The fins about one edge, as the engine reads them.
pub trait RadialFins {
    /// How many fins the edge has.
    fn fin_count(&self) -> usize;
    /// The point and unit tangent at fraction `s` of the edge, the
    /// tangent along the direction the order is taken about.
    fn station(&self, s: f64) -> Option<(Point3, Vec3)>;
    /// Fin `fin`'s face at station `(p, t)`.
    fn frame(&self, fin: usize, p: Point3, t: Vec3) -> Option<RadialFrame<'_>>;
    /// How far fin `fin`'s face reaches from `p` before another of its
    /// edges (infinite when unknown): tied faces are probed within it.
    fn clearance(&self, fin: usize, p: Point3) -> f64;
    /// The distance between the edge's end vertices.
    fn chord(&self) -> f64;
    /// Whether every fin's face lies on a plane.
    fn all_planar(&self) -> bool;
}

/// The tolerances an order is taken at.
#[derive(Debug, Clone, Copy)]
pub struct RadialTolerance {
    /// Angular tolerance of a first-order tie (radians).
    pub tie_angle: f64,
    /// The operation's linear tolerance: the probe distance is the
    /// edge's length × 1e-2 floored at a hundred of these, never an
    /// absolute constant.
    pub linear: f64,
    /// The ladder rungs in use.
    pub ladder: TieLadder,
}

/// Why radial ordering declined. Fins are indices into the source.
#[derive(Debug, Clone, PartialEq)]
pub enum RadialOrderDecline {
    /// Two or more fins leave the edge in the same direction and no
    /// rung of the tie ladder separates them.
    Ambiguous {
        /// Fins that stayed tied.
        fins: Vec<usize>,
        /// The station (fraction of the edge from its start).
        station: f64,
        /// What the ladder measured for every fin of the tie.
        evidence: Vec<TieEvidence>,
    },
    /// The cyclic order differs between stations along the edge.
    Inconsistent {
        /// The first station's order and the first that disagrees.
        stations: Vec<(f64, Vec<usize>)>,
    },
    /// A face's normal or the edge tangent could not be evaluated.
    Unevaluable {
        /// The fin whose face failed, if the edge itself did not.
        fin: Option<usize>,
    },
}

/// The scales one edge is ordered at.
#[derive(Debug, Clone, Copy)]
struct EdgeScale {
    tie_angle: f64,
    probe: f64,
    /// Curvatures closer than `tie_angle·1e-3 / feature` turn the faces
    /// apart by less than a thousandth of the angular tolerance across
    /// the edge's feature length.
    feature: f64,
    ladder: TieLadder,
    /// The shortest probe the ladder takes (a few linear tolerances).
    floor: f64,
    /// The operation's linear tolerance.
    linear: f64,
}

/// One fin's sort keys at a station.
#[derive(Debug, Clone, Copy)]
struct Key {
    fin: usize,
    angle: f64,
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

/// The fins of `src` in radial order (right-handed about the station
/// tangent), starting from fin 0, agreeing at every station.
///
/// # Errors
///
/// [`RadialOrderDecline`] when a tie survives the ladder, the order
/// changes along the edge, or a face cannot be evaluated.
pub fn order_fins(
    src: &impl RadialFins,
    tolerance: RadialTolerance,
) -> Result<Vec<usize>, RadialOrderDecline> {
    let m = src.fin_count();
    if m < 3 {
        // Any cyclic order of one or two fins is the order.
        return Ok((0..m).collect());
    }
    let length = edge_length(src);
    let scale = EdgeScale {
        tie_angle: tolerance.tie_angle,
        probe: (length * 1.0e-2).max(tolerance.linear * 100.0),
        feature: length.max(tolerance.linear * 100.0),
        ladder: tolerance.ladder,
        floor: tolerance.linear * 10.0,
        linear: tolerance.linear,
    };
    let stations = stations(src, length, tolerance.linear);
    let first = order_at(src, stations[0], &scale)?;
    for &s in &stations[1..] {
        let other = order_at(src, s, &scale)?;
        if !same_cycle(&first, &other) {
            return Err(RadialOrderDecline::Inconsistent {
                stations: vec![(stations[0], first), (s, other)],
            });
        }
    }
    let start = first.iter().position(|&f| f == 0).unwrap_or(0);
    Ok((0..m).map(|k| first[(start + k) % m]).collect())
}

/// The stations (edge fractions) the order is taken at, the first
/// deciding it. A straight edge between planes gets one: each fin's
/// inward direction is then the same along the whole edge.
fn stations(src: &impl RadialFins, length: f64, linear: f64) -> &'static [f64] {
    if length > src.chord() + linear {
        &[0.5, 0.25, 0.75, 0.125, 0.875]
    } else if src.all_planar() {
        &[0.5]
    } else {
        &[0.5, 0.25, 0.75]
    }
}

/// Length of the edge, sampled (closed edges have no chord).
fn edge_length(src: &impl RadialFins) -> f64 {
    let pts: Vec<Point3> = (0..=8)
        .filter_map(|k| src.station(f64::from(k) / 8.0).map(|(p, _)| p))
        .collect();
    pts.windows(2).map(|w| w[0].distance_to(&w[1])).sum()
}

/// Whether two cyclic orders are rotations of each other.
pub fn same_cycle(a: &[usize], b: &[usize]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let Some(&a0) = a.first() else {
        return true;
    };
    let Some(start) = b.iter().position(|&x| x == a0) else {
        return false;
    };
    (0..a.len()).all(|k| a[k] == b[(start + k) % b.len()])
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

/// The signed angle, in (-π, π], of `f`'s face a distance `dist` along
/// its first-order ray from `p`, measured from that ray about `t`.
fn probe_angle(f: &RadialFrame<'_>, (p, t): (Point3, Vec3), dist: f64) -> Option<f64> {
    let target = p + f.d * dist;
    let q = crate::geometry::uv_projection::project_point_to_uv_near(
        f.surface, &target, f.uv.0, f.uv.1,
    )
    .or_else(|| f.surface.closest_point(&target))
    .map(|(u, v)| f.surface.point_at(u, v))?;
    let w = q - p;
    let w = w - t * w.dot(&t);
    Some(
        w.normalized()
            .map_or(0.0, |w| f.d.cross(&w).dot(&t).patan2(f.d.dot(&w))),
    )
}

/// Order a first-order tie down the ladder.
fn break_tie(
    src: &impl RadialFins,
    mut cluster: Vec<Key>,
    frames: &[RadialFrame<'_>],
    (p, t): (Point3, Vec3),
    (station, beyond): (f64, f64),
    scale: &EdgeScale,
) -> Result<Vec<Key>, RadialOrderDecline> {
    let unevaluable = |fin| RadialOrderDecline::Unevaluable { fin: Some(fin) };
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
        .map(|k| src.clearance(k.fin, p))
        .fold(f64::INFINITY, f64::min);
    let near = (0.5 * reach).max(scale.floor).min(scale.probe).max(beyond);
    if scale.ladder.offset_probe {
        let mut narrow = cluster.clone();
        for k in &mut narrow {
            k.probe =
                Some(probe_angle(&frames[k.fin], (p, t), near).ok_or_else(|| unevaluable(k.fin))?);
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
            let f = &frames[k.fin];
            let kappa = cross_section_curvature(f.surface, f.uv, f.d, t)
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
                probe_angle(&frames[k.fin], (p, t), scale.probe)
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
        return Err(RadialOrderDecline::Ambiguous {
            fins: ordered.iter().map(|k| k.fin).collect(),
            station,
            evidence: ordered.iter().map(Key::evidence).collect(),
        });
    }
    Ok(ordered)
}

/// The fins in order about the tangent at station `s`.
fn order_at(
    src: &impl RadialFins,
    s: f64,
    scale: &EdgeScale,
) -> Result<Vec<usize>, RadialOrderDecline> {
    let (p, t) = src
        .station(s)
        .ok_or(RadialOrderDecline::Unevaluable { fin: None })?;
    let frames = (0..src.fin_count())
        .map(|fin| {
            src.frame(fin, p, t)
                .ok_or(RadialOrderDecline::Unevaluable { fin: Some(fin) })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reference = frames[0].d;
    let mut keyed: Vec<Key> = frames
        .iter()
        .enumerate()
        .map(|(fin, f)| {
            let mut a = angle(reference, f.d, t);
            // A direction just below 2π is tied with the reference at
            // 0: wrap it below zero so the tie forms one cluster.
            if a > std::f64::consts::TAU - scale.tie_angle {
                a -= std::f64::consts::TAU;
            }
            Key {
                fin,
                angle: a,
                curvature: None,
                probe: None,
            }
        })
        .collect();
    keyed.sort_by(|a, b| a.angle.total_cmp(&b.angle));
    let kappa: Vec<Option<f64>> = frames
        .iter()
        .map(|f| ladder::cross_section_curvature(f.surface, f.uv, f.d, t))
        .collect();
    let mut out = Vec::with_capacity(keyed.len());
    for (cluster, beyond) in tie_runs(keyed, &kappa, scale) {
        let cluster = if cluster.len() > 1 {
            break_tie(src, cluster, &frames, (p, t), (s, beyond), scale)?
        } else {
            cluster
        };
        out.extend(cluster.iter().map(|k| k.fin));
    }
    Ok(out)
}

/// The first-order ties among `keyed` (sorted by angle), each with the
/// distance its probe must reach. Neighbours `θ` apart are tied where
/// the gap they open within the edge's feature length is within the
/// linear tolerance ([`tangent_gap`]; a bare angular tie ordered two
/// balls of radius 1000 crossing at 9e-5 rad by curvature, as if they
/// touched, and their lens two tolerances thick read as in neither,
/// SS-E 695): where their cross-section curvatures, `Δκ` apart, bring
/// them back across each other, the lens between them is at most
/// `θ² / 2Δκ` thick, and where that is within the tolerance the faces
/// are tangent at the model's resolution (two cylinders on a ruling overlapping by 5e-7
/// cross at 1.4e-3 rad and recross 1.3e-3 away; ordered by the crossing,
/// the wedge between them lay in both solids). A tie is probed beyond
/// the lens its faces can enclose, twice its width with the faces up to
/// the tolerance apart at the edge (the glued edge of those cylinders sat
/// on the overlap's midline, the faces parallel there and 5e-7 apart:
/// probed 1e-5 out, they read the overlap's order), at most the edge's
/// length. A real crossing's lens (5.9e-4 rad into a 1e-5 lens) stays
/// untied. The run wrapping past 2π joins the first.
fn tie_runs(keyed: Vec<Key>, kappa: &[Option<f64>], scale: &EdgeScale) -> Vec<(Vec<Key>, f64)> {
    // The lens's width when the faces lie up to the tolerance apart at
    // the edge: they separate as `θ r + Δκ r² / 2` against the gap.
    let width = |gap: f64, dk: f64| (gap + (gap * gap + 2.0 * dk * scale.linear).sqrt()) / dk;
    let lens = |a: &Key, b: &Key, gap: f64| -> Option<f64> {
        let dk = match (kappa[a.fin], kappa[b.fin]) {
            (Some(ka), Some(kb)) => (ka - kb).abs(),
            _ => 0.0,
        };
        let beyond = if dk > 0.0 {
            (2.0 * width(gap.max(0.0), dk)).min(scale.feature)
        } else {
            0.0
        };
        // A tangent touch by distance: inside the open tangency band, so a
        // lens exactly a tolerance thick (the near-tangent cap casting's
        // boss 1e-6 off the rim) is a crossing.
        let band = crate::tolerance::Tolerance::default()
            .with_linear(scale.linear)
            .tangency_band(0.0);
        (tangent_gap(gap, dk, scale.feature) < band).then_some(beyond)
    };
    let mut out: Vec<(Vec<Key>, f64)> = Vec::new();
    for k in keyed {
        let joined = out.last_mut().and_then(|(run, beyond)| {
            let last = run.last()?;
            let reach = lens(last, &k, k.angle - last.angle)?;
            *beyond = beyond.max(reach);
            Some(run)
        });
        match joined {
            Some(run) => run.push(k),
            None => out.push((vec![k], 0.0)),
        }
    }
    if out.len() > 1
        && let (Some(first), Some(last)) = (
            out.first().and_then(|r| r.0.first()).copied(),
            out.last().and_then(|r| r.0.last()).copied(),
        )
        && let Some(reach) = lens(
            &last,
            &first,
            first.angle + std::f64::consts::TAU - last.angle,
        )
        && let Some((mut tail, beyond)) = out.pop()
    {
        for k in &mut tail {
            k.angle -= std::f64::consts::TAU;
        }
        let head = &mut out[0];
        tail.append(&mut head.0);
        head.0 = tail;
        head.1 = head.1.max(beyond).max(reach);
    }
    out
}
