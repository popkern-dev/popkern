//! One refinement level: the hard system, the fixed rim rows, the rim
//! samples with their neighbour jets, and the solves over the null
//! space.

use std::sync::Arc;

use super::jets::{NeighbourJet, NeighbourSampler};
use super::measure::RimSample;
use super::nullspace::{HardSystem, hard_system};
use super::rim::side_rows;
use super::space::{SectorSpace, gauss_legendre};
use super::{RimHalf, VariationalCapError, VariationalCapSide};
use crate::geometry::{Point3, Vec3};

/// Model ↔ normalised frame: `p_model = origin + scale · p`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Frame {
    pub origin: Point3,
    pub scale: f64,
}

impl Frame {
    pub fn to_local(self, p: Point3) -> Vec3 {
        Vec3::new(
            (p.x - self.origin.x) / self.scale,
            (p.y - self.origin.y) / self.scale,
            (p.z - self.origin.z) / self.scale,
        )
    }
}

/// The per-level state.
pub(super) struct Level {
    pub space: SectorSpace,
    pub hs: Arc<HardSystem>,
    /// Fixed (rim) values per coordinate, in `hs.fixed` order.
    pub fixed: [Vec<f64>; 3],
    /// Particular free values per coordinate.
    pub xp: [Vec<f64>; 3],
    /// Rim samples with neighbour jets (normalised frame).
    pub samples: Vec<RimSample>,
}

/// Build a level: rim rows, consistency check, samples.
pub(super) fn build_level(
    space: SectorSpace,
    sides: &[VariationalCapSide],
    frame: &Frame,
    splits: &[f64],
) -> Result<Level, VariationalCapError> {
    let n = space.n;
    let m = space.m;
    let hs = hard_system(&space).ok_or(VariationalCapError::Singular {
        stage: "spoke null space",
    })?;
    let mut full: [Vec<f64>; 3] = std::array::from_fn(|_| vec![f64::NAN; space.total()]);
    fn put(full: &mut [Vec<f64>; 3], idx: usize, p: Vec3) {
        for (c, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            full[c][idx] = v;
        }
    }
    let mut params = Vec::with_capacity(n);
    for (s, side) in sides.iter().enumerate() {
        let rows = side_rows(&space, side, s, splits[s])?;
        let prev = (s + n - 1) % n;
        for (j, &p) in rows.trailing.iter().enumerate() {
            put(&mut full, space.index(s, m - 1, j), frame.to_local(p));
        }
        // The leading row's last point is corner c_{s−1}, written by side
        // s−1's trailing row (closure was checked against the tolerance
        // before any level was built).
        for (i, &p) in rows.leading.iter().enumerate().take(m - 1) {
            put(&mut full, space.index(prev, i, m - 1), frame.to_local(p));
        }
        params.push(rows.params);
    }
    let fixed: [Vec<f64>; 3] =
        std::array::from_fn(|c| hs.fixed.iter().map(|&i| full[c][i]).collect());
    if fixed.iter().flatten().any(|v| !v.is_finite()) {
        return Err(VariationalCapError::Inconsistent {
            stage: "rim rows incomplete",
            residual: f64::INFINITY,
        });
    }
    let mut xp: [Vec<f64>; 3] = Default::default();
    for c in 0..3 {
        let (x, residual) = hs.particular(&fixed[c]);
        if residual > 1.0e-9 {
            return Err(VariationalCapError::Inconsistent {
                stage: "spoke rows against rim rows",
                residual,
            });
        }
        xp[c] = x;
    }
    let samples = rim_samples(&space, sides, frame, &params)?;
    Ok(Level {
        space,
        hs,
        fixed,
        xp,
        samples,
    })
}

/// Gauss samples on every rim half (`degree + 3` per span) plus both
/// ends, each with its neighbour jet in the normalised frame.
fn rim_samples(
    space: &SectorSpace,
    sides: &[VariationalCapSide],
    frame: &Frame,
    params: &[(f64, f64, f64)],
) -> Result<Vec<RimSample>, VariationalCapError> {
    let n = space.n;
    let nodes = gauss_legendre(space.degree + 3);
    let mut xs: Vec<f64> = vec![0.0];
    for (a, b) in space.spans() {
        xs.extend(nodes.iter().map(|&(x, _)| a + (b - a) * x));
    }
    xs.push(1.0);
    let mut out = Vec::new();
    for (s, side) in sides.iter().enumerate() {
        let (t0, tm, t1) = params[s];
        for half in [RimHalf::Leading, RimHalf::Trailing] {
            let mut sampler = NeighbourSampler::new(&side.neighbour, side.neighbour_reversed)
                .map_err(|reason| VariationalCapError::NeighbourJet {
                    side: s,
                    at: tm,
                    reason,
                })?;
            let map = super::rim::HalfMap::new(t0, tm, t1);
            let (sector, trailing) = match half {
                RimHalf::Trailing => (s, true),
                RimHalf::Leading => ((s + n - 1) % n, false),
            };
            for &x in &xs {
                let t = map.t(trailing, x);
                let p = crate::geometry::curve::CurveEval::point_at(&side.curve, t);
                let model =
                    sampler
                        .jet_at(&p)
                        .map_err(|reason| VariationalCapError::NeighbourJet {
                            side: s,
                            at: t,
                            reason,
                        })?;
                out.push(RimSample {
                    side: s,
                    half,
                    sector,
                    x,
                    t,
                    jet: localise(&model, frame),
                });
            }
        }
    }
    Ok(out)
}

/// A model-frame jet in the normalised frame (curvature × scale).
fn localise(jet: &NeighbourJet, frame: &Frame) -> NeighbourJet {
    let mut form = jet.form;
    for row in &mut form {
        for v in row.iter_mut() {
            *v *= frame.scale;
        }
    }
    let foot = frame.to_local(jet.foot);
    NeighbourJet {
        foot: Point3::new(foot.x, foot.y, foot.z),
        normal: jet.normal,
        form,
    }
}

/// Full coefficient vectors from null-space coordinates.
pub(super) fn assemble(level: &Level, y: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
    let hs = &level.hs;
    let k = hs.dim();
    std::array::from_fn(|c| {
        let mut x = vec![0.0; level.space.total()];
        for (r, &i) in hs.free.iter().enumerate() {
            x[i] = level.xp[c][r] + hs.dot_row(r, &y[c][..k]);
        }
        for (r, &i) in hs.fixed.iter().enumerate() {
            x[i] = level.fixed[c][r];
        }
        x
    })
}
