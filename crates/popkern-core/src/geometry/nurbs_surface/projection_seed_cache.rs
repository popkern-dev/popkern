//! Content-keyed memoization of a NURBS surface's projection seed
//! field (PK-513).
//!
//! The multi-sheet projector ranks every knot-span seed by its 3D
//! distance to the query before refining the best few. On a coiled
//! thread wall that is ~600 surface evaluations per projection, and a
//! chamfer-on-thread march pays tens of thousands of projections
//! (measured: 52k multisheet projections in one M10 stud-tip build,
//! ~31M seed evaluations of points that never change). The seed
//! parameters and their 3D points are a pure function of the surface,
//! so they are computed once per surface content.
//!
//! Design (the PK-273 Gregory cache pattern):
//! - **Key**: degrees, grid sizes and the bit patterns of every knot,
//!   control point and weight, folded by two independently seeded
//!   64-bit mixers plus the word count. Exact content identity (bit
//!   patterns, so NaN/-0.0 are distinct keys).
//! - **Storage**: thread-local; no locks, each thread warms its own.
//! - **Bound**: cleared when it exceeds [`CACHE_CAPACITY`] entries, a
//!   deterministic policy; a re-miss only re-pays the old cost.
//! - **Bit-identical**: the cached points are the values the uncached
//!   ranking would compute, in the same order.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use super::NurbsSurface;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceEval;

/// Entry bound for the per-thread seed-field map (each field is a few
/// tens of KB on a coiled wall).
const CACHE_CAPACITY: usize = 64;

/// A surface's projection seeds and their 3D points, index-aligned,
/// plus the axis-aligned boxes of its knot spans' control nets.
#[derive(Debug)]
pub struct ProjectionSeedField {
    /// Seed parameters, as [`NurbsSurface::projection_seeds`] returns.
    pub seeds: Vec<(f64, f64)>,
    /// `surface.point_at(seed)` for every seed.
    pub points: Vec<Point3>,
    /// `(min, max)` of each live knot span's local control net. With
    /// positive weights the surface lies in their union; empty when a
    /// weight is not positive (no containment claim).
    pub span_boxes: Vec<(Point3, Point3)>,
}

impl ProjectionSeedField {
    /// Whether `p` may lie within `pad` of the surface: inside some
    /// padded span box. `false` is a proof that every surface point is
    /// farther than `pad` from `p`; `true` claims nothing.
    pub fn may_reach(&self, p: &Point3, pad: f64) -> bool {
        if self.span_boxes.is_empty() {
            return true;
        }
        self.span_boxes.iter().any(|(lo, hi)| {
            p.x >= lo.x - pad
                && p.x <= hi.x + pad
                && p.y >= lo.y - pad
                && p.y <= hi.y + pad
                && p.z >= lo.z - pad
                && p.z <= hi.z + pad
        })
    }
}

/// The control-net box of every live knot span (u-major storage:
/// grid position `(i_u, i_v)` is `i_u * num_v + i_v`).
fn span_boxes(s: &NurbsSurface) -> Vec<(Point3, Point3)> {
    let (p, q) = (s.degree_u, s.degree_v);
    if !s.weights.iter().all(|w| *w > 0.0)
        || s.control_points.len() != s.num_u * s.num_v
        || s.num_u <= p
        || s.num_v <= q
    {
        return Vec::new();
    }
    let live = |knots: &[f64], i: usize| {
        knots
            .get(i + 1)
            .zip(knots.get(i))
            .is_some_and(|(b, a)| b > a)
    };
    let mut boxes = Vec::new();
    for i in p..s.num_u {
        if !live(&s.knots_u, i) {
            continue;
        }
        for j in q..s.num_v {
            if !live(&s.knots_v, j) {
                continue;
            }
            let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
            let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for ii in i - p..=i {
                for jj in j - q..=j {
                    let c = s.control_points[ii * s.num_v + jj];
                    lo = Point3::new(lo.x.min(c.x), lo.y.min(c.y), lo.z.min(c.z));
                    hi = Point3::new(hi.x.max(c.x), hi.y.max(c.y), hi.z.max(c.z));
                }
            }
            boxes.push((lo, hi));
        }
    }
    boxes
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct CacheKey {
    words: u64,
    hash_a: u64,
    hash_b: u64,
}

thread_local! {
    static CACHE: RefCell<HashMap<CacheKey, Option<Rc<ProjectionSeedField>>>> =
        RefCell::new(HashMap::new());
}

/// Two independently seeded multiply-rotate mixers over a word stream.
struct DualMixer {
    a: u64,
    b: u64,
    words: u64,
}

impl DualMixer {
    fn new() -> Self {
        Self {
            a: 0x9e37_79b9_7f4a_7c15,
            b: 0x2545_f491_4f6c_dd1d,
            words: 0,
        }
    }

    fn word(&mut self, w: u64) {
        self.a = (self.a ^ w)
            .wrapping_mul(0xff51_afd7_ed55_8ccd)
            .rotate_left(29);
        self.b = (self.b ^ w.rotate_left(17))
            .wrapping_mul(0xc4ce_b9fe_1a85_ec53)
            .rotate_left(31);
        self.words += 1;
    }

    fn key(self) -> CacheKey {
        CacheKey {
            words: self.words,
            hash_a: self.a,
            hash_b: self.b,
        }
    }
}

fn content_key(surface: &NurbsSurface) -> CacheKey {
    let mut m = DualMixer::new();
    for n in [
        surface.degree_u,
        surface.degree_v,
        surface.num_u,
        surface.num_v,
        surface.knots_u.len(),
        surface.knots_v.len(),
    ] {
        m.word(n as u64);
    }
    for x in surface
        .knots_u
        .iter()
        .chain(&surface.knots_v)
        .chain(&surface.weights)
    {
        m.word(x.to_bits());
    }
    for p in &surface.control_points {
        m.word(p.x.to_bits());
        m.word(p.y.to_bits());
        m.word(p.z.to_bits());
    }
    m.key()
}

impl NurbsSurface {
    /// Exact content fingerprint (degrees, grid sizes, and the bit
    /// patterns of every knot, control point and weight) for
    /// content-keyed memoization by geometry consumers.
    pub(crate) fn content_fingerprint(&self) -> (u64, u64, u64) {
        let k = content_key(self);
        (k.words, k.hash_a, k.hash_b)
    }

    /// The projection seeds with their 3D points, memoized by surface
    /// content (see the module docs). `None` when the chart has too
    /// few spans to need seeds (as [`NurbsSurface::projection_seeds`]).
    pub fn projection_seed_field(&self) -> Option<Rc<ProjectionSeedField>> {
        let key = content_key(self);
        if let Some(hit) = CACHE.with(|c| c.borrow().get(&key).cloned()) {
            return hit;
        }
        let field = self.projection_seeds().map(|seeds| {
            let points = seeds.iter().map(|&(u, v)| self.point_at(u, v)).collect();
            Rc::new(ProjectionSeedField {
                seeds,
                points,
                span_boxes: span_boxes(self),
            })
        });
        CACHE.with(|c| {
            let mut map = c.borrow_mut();
            if map.len() >= CACHE_CAPACITY {
                map.clear();
            }
            map.insert(key, field.clone());
        });
        field
    }
}

#[cfg(test)]
#[path = "projection_seed_cache_tests.rs"]
mod tests;
