//! Surface jet enclosures: every partial up to [`MAX_JET_ORDER`] over a
//! parameter box, dispatched per surface kind.

use super::analytic_jet::{AnalyticJetSource, frame_of};
use super::interval::Iv;
use super::interval_vec::IvVec3;
use super::nurbs_jet::NurbsJetSource;
use super::{CertifyDecline, UnsupportedGeometry};
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::SurfaceGeometry;

/// Highest total partial order a jet carries: the order-4 Taylor form
/// of the second fundamental form needs `S` to order 6.
pub const MAX_JET_ORDER: usize = 6;

/// Number of partials of total order `≤` [`MAX_JET_ORDER`].
pub const JET_SLOTS: usize = (MAX_JET_ORDER + 1) * (MAX_JET_ORDER + 2) / 2;

/// The partial orders `(∂u, ∂v)` a jet carries, by increasing total
/// order and ∂u first within a total (the Leibniz quotient relies on
/// this ordering).
pub const JET_ORDERS: [(usize, usize); JET_SLOTS] = jet_orders();

const fn jet_orders() -> [(usize, usize); JET_SLOTS] {
    let mut out = [(0, 0); JET_SLOTS];
    let (mut total, mut k) = (0, 0);
    while total <= MAX_JET_ORDER {
        let mut b = 0;
        while b <= total {
            out[k] = (total - b, b);
            k += 1;
            b += 1;
        }
        total += 1;
    }
    out
}

/// Slot of `(∂u^a, ∂v^b)` in [`JET_ORDERS`].
pub fn jet_slot(a: usize, b: usize) -> usize {
    let total = a + b;
    // Orders of total t start at t(t+1)/2 and list ∂u first.
    total * (total + 1) / 2 + b
}

/// Enclosures of `S` and its partials up to [`MAX_JET_ORDER`] over a
/// box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceJet {
    partials: [IvVec3; JET_SLOTS],
}

impl SurfaceJet {
    /// A jet with every partial unbounded.
    pub fn entire() -> Self {
        Self {
            partials: [IvVec3::ENTIRE; JET_SLOTS],
        }
    }

    /// A jet with every partial exactly zero.
    pub fn zero() -> Self {
        Self {
            partials: [IvVec3::ZERO; JET_SLOTS],
        }
    }

    /// Enclosure of `∂u^a ∂v^b S` (`a + b ≤` [`MAX_JET_ORDER`]; higher
    /// orders read as unbounded).
    pub fn partial(&self, a: usize, b: usize) -> IvVec3 {
        if a + b > MAX_JET_ORDER {
            return IvVec3::ENTIRE;
        }
        self.partials[jet_slot(a, b)]
    }

    /// Set the enclosure of `∂u^a ∂v^b S`.
    pub fn set(&mut self, a: usize, b: usize, value: IvVec3) {
        if a + b <= MAX_JET_ORDER {
            self.partials[jet_slot(a, b)] = value;
        }
    }

    /// Component-wise hull of two jets.
    pub fn hull(&self, other: &SurfaceJet) -> SurfaceJet {
        let mut out = *self;
        for (a, b) in out.partials.iter_mut().zip(other.partials.iter()) {
            *a = a.hull(*b);
        }
        out
    }

    /// True when every partial is finite.
    pub fn is_finite(&self) -> bool {
        self.partials.iter().all(|p| p.is_finite())
    }
}

/// A surface prepared for rigorous jet enclosure.
#[derive(Debug, Clone)]
pub enum SurfaceJetSource {
    /// Closed-form analytic kind.
    Analytic(AnalyticJetSource),
    /// Rational tensor-product B-spline (NURBS, Gregory exact forms,
    /// polyhedral-spline Bézier nets).
    Spline(Box<NurbsJetSource>),
}

/// Relative weight below which a Gregory exact form's entry is an
/// ε-regularized corner (the exact form has weight 0 there).
const REGULARIZED_CORNER_WEIGHT: f64 = 1e-12;

impl SurfaceJetSource {
    /// Prepare `surface`, or say why it cannot be certified.
    pub fn new(surface: &SurfaceGeometry) -> Result<Self, CertifyDecline> {
        Ok(match surface {
            SurfaceGeometry::Plane {
                origin,
                u_axis,
                v_axis,
                ..
            } => SurfaceJetSource::Analytic(AnalyticJetSource::Plane {
                origin: *origin,
                u_axis: *u_axis,
                v_axis: *v_axis,
            }),
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => SurfaceJetSource::Analytic(AnalyticJetSource::Cylinder {
                origin: *origin,
                axis: *axis,
                frame: frame_of(*axis),
                radius: *radius,
            }),
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } => SurfaceJetSource::Analytic(AnalyticJetSource::Cone {
                apex: *apex,
                axis: *axis,
                frame: frame_of(*axis),
                half_angle: *half_angle,
            }),
            SurfaceGeometry::Sphere { center, radius } => {
                SurfaceJetSource::Analytic(AnalyticJetSource::Sphere {
                    center: *center,
                    radius: *radius,
                })
            }
            SurfaceGeometry::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            } => SurfaceJetSource::Analytic(AnalyticJetSource::Torus {
                center: *center,
                axis: *axis,
                frame: frame_of(*axis),
                major: *major_radius,
                minor: *minor_radius,
            }),
            SurfaceGeometry::Nurbs(nurbs) => {
                SurfaceJetSource::Spline(Box::new(NurbsJetSource::new(nurbs, &[])?))
            }
            SurfaceGeometry::Gregory(patch) => {
                let exact = gregory_exact_form(patch)?;
                let max_weight = exact.weights.iter().copied().fold(0.0_f64, f64::max);
                let regularized: Vec<usize> = exact
                    .weights
                    .iter()
                    .enumerate()
                    .filter(|&(_, &w)| w > 0.0 && w <= REGULARIZED_CORNER_WEIGHT * max_weight)
                    .map(|(k, _)| k)
                    .collect();
                SurfaceJetSource::Spline(Box::new(NurbsJetSource::new(&exact, &regularized)?))
            }
            SurfaceGeometry::PolyhedralSpline(patch) => {
                let order = patch.net.len();
                let degree = order - 1;
                let mut control_points = Vec::with_capacity(order * order);
                for i_u in 0..order {
                    for i_v in 0..order {
                        control_points.push(patch.net[i_v][i_u]);
                    }
                }
                let knots: Vec<f64> = (0..2 * order)
                    .map(|k| if k < order { 0.0 } else { 1.0 })
                    .collect();
                let bezier = NurbsSurface {
                    degree_u: degree,
                    degree_v: degree,
                    knots_u: knots.clone(),
                    knots_v: knots,
                    control_points,
                    num_u: order,
                    num_v: order,
                    weights: vec![1.0; order * order],
                };
                SurfaceJetSource::Spline(Box::new(NurbsJetSource::new(&bezier, &[])?))
            }
        })
    }

    /// Interior parameter breaks (knots) in u and v; empty for the
    /// analytic kinds.
    pub fn breaks(&self) -> (Vec<f64>, Vec<f64>) {
        match self {
            SurfaceJetSource::Analytic(_) => (Vec::new(), Vec::new()),
            SurfaceJetSource::Spline(source) => source.breaks(),
        }
    }

    /// Enclosure of every partial up to [`MAX_JET_ORDER`] over the box
    /// `u × v`,
    /// which must lie in one knot span.
    pub fn jet(&mut self, u: Iv, v: Iv) -> Result<SurfaceJet, CertifyDecline> {
        let jet = match self {
            SurfaceJetSource::Analytic(source) => source.jet(u, v),
            SurfaceJetSource::Spline(source) => source.jet(u, v)?,
        };
        if jet.is_finite() {
            Ok(jet)
        } else {
            Err(CertifyDecline::NonFinite)
        }
    }

    /// Knot spans the box meets (`None` for the analytic kinds, which
    /// have a single smooth branch).
    pub fn spans_touching(
        &self,
        u: Iv,
        v: Iv,
    ) -> Result<Vec<Option<(usize, usize)>>, CertifyDecline> {
        match self {
            SurfaceJetSource::Analytic(_) => Ok(vec![None]),
            SurfaceJetSource::Spline(source) => {
                Ok(source.spans_touching(u, v)?.into_iter().map(Some).collect())
            }
        }
    }

    /// Enclosure over the box of one span's (extended) polynomial, or
    /// of the analytic surface.
    pub fn jet_in(
        &mut self,
        u: Iv,
        v: Iv,
        span: Option<(usize, usize)>,
    ) -> Result<SurfaceJet, CertifyDecline> {
        let jet = match (self, span) {
            (SurfaceJetSource::Analytic(source), _) => source.jet(u, v),
            (SurfaceJetSource::Spline(source), Some(span)) => source.jet_in_span(u, v, span)?,
            (SurfaceJetSource::Spline(source), None) => source.jet(u, v)?,
        };
        if jet.is_finite() {
            Ok(jet)
        } else {
            Err(CertifyDecline::NonFinite)
        }
    }

    /// Enclosure over any box, hulled across knot spans; the flag
    /// says whether the box straddled a knot (see
    /// [`NurbsJetSource::jet_piecewise`]).
    pub fn jet_piecewise(&mut self, u: Iv, v: Iv) -> Result<(SurfaceJet, bool), CertifyDecline> {
        let (jet, straddled) = match self {
            SurfaceJetSource::Analytic(source) => (source.jet(u, v), false),
            SurfaceJetSource::Spline(source) => source.jet_piecewise(u, v)?,
        };
        if jet.is_finite() {
            Ok((jet, straddled))
        } else {
            Err(CertifyDecline::NonFinite)
        }
    }
}

/// The exact rational form of a Gregory patch in its own `(u, v)`
/// chart. Variants whose conversion is not exact in that chart are
/// declined: the dual-twist containment conversion is a midpoint
/// approximation, and the triangular conversion uses a different
/// parameter mapping from the patch's evaluator.
fn gregory_exact_form(patch: &GregoryPatch) -> Result<NurbsSurface, CertifyDecline> {
    Ok(match patch {
        GregoryPatch::Rectangular(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularDeg5(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularBiquinticFull(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularBiseptic(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularTwoStripe(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularBiNonic(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularBiNonicDualJet(r) => r.to_nurbs_exact(),
        GregoryPatch::RectangularBiNonicDualTwist(_) => {
            return Err(CertifyDecline::UnsupportedSurface(
                UnsupportedGeometry::DualTwistGregory,
            ));
        }
        GregoryPatch::Triangular(_) => {
            return Err(CertifyDecline::UnsupportedSurface(
                UnsupportedGeometry::TriangularGregory,
            ));
        }
    })
}
