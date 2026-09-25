//! Surface jet enclosures: every partial up to order 3 over a
//! parameter box, dispatched per surface kind.

use super::CertifyDecline;
use super::analytic_jet::{AnalyticJetSource, frame_of};
use super::interval::Iv;
use super::interval_vec::IvVec3;
use super::nurbs_jet::NurbsJetSource;
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::SurfaceGeometry;

/// The partial orders `(∂u, ∂v)` a jet carries, by increasing total
/// order (the Leibniz quotient relies on this ordering).
pub const JET_ORDERS: [(usize, usize); 10] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (2, 0),
    (1, 1),
    (0, 2),
    (3, 0),
    (2, 1),
    (1, 2),
    (0, 3),
];

/// Slot of `(∂u^a, ∂v^b)` in [`JET_ORDERS`].
pub fn jet_slot(a: usize, b: usize) -> usize {
    let total = a + b;
    // Orders of total t start at t(t+1)/2 and list ∂u first.
    total * (total + 1) / 2 + b
}

/// Enclosures of `S` and its partials up to order 3 over a box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceJet {
    partials: [IvVec3; 10],
}

impl SurfaceJet {
    /// A jet with every partial unbounded.
    pub fn entire() -> Self {
        Self {
            partials: [IvVec3::ENTIRE; 10],
        }
    }

    /// A jet with every partial exactly zero.
    pub fn zero() -> Self {
        Self {
            partials: [IvVec3::ZERO; 10],
        }
    }

    /// Enclosure of `∂u^a ∂v^b S` (`a + b ≤ 3`; higher orders read as
    /// unbounded).
    pub fn partial(&self, a: usize, b: usize) -> IvVec3 {
        if a + b > 3 {
            return IvVec3::ENTIRE;
        }
        self.partials[jet_slot(a, b)]
    }

    /// Set the enclosure of `∂u^a ∂v^b S`.
    pub fn set(&mut self, a: usize, b: usize, value: IvVec3) {
        if a + b <= 3 {
            self.partials[jet_slot(a, b)] = value;
        }
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

    /// Enclosure of every partial up to order 3 over the box `u × v`.
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
                "dual-twist Gregory (midpoint containment form)",
            ));
        }
        GregoryPatch::Triangular(_) => {
            return Err(CertifyDecline::UnsupportedSurface(
                "triangular Gregory (exact form uses another chart)",
            ));
        }
    })
}
