//! Conservative uv boxes around pcurves, for exact containment
//! prescreens (PK-743).
//!
//! A point outside a contour's box is outside the contour and outside
//! any polygon sampled from it, so a point-in-contour test can skip the
//! contour without evaluating it. Each box contains the pcurve over its
//! parameter range by construction (no sampling): a segment's ends, an
//! ellipse's centre plus its larger radius, a sinusoid's linear `u` span
//! and `v` amplitude, a positive-weight NURBS curve's control points.

use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};

/// A uv box `((u_min, v_min), (u_max, v_max))`.
pub(crate) type UvBox = ((f64, f64), (f64, f64));

/// A box containing `pc` over its parameter range, or `None` for a form
/// without a cheap conservative bound.
pub(crate) fn pcurve_uv_box(pc: &Pcurve) -> Option<UvBox> {
    let from = |pts: &mut dyn Iterator<Item = (f64, f64)>| -> Option<UvBox> {
        let first = pts.next()?;
        let (mut lo, mut hi) = (first, first);
        for (u, v) in pts {
            lo = (lo.0.min(u), lo.1.min(v));
            hi = (hi.0.max(u), hi.1.max(v));
        }
        Some((lo, hi))
    };
    match &pc.geometry {
        PcurveGeometry::Linear { start, end } => from(&mut [*start, *end].into_iter()),
        PcurveGeometry::Elliptical {
            center,
            u_radius,
            v_radius,
            ..
        } => {
            let r = u_radius.abs().max(v_radius.abs());
            Some(((center.0 - r, center.1 - r), (center.0 + r, center.1 + r)))
        }
        PcurveGeometry::Sinusoidal {
            u_offset,
            v_const,
            v_cos,
            v_sin,
        } => {
            let amp = (v_cos * v_cos + v_sin * v_sin).sqrt();
            let (u0, u1) = (u_offset + pc.t_start, u_offset + pc.t_end);
            Some(((u0.min(u1), v_const - amp), (u0.max(u1), v_const + amp)))
        }
        PcurveGeometry::Nurbs(nurbs) => {
            if nurbs.weights.iter().any(|&w| w <= 0.0) {
                return None;
            }
            from(&mut nurbs.control_points.iter().map(|p| (p.x, p.y)))
        }
        _ => None,
    }
}

/// The union of the boxes of `pcurves`, or `None` when any has none.
pub(crate) fn contour_uv_box<'a>(pcurves: impl IntoIterator<Item = &'a Pcurve>) -> Option<UvBox> {
    let mut acc: Option<UvBox> = None;
    for pc in pcurves {
        let (lo, hi) = pcurve_uv_box(pc)?;
        acc = Some(match acc {
            None => (lo, hi),
            Some((a, b)) => (
                (a.0.min(lo.0), a.1.min(lo.1)),
                (b.0.max(hi.0), b.1.max(hi.1)),
            ),
        });
    }
    acc
}

/// Whether `(u, v)` lies outside `bx` by more than `pad`.
pub(crate) fn outside(bx: &UvBox, u: f64, v: f64, pad: f64) -> bool {
    u < bx.0.0 - pad || v < bx.0.1 - pad || u > bx.1.0 + pad || v > bx.1.1 + pad
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    fn assert_holds(pc: &Pcurve) {
        let bx = pcurve_uv_box(pc).expect("box");
        for k in 0..=1000 {
            let t = pc.t_start + (pc.t_end - pc.t_start) * k as f64 / 1000.0;
            let p = pc.point_at(t);
            assert!(
                !outside(&bx, p.x, p.y, 1e-12),
                "{pc:?} leaves its box at t = {t}"
            );
        }
    }

    #[test]
    fn boxes_contain_their_pcurves() {
        assert_holds(&Pcurve::linear((1.0, -2.0), (-3.0, 4.0), 0.5, 2.0));
        assert_holds(&Pcurve {
            geometry: PcurveGeometry::Elliptical {
                center: (1.0, 2.0),
                u_radius: 3.0,
                v_radius: 1.0,
                rotation: 0.7,
                start_angle: -0.4,
                end_angle: 2.9,
            },
            t_start: 0.0,
            t_end: 1.0,
        });
        assert_holds(&Pcurve {
            geometry: PcurveGeometry::Sinusoidal {
                u_offset: 0.3,
                v_const: 5.0,
                v_cos: 1.5,
                v_sin: -2.0,
            },
            t_start: 0.2,
            t_end: 4.0,
        });
        let nurbs = popkern_core::geometry::NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(2.0, 3.0, 0.0),
                Point3::new(4.0, -1.0, 0.0),
            ],
            weights: vec![1.0, 0.4, 1.0],
        };
        assert_holds(&Pcurve::from_nurbs(nurbs));
    }
}
