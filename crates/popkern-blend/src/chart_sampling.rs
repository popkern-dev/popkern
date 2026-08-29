//! Period-honest interior chart sampling shared by every cap/support
//! normal measure.
//!
//! Periodic chart directions must WRAP a parameter into the domain
//! before the interior clamp. A support pcurve legally chained onto an
//! adjacent period branch (u = 7.78 or u = −1.50 on a cylinder chart
//! of [0, 2π]) evaluates exactly, but clamping such a parameter to the
//! domain edge samples the normal up to a full period away. This is
//! the PK-125 defect class, re-proven twice during PK-255: a 0.888
//! live-equivalence decline (fixed at one author, cd866e59) and a
//! 0.928 acceptance decline that blocked a crossing-free pole-search
//! serve while a dozen sibling copies of the same clamp survived.
//! Every measure that evaluates surface normals at pcurve-derived
//! parameters must sample through this module instead of a local
//! clamp.

use popkern_core::geometry::surface::SurfaceEval;

/// Domain-relative margin pulling a boundary sample just inside the
/// chart so edge-degenerate partials do not poison the normal.
const INTERIOR_MARGIN: f64 = 1.0e-5;

/// Wrap periodic directions into the chart domain, then clamp both
/// directions to the interior margin.
pub(crate) fn periodic_interior_uv(surface: &impl SurfaceEval, u: f64, v: f64) -> (f64, f64) {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let (period_u, period_v) = surface.parameter_periodicity();
    (
        interior_parameter(wrap_periodic(u, u_min, u_max, period_u), u_min, u_max),
        interior_parameter(wrap_periodic(v, v_min, v_max, period_v), v_min, v_max),
    )
}

/// Wrap a single periodic parameter into `[min, max]`; non-periodic
/// directions pass through unchanged.
pub(crate) fn wrap_periodic(value: f64, min: f64, max: f64, period: Option<f64>) -> f64 {
    let Some(period) = period else {
        return value;
    };
    if !(value.is_finite() && min.is_finite() && max.is_finite()) || period <= 0.0 {
        return value;
    }
    let mut wrapped = value;
    while wrapped < min {
        wrapped += period;
    }
    while wrapped > max {
        wrapped -= period;
    }
    wrapped
}

/// Clamp a parameter to the chart interior by the domain-relative
/// margin. Infinite or degenerate domains pass through unchanged.
pub(crate) fn interior_parameter(value: f64, min: f64, max: f64) -> f64 {
    if !value.is_finite() || !min.is_finite() || !max.is_finite() || max <= min {
        return value;
    }
    let margin = ((max - min).abs() * INTERIOR_MARGIN).max(INTERIOR_MARGIN);
    value.clamp(min + margin, max - margin)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};

    fn cylinder() -> SurfaceGeometry {
        SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 0.05,
        }
    }

    #[test]
    fn adjacent_branch_parameters_wrap_before_clamping() {
        let surface = cylinder();
        let ((u_min, u_max), _) = surface.parameter_range();
        if !(u_min.is_finite() && u_max.is_finite()) {
            // Chart is unbounded in u; wrapping is a pass-through and
            // the clamp defect cannot occur on this representation.
            return;
        }
        let period = u_max - u_min;
        let (u_wrapped, _) = periodic_interior_uv(&surface, u_min - 1.499, 0.5);
        let expected = u_min - 1.499 + period;
        assert!(
            (u_wrapped - expected).abs() <= 1.0e-9,
            "adjacent-branch u must wrap into the chart, not clamp: got {u_wrapped}, expected {expected}"
        );
    }

    #[test]
    fn non_periodic_directions_still_clamp() {
        let surface = cylinder();
        let (_, (v_min, v_max)) = surface.parameter_range();
        if !(v_min.is_finite() && v_max.is_finite()) {
            return;
        }
        let (_, v) = periodic_interior_uv(&surface, 0.5, v_max + 10.0);
        assert!(v < v_max, "out-of-domain non-periodic v must clamp inside");
    }
}
