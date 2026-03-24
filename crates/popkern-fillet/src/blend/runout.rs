//! Fillet runout — taper radius to zero at unfilleted endpoints.
//!
//! When a filleted edge meets an unfilleted edge, the fillet must taper
//! smoothly to zero radius. This module provides the radius function
//! that can be fed into [`crate::center_curve::variable_radius_center_curve`].

/// Which end of the fillet tapers to zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunoutEnd {
    /// Taper to zero at the start (t = 0).
    Start,
    /// Taper to zero at the end (t = 1).
    End,
    /// Taper to zero at both ends.
    Both,
}

/// Compute the runout radius function for a fillet that tapers to zero.
///
/// Returns a function `f(t) -> radius` where `t` ranges from 0.0 to 1.0.
/// The radius linearly tapers from `base_radius` to 0.0 at the specified
/// end(s).
///
/// # Examples
///
/// ```
/// use popkern_fillet::blend::runout::{RunoutEnd, runout_radius_fn};
///
/// let f = runout_radius_fn(2.0, RunoutEnd::Start);
/// assert!((f(0.0) - 0.0).abs() < 1e-10);
/// assert!((f(1.0) - 2.0).abs() < 1e-10);
/// ```
pub fn runout_radius_fn(base_radius: f64, runout_end: RunoutEnd) -> Box<dyn Fn(f64) -> f64> {
    match runout_end {
        RunoutEnd::Start => {
            // Radius goes from 0 at t=0 to base_radius at t=1.
            Box::new(move |t: f64| base_radius * t.clamp(0.0, 1.0))
        }
        RunoutEnd::End => {
            // Radius goes from base_radius at t=0 to 0 at t=1.
            Box::new(move |t: f64| base_radius * (1.0 - t.clamp(0.0, 1.0)))
        }
        RunoutEnd::Both => {
            // Radius is 0 at both ends, peaks at base_radius at t=0.5.
            // Uses a parabolic profile: r(t) = 4 * base_radius * t * (1 - t).
            Box::new(move |t: f64| {
                let t = t.clamp(0.0, 1.0);
                4.0 * base_radius * t * (1.0 - t)
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runout_start() {
        let f = runout_radius_fn(2.0, RunoutEnd::Start);
        assert!((f(0.0) - 0.0).abs() < 1e-10, "r(0) should be 0");
        assert!((f(1.0) - 2.0).abs() < 1e-10, "r(1) should be base_radius");
        assert!((f(0.5) - 1.0).abs() < 1e-10, "r(0.5) should be 1.0");
    }

    #[test]
    fn runout_end() {
        let f = runout_radius_fn(3.0, RunoutEnd::End);
        assert!((f(0.0) - 3.0).abs() < 1e-10, "r(0) should be base_radius");
        assert!((f(1.0) - 0.0).abs() < 1e-10, "r(1) should be 0");
        assert!((f(0.5) - 1.5).abs() < 1e-10, "r(0.5) should be 1.5");
    }

    #[test]
    fn runout_both() {
        let f = runout_radius_fn(2.0, RunoutEnd::Both);
        assert!((f(0.0) - 0.0).abs() < 1e-10, "r(0) should be 0");
        assert!((f(1.0) - 0.0).abs() < 1e-10, "r(1) should be 0");
        // Peak at t=0.5: 4 * 2.0 * 0.25 = 2.0
        assert!((f(0.5) - 2.0).abs() < 1e-10, "r(0.5) should be base_radius");
    }

    #[test]
    fn runout_clamped_outside_range() {
        let f = runout_radius_fn(1.0, RunoutEnd::Start);
        // Values outside [0,1] should be clamped.
        assert!((f(-0.5) - 0.0).abs() < 1e-10, "r(-0.5) should clamp to 0");
        assert!((f(1.5) - 1.0).abs() < 1e-10, "r(1.5) should clamp to 1.0");
    }

    #[test]
    fn runout_zero_radius() {
        let f = runout_radius_fn(0.0, RunoutEnd::End);
        assert!((f(0.0) - 0.0).abs() < 1e-10);
        assert!((f(0.5) - 0.0).abs() < 1e-10);
        assert!((f(1.0) - 0.0).abs() < 1e-10);
    }
}
