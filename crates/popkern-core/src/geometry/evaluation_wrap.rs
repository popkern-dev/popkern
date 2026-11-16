//! The parameter wrap a surface applies inside its own evaluation,
//! computed once per query instead of once per evaluation (PK-1223).
//!
//! A NURBS surface in periodic form wraps an out-of-domain parameter
//! back into its knot domain before evaluating, and decides whether it
//! is periodic from its representation: a scan of its knots and wrapped
//! control rows, linear in the net. A Gauss–Newton walk on a free
//! periodic axis steps out of the domain on nearly every evaluation, so
//! on a 6.5-turn coil wall the scan ran 53 million times in one grind
//! (77% of its 190 s). A walk reads the wrap once and applies it before
//! each evaluation; the surface then sees an in-domain parameter and
//! leaves it alone, so the evaluation is bit-identical to the unwrapped
//! one.

/// One axis's wrap: the evaluation domain `[lo, hi]` and its period.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxisWrap {
    /// Low end of the evaluation domain.
    pub lo: f64,
    /// High end of the evaluation domain.
    pub hi: f64,
    /// The period an out-of-domain parameter is wrapped by.
    pub period: f64,
}

impl AxisWrap {
    /// `t` brought into `[lo, hi]` by whole periods; `t` unchanged when
    /// it is inside already or not finite. A wrap that rounds past the
    /// domain returns `t` too, leaving the surface's own wrap to act, so
    /// evaluating at the result is always the surface's evaluation at
    /// `t`.
    pub fn apply(&self, t: f64) -> f64 {
        let Some(wrapped) = self.wrap(t) else {
            return t;
        };
        if (self.lo..=self.hi).contains(&wrapped) {
            wrapped
        } else {
            t
        }
    }

    /// The surface's own wrap of `t`: `None` when `t` is inside the
    /// domain or not finite, else `lo + (t − lo) mod period`.
    pub fn wrap(&self, t: f64) -> Option<f64> {
        if (self.lo..=self.hi).contains(&t) || !t.is_finite() {
            return None;
        }
        Some(self.lo + (t - self.lo).rem_euclid(self.period))
    }
}

/// The wrap a surface's evaluation applies to `(u, v)`: per axis, `None`
/// where evaluation takes the parameter as given.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct EvaluationWrap {
    /// The wrap along `u`.
    pub u: Option<AxisWrap>,
    /// The wrap along `v`.
    pub v: Option<AxisWrap>,
}

impl EvaluationWrap {
    /// `(u, v)` with each axis's wrap applied (see [`AxisWrap::apply`]).
    pub fn apply(&self, u: f64, v: f64) -> (f64, f64) {
        (
            self.u.map_or(u, |w| w.apply(u)),
            self.v.map_or(v, |w| w.apply(v)),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inside_and_non_finite_parameters_are_left_alone() {
        let w = AxisWrap {
            lo: 0.0,
            hi: 1.0,
            period: 1.0,
        };
        assert_eq!(w.apply(0.25), 0.25);
        assert_eq!(w.apply(1.0), 1.0);
        assert!(w.apply(f64::NAN).is_nan());
        assert_eq!(w.apply(f64::INFINITY), f64::INFINITY);
    }

    #[test]
    fn outside_parameters_wrap_by_whole_periods() {
        let w = AxisWrap {
            lo: 2.0,
            hi: 5.0,
            period: 3.0,
        };
        assert_eq!(w.apply(6.5), 3.5);
        assert_eq!(w.apply(-0.5), 2.5);
        let both = EvaluationWrap {
            u: Some(w),
            v: None,
        };
        assert_eq!(both.apply(6.5, 9.0), (3.5, 9.0));
    }
}
