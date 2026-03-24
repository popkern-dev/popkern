//! Arc-length reparameterization for path curves.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;

/// Lookup table mapping arc-length to curve parameter.
#[derive(Debug)]
pub struct ArcLengthTable {
    /// Total arc length of the curve.
    pub total_length: f64,
    /// (cumulative_arc_length, curve_parameter) pairs.
    samples: Vec<(f64, f64)>,
}

impl ArcLengthTable {
    /// Build from a NURBS curve with `n_samples` uniform parameter samples.
    pub fn build(curve: &NurbsCurve, n_samples: usize) -> Self {
        let (t_min, t_max) = curve.parameter_range();
        let mut samples = Vec::with_capacity(n_samples + 1);
        let mut cumulative = 0.0;
        samples.push((0.0, t_min));

        for i in 1..=n_samples {
            let t_prev = t_min + (i - 1) as f64 / n_samples as f64 * (t_max - t_min);
            let t_curr = t_min + i as f64 / n_samples as f64 * (t_max - t_min);
            cumulative += curve.length_between(t_prev, t_curr);
            samples.push((cumulative, t_curr));
        }

        Self {
            total_length: cumulative,
            samples,
        }
    }

    /// Convert arc-length `s` to curve parameter `t` via binary search and linear interpolation.
    pub fn s_to_t(&self, s: f64) -> f64 {
        if s <= 0.0 {
            return self.samples.first().map(|p| p.1).unwrap_or(0.0);
        }
        if s >= self.total_length {
            return self.samples.last().map(|p| p.1).unwrap_or(1.0);
        }

        let idx = self.samples.partition_point(|&(arc, _)| arc < s);
        if idx == 0 {
            return self.samples[0].1;
        }

        let (s0, t0) = self.samples[idx - 1];
        let (s1, t1) = self.samples[idx];
        let frac = if (s1 - s0).abs() < 1e-15 {
            0.0
        } else {
            (s - s0) / (s1 - s0)
        };
        t0 + frac * (t1 - t0)
    }

    /// Convert curve parameter `t` to arc-length `s` via binary search and linear interpolation.
    pub fn t_to_s(&self, t: f64) -> f64 {
        let idx = self.samples.partition_point(|&(_, param)| param < t);
        if idx == 0 {
            return 0.0;
        }
        if idx >= self.samples.len() {
            return self.total_length;
        }

        let (s0, t0) = self.samples[idx - 1];
        let (s1, t1) = self.samples[idx];
        let frac = if (t1 - t0).abs() < 1e-15 {
            0.0
        } else {
            (t - t0) / (t1 - t0)
        };
        s0 + frac * (s1 - s0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    fn straight_line() -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn straight_line_total_length() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        // Distance from (0,0,0) to (3,4,0) = 5.0
        assert!(
            (table.total_length - 5.0).abs() < 0.01,
            "total_length should be ~5.0, got {}",
            table.total_length
        );
    }

    #[test]
    fn straight_line_s_to_t_is_linear() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);

        // For a straight line, s_to_t should be linear.
        for i in 0..=10 {
            let s = table.total_length * i as f64 / 10.0;
            let t = table.s_to_t(s);
            let expected_t = i as f64 / 10.0;
            assert!(
                (t - expected_t).abs() < 0.01,
                "s_to_t({s}) = {t}, expected ~{expected_t}"
            );
        }
    }

    #[test]
    fn round_trip_s_to_t_to_s() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 200);

        for i in 1..10 {
            let t_orig = i as f64 / 10.0;
            let s = table.t_to_s(t_orig);
            let t_back = table.s_to_t(s);
            assert!(
                (t_back - t_orig).abs() < 0.02,
                "round-trip failed: t={t_orig} -> s={s} -> t={t_back}"
            );
        }
    }

    #[test]
    fn boundary_s_zero_gives_t_min() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let t = table.s_to_t(0.0);
        assert!(t.abs() < 1e-10, "s=0 should give t_min=0.0, got {t}");
    }

    #[test]
    fn boundary_s_total_gives_t_max() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let t = table.s_to_t(table.total_length);
        assert!(
            (t - 1.0).abs() < 1e-10,
            "s=total should give t_max=1.0, got {t}"
        );
    }

    #[test]
    fn negative_s_clamps_to_t_min() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let t = table.s_to_t(-1.0);
        assert!(t.abs() < 1e-10);
    }

    #[test]
    fn excessive_s_clamps_to_t_max() {
        let curve = straight_line();
        let table = ArcLengthTable::build(&curve, 100);
        let t = table.s_to_t(999.0);
        assert!((t - 1.0).abs() < 1e-10);
    }
}
