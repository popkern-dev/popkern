//! Read-only query accessors for [`Sketch`] elements.
//!
//! Provides iterators, element enumeration, and curve evaluation
//! methods used by downstream crates (e.g., `popkern-sketch`).
//!
//! Split from `sketch_ops.rs` to stay within the 500-line file limit.

use crate::Sketch;
use crate::elements::{
    ArcHandle, ElementHandle, EllipseHandle, HyperbolaHandle, LineHandle, ParabolaHandle,
    PointHandle, SplineHandle,
};
use crate::equations::spline_eval;

impl Sketch {
    /// Iterate over all line handles in the sketch.
    pub fn line_handles(&self) -> impl Iterator<Item = LineHandle> + '_ {
        self.lines.keys()
    }

    /// Iterate over all arc handles in the sketch.
    pub fn arc_handles(&self) -> impl Iterator<Item = ArcHandle> + '_ {
        self.arcs.keys()
    }

    /// Iterate over all ellipse handles in the sketch.
    pub fn ellipse_handles(&self) -> impl Iterator<Item = EllipseHandle> + '_ {
        self.ellipses.keys()
    }

    /// Iterate over all parabola handles in the sketch.
    pub fn parabola_handles(&self) -> impl Iterator<Item = ParabolaHandle> + '_ {
        self.parabolas.keys()
    }

    /// Iterate over all hyperbola handles in the sketch.
    pub fn hyperbola_handles(&self) -> impl Iterator<Item = HyperbolaHandle> + '_ {
        self.hyperbolas.keys()
    }

    /// Iterate over all spline handles in the sketch.
    pub fn spline_handles(&self) -> impl Iterator<Item = SplineHandle> + '_ {
        self.splines.keys()
    }

    /// Return `true` if the given element is marked as construction geometry.
    ///
    /// Returns `None` if the handle is invalid.
    pub fn is_construction(&self, element: &ElementHandle) -> Option<bool> {
        match element {
            ElementHandle::Point(h) => self.points.get(*h).map(|p| p.construction),
            ElementHandle::Line(h) => self.lines.get(*h).map(|l| l.construction),
            ElementHandle::Arc(h) => self.arcs.get(*h).map(|a| a.construction),
            ElementHandle::Ellipse(h) => self.ellipses.get(*h).map(|e| e.construction),
            ElementHandle::Parabola(h) => self.parabolas.get(*h).map(|p| p.construction),
            ElementHandle::Hyperbola(h) => self.hyperbolas.get(*h).map(|h| h.construction),
            ElementHandle::Spline(h) => self.splines.get(*h).map(|s| s.construction),
        }
    }

    /// Evaluate a spline at parameter `t`, returning `(x, y)`.
    ///
    /// Returns `None` if the handle is invalid.
    pub fn spline_point_at(&self, handle: SplineHandle, t: f64) -> Option<(f64, f64)> {
        let spline = self.splines.get(handle)?;
        let result = spline_eval::eval_point(spline, t, self);
        Some((result.x, result.y))
    }

    /// Return the knot vector for a spline, or `None` if the handle is
    /// invalid.
    pub fn spline_knots(&self, handle: SplineHandle) -> Option<&[f64]> {
        self.splines.get(handle).map(|s| s.knots.as_slice())
    }

    /// Return the degree of a spline, or `None` if the handle is invalid.
    pub fn spline_degree(&self, handle: SplineHandle) -> Option<usize> {
        self.splines.get(handle).map(|s| s.degree)
    }

    /// Return whether a spline is periodic, or `None` if the handle is
    /// invalid.
    pub fn spline_periodic(&self, handle: SplineHandle) -> Option<bool> {
        self.splines.get(handle).map(|s| s.periodic)
    }

    /// Return the weight vector for a spline, or `None` if the handle is
    /// invalid.
    pub fn spline_weights(&self, handle: SplineHandle) -> Option<&[f64]> {
        self.splines.get(handle).map(|s| s.weights.as_slice())
    }

    /// Iterate over all non-construction curve element handles.
    ///
    /// Returns lines, arcs, ellipses, parabolas, hyperbolas, and splines
    /// that are not marked as construction geometry.
    pub fn non_construction_curves(&self) -> Vec<ElementHandle> {
        let mut out = Vec::new();
        for (h, l) in &self.lines {
            if !l.construction {
                out.push(ElementHandle::Line(h));
            }
        }
        for (h, a) in &self.arcs {
            if !a.construction {
                out.push(ElementHandle::Arc(h));
            }
        }
        for (h, e) in &self.ellipses {
            if !e.construction {
                out.push(ElementHandle::Ellipse(h));
            }
        }
        for (h, p) in &self.parabolas {
            if !p.construction {
                out.push(ElementHandle::Parabola(h));
            }
        }
        for (h, hp) in &self.hyperbolas {
            if !hp.construction {
                out.push(ElementHandle::Hyperbola(h));
            }
        }
        for (h, s) in &self.splines {
            if !s.construction {
                out.push(ElementHandle::Spline(h));
            }
        }
        out
    }

    /// Return all point handles in the sketch.
    pub fn point_handles(&self) -> impl Iterator<Item = PointHandle> + '_ {
        self.points.keys()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_handles_iterator() {
        let mut sk = Sketch::new();
        let p1 = sk.add_point(0.0, 0.0);
        let p2 = sk.add_point(1.0, 0.0);
        let _lh = sk.add_line(p1, p2);
        assert_eq!(sk.line_handles().count(), 1);
    }

    #[test]
    fn is_construction_returns_correct_value() {
        let mut sk = Sketch::new();
        let p1 = sk.add_point(0.0, 0.0);
        let p2 = sk.add_point(1.0, 0.0);
        let lh = sk.add_line(p1, p2);
        assert_eq!(sk.is_construction(&ElementHandle::Line(lh)), Some(false));
        sk.set_construction(ElementHandle::Line(lh), true);
        assert_eq!(sk.is_construction(&ElementHandle::Line(lh)), Some(true));
    }

    #[test]
    fn non_construction_curves_filters_construction() {
        let mut sk = Sketch::new();
        let p1 = sk.add_point(0.0, 0.0);
        let p2 = sk.add_point(1.0, 0.0);
        let p3 = sk.add_point(1.0, 1.0);
        let l1 = sk.add_line(p1, p2);
        let l2 = sk.add_line(p2, p3);
        sk.set_construction(ElementHandle::Line(l1), true);
        let curves = sk.non_construction_curves();
        assert_eq!(curves.len(), 1);
        assert_eq!(curves[0], ElementHandle::Line(l2));
    }

    #[test]
    fn spline_point_at_evaluates() {
        let mut sk = Sketch::new();
        let p0 = sk.add_point(0.0, 0.0);
        let p1 = sk.add_point(1.0, 1.0);
        let p2 = sk.add_point(2.0, 0.0);
        let sh = sk.add_spline(&[p0, p1, p2], 2);

        // Evaluate at start and end of the knot range
        let start = sk.spline_point_at(sh, 0.0);
        assert!(start.is_some());
        let (sx, sy) = start.unwrap();
        assert!((sx - 0.0).abs() < 1e-10);
        assert!((sy - 0.0).abs() < 1e-10);

        let end = sk.spline_point_at(sh, 1.0);
        assert!(end.is_some());
        let (ex, ey) = end.unwrap();
        assert!((ex - 2.0).abs() < 1e-10);
        assert!((ey - 0.0).abs() < 1e-10);
    }
}
