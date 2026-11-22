//! The clip's boundary band: a length in space, read in a face's chart
//! by the chart's speed along each axis (tolerances are lengths;
//! chart-service design §5, operation 3).
//!
//! The band was a uv distance, a thousandth of the record's uv extent:
//! on a cylinder a millimetre across, u in radians and v in length, the
//! extent was the 2π of u and the band 6.3e-3 — wider than the face's
//! whole height of 2.5e-3 in v. A cap's ellipse on the cylinder's
//! extension, 60 tolerances past its rim, rode the boundary "everywhere"
//! and was kept whole on the face (corpus scene 1542, two millimetre
//! cylinders crossing tangentially). Every uv gap the band judges is now
//! measured in space by the speeds `|S_u|`, `|S_v|` (a unit-scale chart's
//! speeds are near one, so the band there is unchanged).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

/// Samples of the record read for its extent and the chart's speeds.
const SAMPLES: usize = 8;

/// The boundary band of one record on one face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct UvBand {
    /// The band as a length in space.
    pub(super) length: f64,
    /// The chart's mean speed along u over the record.
    ku: f64,
    /// The chart's mean speed along v over the record.
    kv: f64,
}

impl UvBand {
    /// A thousandth of the record's extent in space, at least ten
    /// tolerances (record-fit slop), read by the chart's speeds over
    /// the record.
    pub(super) fn of(brep: &BRep, face: FaceHandle, curve_uv: &NurbsCurve) -> UvBand {
        let floor = brep.tolerance.linear() * 10.0;
        let surface = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry);
        let Some(surface) = surface else {
            return UvBand {
                length: floor,
                ku: 1.0,
                kv: 1.0,
            };
        };
        let (t0, t1) = curve_uv.parameter_range();
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        let (mut ku, mut kv) = (0.0, 0.0);
        for k in 0..=SAMPLES {
            let q = curve_uv.point_at(t0 + (t1 - t0) * k as f64 / SAMPLES as f64);
            let p = surface.point_at(q.x, q.y);
            for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[i] = lo[i].min(c);
                hi[i] = hi[i].max(c);
            }
            let (su, sv) = surface.first_partials(q.x, q.y);
            ku += su.length();
            kv += sv.length();
        }
        let n = (SAMPLES + 1) as f64;
        let extent = (0..3).map(|i| (hi[i] - lo[i]).powi(2)).sum::<f64>().sqrt();
        let speed = |k: f64| if k.is_finite() && k > 0.0 { k / n } else { 1.0 };
        UvBand {
            length: (extent * 1.0e-3).max(floor),
            ku: speed(ku),
            kv: speed(kv),
        }
    }

    /// The same chart reading with the length scaled by `f`.
    pub(super) fn scaled(self, f: f64) -> UvBand {
        UvBand {
            length: self.length * f,
            ..self
        }
    }

    /// A uv point in the chart scaled to space: gaps between scaled
    /// points are lengths.
    pub(super) fn in_space(&self, u: f64, v: f64) -> (f64, f64) {
        (u * self.ku, v * self.kv)
    }

    /// The gap in space between two chart points.
    pub(super) fn gap(&self, a: (f64, f64), b: (f64, f64)) -> f64 {
        let (du, dv) = ((a.0 - b.0) * self.ku, (a.1 - b.1) * self.kv);
        (du * du + dv * dv).sqrt()
    }

    /// Whether two chart points lie within the band of each other.
    pub(super) fn holds(&self, a: (f64, f64), b: (f64, f64)) -> bool {
        self.gap(a, b) < self.length
    }

    /// The band's reach along u and along v in chart units.
    pub(super) fn per_axis(&self) -> (f64, f64) {
        (self.length / self.ku, self.length / self.kv)
    }
}
