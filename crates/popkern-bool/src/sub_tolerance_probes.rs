//! Probe points inside a face, evaluated on demand (PK-1131, PK-1161).
//!
//! A face's probes are its interior point and the stations of a grid
//! over its parameter box that lie inside its trim, or the middle of
//! its boundary samples when it is too narrow to hold any. One point a
//! face is not enough: a cube less its inscribed ball has its faces'
//! middles where the ball touches them, each within the tolerance of
//! the ball, though the corners between are far from it.
//!
//! The sub-tolerance tests run on every boolean, and most of them end
//! at the first probe (an ordinary body's first face is far from the
//! rest) or skip most stations (a station far from the other face's
//! box). The grid is laid out only when the interior point does not
//! end the search, and a station's trim containment, the costly part,
//! is tested only for a station the caller keeps.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

/// Grid stations per parameter direction probed inside each face.
const GRID: usize = 4;

/// A grid station: its point, and its trim containment once tested.
struct Station {
    point: Point3,
    uv: (f64, f64),
    inside: Option<bool>,
}

/// The probe points of one face, in a fixed order: the interior point,
/// the stations inside the trim, else the boundary samples' middle.
pub(super) struct FaceProbes<'a> {
    brep: &'a BRep,
    face: FaceHandle,
    interior: Option<Option<Point3>>,
    stations: Option<Vec<Station>>,
}

impl<'a> FaceProbes<'a> {
    /// The probes of `face`, none evaluated yet.
    pub(super) fn new(brep: &'a BRep, face: FaceHandle) -> Self {
        Self {
            brep,
            face,
            interior: None,
            stations: None,
        }
    }

    fn interior(&mut self) -> Option<Point3> {
        let (brep, face) = (self.brep, self.face);
        *self.interior.get_or_insert_with(|| {
            crate::entity_lineage::interior_point(brep, face).map(|(p, _)| p)
        })
    }

    fn stations(&mut self) -> &mut Vec<Station> {
        let (brep, face) = (self.brep, self.face);
        self.stations.get_or_insert_with(|| {
            let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                return Vec::new();
            };
            let s = &surface.geometry;
            let samples = brep.face_outer_loop_sample_positions(face, 4);
            let ((u0, u1), (v0, v1)) =
                popkern_core::geometry::uv_projection::exact_bounded_parameter_range(&samples, s);
            let mut out = Vec::with_capacity(GRID * GRID);
            for i in 0..GRID {
                for j in 0..GRID {
                    let u = u0 + (u1 - u0) * (i as f64 + 0.5) / GRID as f64;
                    let v = v0 + (v1 - v0) * (j as f64 + 0.5) / GRID as f64;
                    out.push(Station {
                        point: s.point_at(u, v),
                        uv: (u, v),
                        inside: None,
                    });
                }
            }
            out
        })
    }

    /// Whether station `k` lies inside the face's trim (tested once).
    fn station_inside(&mut self, k: usize) -> bool {
        let (brep, face) = (self.brep, self.face);
        let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            return false;
        };
        let station = &mut self.stations()[k];
        *station.inside.get_or_insert_with(|| {
            brep.face_uv_trim_containment(face, &surface.geometry, station.uv.0, station.uv.1)
                == Some(FaceUvContainment::Interior)
        })
    }

    /// The middle of the face's boundary samples.
    fn fallback(&self) -> Option<Point3> {
        let samples = self.brep.face_outer_loop_sample_positions(self.face, 4);
        if samples.is_empty() {
            return None;
        }
        let n = samples.len() as f64;
        let sum = samples.iter().fold(Point3::new(0.0, 0.0, 0.0), |acc, q| {
            Point3::new(acc.x + q.x, acc.y + q.y, acc.z + q.z)
        });
        Some(Point3::new(sum.x / n, sum.y / n, sum.z / n))
    }

    /// Whether the face has no probe at all.
    pub(super) fn is_empty(&mut self) -> bool {
        if self.interior().is_some() {
            return false;
        }
        let n = self.stations().len();
        !(0..n).any(|k| self.station_inside(k)) && self.fallback().is_none()
    }

    /// The first `f(p)` that is `Some` over the probes `p` that `keep`
    /// accepts, in the probes' order. `keep` is asked before a
    /// station's trim containment is tested, so a cheap rejection
    /// spares it.
    pub(super) fn find_map<T>(
        &mut self,
        keep: impl Fn(&Point3) -> bool,
        mut f: impl FnMut(Point3) -> Option<T>,
    ) -> Option<T> {
        let interior = self.interior();
        if let Some(p) = interior
            && keep(&p)
            && let Some(t) = f(p)
        {
            return Some(t);
        }
        let n = self.stations().len();
        let mut any_inside = false;
        for k in 0..n {
            let p = self.stations()[k].point;
            // The fallback stands only when no station is inside, so
            // with no interior point every station's containment is
            // needed; otherwise a rejected station is never tested.
            if !keep(&p) && interior.is_some() {
                continue;
            }
            if !self.station_inside(k) {
                continue;
            }
            any_inside = true;
            if keep(&p)
                && let Some(t) = f(p)
            {
                return Some(t);
            }
        }
        if interior.is_none() && !any_inside {
            return self.fallback().filter(|p| keep(p)).and_then(f);
        }
        None
    }
}
