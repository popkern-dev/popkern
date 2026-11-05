//! Polygon frustum constructor: two parallel polygons with the same
//! number of corners, joined corner to corner by planar side faces (a
//! pyramid's frustum, a prism, a wedge's ramp).

use crate::error::CoreError;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// The Newell normal of a polygon (its area-weighted normal, any
/// winding).
fn newell(points: &[Point3]) -> Vec3 {
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    for (i, p) in points.iter().enumerate() {
        let q = points[(i + 1) % points.len()];
        n = n + Vec3::new(
            (p.y - q.y) * (p.z + q.z),
            (p.z - q.z) * (p.x + q.x),
            (p.x - q.x) * (p.y + q.y),
        );
    }
    n
}

fn centroid(points: &[Point3]) -> Point3 {
    let s = points.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| {
        acc + (*p - Point3::new(0.0, 0.0, 0.0))
    });
    Point3::new(0.0, 0.0, 0.0) + s * (1.0 / points.len() as f64)
}

impl BRep {
    /// Build the solid between a `bottom` and a `top` polygon, corner
    /// `i` of one joined to corner `i` of the other.
    ///
    /// Built as an Euler extrude of the bottom face whose cap corners
    /// are then moved onto the top polygon, so the topology is a
    /// prism's (V = 2n, E = 3n, F = n + 2); every face is then given
    /// its plane. Either winding is accepted: the corners are taken in
    /// the order whose normal points from the bottom to the top, so the
    /// shell faces outward.
    ///
    /// # Errors
    ///
    /// [`CoreError::PreconditionFailed`] when the polygons differ in
    /// corner count or have fewer than 3, when the top does not stand
    /// off the bottom's plane, or when a side face (or either polygon)
    /// is not planar within `tolerance` (each is checked, so a twisted
    /// side is reported rather than served with a fitted plane);
    /// [`CoreError::Validation`] if the result fails validation.
    pub fn polygon_frustum(
        bottom: &[Point3],
        top: &[Point3],
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        let n = bottom.len();
        if n < 3 || top.len() != n {
            return Err(CoreError::PreconditionFailed(format!(
                "polygon_frustum needs two polygons of the same corner count (at least 3), got {} and {}",
                bottom.len(),
                top.len()
            )));
        }
        let rise = centroid(top) - centroid(bottom);
        let normal = newell(bottom);
        let along = normal.dot(&rise);
        if along.abs() <= tolerance.linear() * normal.length() {
            return Err(CoreError::PreconditionFailed(
                "polygon_frustum: the top polygon does not stand off the bottom's plane".into(),
            ));
        }
        let (bottom, top): (Vec<Point3>, Vec<Point3>) = if along > 0.0 {
            (bottom.to_vec(), top.to_vec())
        } else {
            (
                bottom.iter().rev().copied().collect(),
                top.iter().rev().copied().collect(),
            )
        };
        let planar = |quad: &[Point3], what: &str| -> Result<(), CoreError> {
            let Some(unit) = newell(quad).normalized() else {
                return Err(CoreError::PreconditionFailed(format!(
                    "polygon_frustum: the {what} is degenerate"
                )));
            };
            let c = centroid(quad);
            match quad
                .iter()
                .map(|p| (*p - c).dot(&unit).abs())
                .fold(0.0, f64::max)
            {
                d if d <= tolerance.linear() => Ok(()),
                d => Err(CoreError::PreconditionFailed(format!(
                    "polygon_frustum: the {what} is not planar (a corner {d:e} off its plane)"
                ))),
            }
        };
        planar(&bottom, "bottom")?;
        planar(&top, "top")?;
        for i in 0..n {
            let j = (i + 1) % n;
            planar(
                &[bottom[i], bottom[j], top[j], top[i]],
                &format!("side face {i}"),
            )?;
        }

        let mut brep = BRep::new(tolerance);
        let mvfs = brep.mvfs(bottom[0])?;
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, bottom[1])?;
        let mut last_twin = mev1.half_edge_twin;
        for &point in &bottom[2..] {
            let mev = brep.mev(last_twin, point)?;
            last_twin = mev.half_edge_twin;
        }
        let mef = brep.mef(last_twin, mev1.half_edge)?;
        let unit = rise
            .normalized()
            .ok_or_else(|| CoreError::PreconditionFailed("polygon_frustum: no rise".into()))?;
        let distance = rise.length();
        let extruded = brep.extrude(mef.face, unit, distance, None)?;
        // Each cap corner sits above its bottom corner by the rise;
        // move it onto its top corner.
        for v in extruded.new_vertices {
            let p = brep.vertices[v].position;
            let Some(i) = bottom
                .iter()
                .position(|b| (*b + unit * distance).distance_to(&p) <= tolerance.linear())
            else {
                return Err(CoreError::PreconditionFailed(
                    "polygon_frustum: an extruded corner matches no bottom corner".into(),
                ));
            };
            brep.vertices[v].position = top[i];
        }
        brep.attach_geometry()?;
        brep.attach_planar_pcurves()?;
        brep.orient_planar_outward()?;
        brep.validate().map_err(CoreError::Validation)?;
        Ok(brep)
    }

    /// Turn a closed body of planar faces outward: when its faces point
    /// into the material (the Euler constructors' prisms, cubes and
    /// pyramids are built that way), reverse the whole body with
    /// [`BRep::reverse_body_orientation`], winding and flags together.
    /// Returns whether it reversed.
    ///
    /// Decided by the enclosed volume Σ area · (n · corner) / 3 over the
    /// faces, n each face's effective normal (surface normal by its
    /// orientation flag): positive for an outward shell, negative for
    /// an inward one, exactly, whatever the body's position.
    ///
    /// # Errors
    ///
    /// [`CoreError::PreconditionFailed`] if a face is not planar, or the
    /// enclosed volume is zero (no material side to decide).
    pub fn orient_planar_outward(&mut self) -> Result<bool, CoreError> {
        use crate::geometry::surface::SurfaceGeometry;
        let mut volume = 0.0;
        for (fh, face) in &self.faces {
            let Some(SurfaceGeometry::Plane { normal, .. }) =
                face.surface.as_ref().map(|s| &s.geometry)
            else {
                return Err(CoreError::PreconditionFailed(
                    "orient_planar_outward: every face must be planar".into(),
                ));
            };
            let flag = face.surface.as_ref().is_some_and(|s| s.orientation);
            let outward = if flag { *normal } else { -*normal };
            let corners = self.face_outer_loop_positions(fh);
            let Some(first) = corners.first() else {
                continue;
            };
            let area = newell(&corners).length() * 0.5;
            volume += area * outward.dot(&(*first - Point3::new(0.0, 0.0, 0.0))) / 3.0;
        }
        if volume == 0.0 || !volume.is_finite() {
            return Err(CoreError::PreconditionFailed(
                "orient_planar_outward: the body encloses no volume".into(),
            ));
        }
        if volume < 0.0 {
            self.reverse_body_orientation();
            return Ok(true);
        }
        Ok(false)
    }
}

#[cfg(test)]
#[path = "frustum_tests.rs"]
mod tests;
