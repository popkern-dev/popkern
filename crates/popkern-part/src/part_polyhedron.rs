//! Polyhedral features on [`Part`]: hull construction, vertex
//! truncation, and the representation-level edit escape hatch.
//! Algorithms live in [`crate::poly`]; this file is the Part-level
//! wiring (validation, feature log threading).

use popkern_core::tolerance::Tolerance;

use crate::error::PartError;
use crate::part::Part;

impl Part {
    /// Build a part as the convex hull of `points` (coplanar hull
    /// triangles merged into polygons). Intended for canonical
    /// solids — the hull is O(n⁴) plane enumeration.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] for degenerate input;
    /// validation errors if wiring fails.
    pub fn polyhedron(
        points: &[popkern_core::geometry::Point3],
        tolerance: Tolerance,
    ) -> Result<Part, PartError> {
        let faces = crate::poly::convex_hull_faces(points);
        if faces.len() < 4 {
            return Err(PartError::NotApplicable(format!(
                "convex hull of {} points produced only {} faces",
                points.len(),
                faces.len()
            )));
        }
        let specs: Vec<crate::poly::PolyFace> = faces
            .into_iter()
            .map(crate::poly::PolyFace::simple)
            .collect();
        let brep = crate::poly::build_polyhedron_brep(points, &specs, tolerance)?;
        Part::from_validated(
            brep,
            tolerance,
            "polyhedron",
            format!("convex hull of {} points", points.len()),
        )
    }

    /// Uniform vertex truncation on a polyhedral part: every vertex
    /// becomes a flat corner facet cut `fraction` (0 < f < 0.5)
    /// along each incident edge.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the part is not polyhedral
    /// (non-planar faces or curved edges) or the fraction is out of
    /// range.
    pub fn truncate_vertices(&self, fraction: f64) -> Result<Part, PartError> {
        if !(fraction > 0.0 && fraction < 0.5) {
            return Err(PartError::NotApplicable(format!(
                "truncation fraction must be in (0, 0.5), got {fraction}"
            )));
        }
        let (points, faces) = crate::poly::extract_polyhedron(self.brep())?;
        if faces.iter().any(|f| !f.inners.is_empty()) {
            return Err(PartError::NotApplicable(
                "truncate_vertices requires faces without inner loops (truncate before engraving)"
                    .into(),
            ));
        }
        let outers: Vec<Vec<usize>> = faces.iter().map(|f| f.outer.clone()).collect();
        let (tp, tf) = crate::poly::truncate_vertices(&points, &outers, fraction);
        let specs: Vec<crate::poly::PolyFace> =
            tf.into_iter().map(crate::poly::PolyFace::simple).collect();
        let brep = crate::poly::build_polyhedron_brep(&tp, &specs, self.tolerance())?;
        let mut part = Part::from_validated(
            brep,
            self.tolerance(),
            "truncate_vertices",
            format!("fraction {fraction}"),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Edit the part as a polyhedron representation: the closure
    /// receives `(points, faces)` (see [`crate::poly`]) and may
    /// engrave pockets or otherwise extend them; the result is
    /// rebuilt and validated. The escape hatch for per-face work
    /// like dice numbering.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the part is not polyhedral;
    /// whatever the closure returns; validation errors on rebuild.
    pub fn edit_polyhedron(
        &self,
        summary: impl Into<String>,
        edit: impl FnOnce(
            &mut Vec<popkern_core::geometry::Point3>,
            &mut Vec<crate::poly::PolyFace>,
        ) -> Result<(), PartError>,
    ) -> Result<Part, PartError> {
        let (mut points, mut faces) = crate::poly::extract_polyhedron(self.brep())?;
        edit(&mut points, &mut faces)?;
        let brep = crate::poly::build_polyhedron_brep(&points, &faces, self.tolerance())?;
        let mut part =
            Part::from_validated(brep, self.tolerance(), "edit_polyhedron", summary.into())?;
        part.merge_history_from(self);
        Ok(part)
    }
}
