//! Exact sections of the built document (PK-856): the page's clipping
//! plane cut through every body by the kernel (`plane_section`), drawn
//! as SVG and DXF and measured per body. The clipping view itself stays
//! display-only; this is what "export this section" asks for.

use popkern_assembly::Document;
use popkern_bool::ops::plane_section::{SectionPlane, plane_section};
use popkern_core::geometry::{Point3, Vec3};
use popkern_io::section_drawing::{SectionArea, SectionBody, SectionDrawing};
use serde::Serialize;

use crate::session::Session;
use popkern_rhai::Outcome;

/// One body's measured section.
#[derive(Debug, Clone, Serialize)]
pub struct SectionMeasure {
    /// The body (`body N`, or an occurrence path).
    pub body: String,
    /// Its section's area, document units squared.
    pub area: f64,
    /// Its centroid in the plane's 2-D frame, document units.
    pub centroid: [f64; 2],
}

/// A section's reply.
#[derive(Debug, Clone, Serialize)]
pub struct SectionReply {
    /// The drawing as SVG (document units, y up in the plane's frame).
    pub svg: String,
    /// The drawing as ASCII DXF.
    pub dxf: String,
    /// Each body with a section, measured.
    pub bodies: Vec<SectionMeasure>,
    /// The plane's 2-D frame: x and y directions in the document.
    pub frame: [[f64; 3]; 2],
}

impl Session {
    /// Cut the built document by the page's clipping plane `{ normal,
    /// offset }` (world SI, the side `normal · x > offset` removed, as
    /// the viewer draws it): every body's exact section, one drawing.
    ///
    /// # Errors
    ///
    /// [`Outcome::GateFailure`] with nothing built, a zero normal, or a
    /// body the section boolean declines.
    pub fn section(&mut self, normal: [f64; 3], offset: f64) -> Result<SectionReply, Outcome> {
        let gate = |detail: String| Outcome::GateFailure {
            step: "section".into(),
            detail,
        };
        let Some(document) = self.document_for_section() else {
            return Err(gate(
                "build the model first: a section cuts the built document".into(),
            ));
        };
        let unit = popkern_rhai::film::FilmContext::from_document(document).unit_scale;
        let n = Vec3::new(normal[0], normal[1], normal[2]);
        let Some(unit_n) = n.normalized() else {
            return Err(gate("a section plane needs a non-zero normal".into()));
        };
        let origin = Point3::new(0.0, 0.0, 0.0) + unit_n * (offset / unit);
        let plane = SectionPlane::new(origin, unit_n)
            .ok_or_else(|| gate("a section plane needs a non-zero normal".into()))?;
        let bodies = bodies_of(document).map_err(gate)?;
        let mut drawing = SectionDrawing::default();
        let mut measures = Vec::new();
        for (name, brep) in bodies {
            // Loops drawn within a thousandth of the body's size.
            let size =
                popkern_core::geometry::intersection::surface_subdivision::Aabb::from_points(
                    &brep
                        .vertices
                        .values()
                        .map(|v| v.position)
                        .collect::<Vec<_>>(),
                )
                .map_or(1.0, |b| b.min.distance_to(&b.max));
            let chord = 1e-3 * size.max(1e-9);
            let section = plane_section(&brep, &plane, chord, &brep.tolerance)
                .map_err(|e| gate(format!("{name}: the section declined: {e}")))?;
            if section.regions.is_empty() {
                continue;
            }
            let (mut area, mut cx, mut cy) = (0.0, 0.0, 0.0);
            for r in &section.regions {
                area += r.area;
                cx += r.area * r.centroid[0];
                cy += r.area * r.centroid[1];
            }
            measures.push(SectionMeasure {
                body: name.clone(),
                area,
                centroid: if area > 0.0 {
                    [cx / area, cy / area]
                } else {
                    [0.0, 0.0]
                },
            });
            drawing.bodies.push(SectionBody {
                name,
                regions: section
                    .regions
                    .into_iter()
                    .map(|r| SectionArea {
                        outer: r.outer.points,
                        holes: r.holes.into_iter().map(|h| h.points).collect(),
                    })
                    .collect(),
            });
        }
        let (u, v) = (plane.u_axis, plane.v_axis());
        Ok(SectionReply {
            svg: drawing.to_svg(),
            dxf: drawing.to_dxf(),
            bodies: measures,
            frame: [[u.x, u.y, u.z], [v.x, v.y, v.z]],
        })
    }
}

/// Every alive body of `document` in world placement, named.
fn bodies_of(
    document: &mut Document,
) -> Result<Vec<(String, popkern_core::topology::BRep)>, String> {
    match document {
        Document::Part(model) => {
            let mut out = Vec::new();
            for (i, body) in model.alive_bodies().into_iter().enumerate() {
                let brep = model.brep_of(body).map_err(|e| e.to_string())?;
                out.push((format!("body {i}"), brep));
            }
            Ok(out)
        }
        Document::Assembly(assembly) => {
            let flat = assembly.flatten().map_err(|e| e.to_string())?;
            let multi = |path: &str| flat.iter().filter(|f| f.path.to_string() == path).count() > 1;
            Ok(flat
                .iter()
                .filter(|f| f.visible)
                .map(|f| {
                    let path = f.path.to_string();
                    let name = if multi(&path) {
                        format!("{path} body {}", f.body)
                    } else {
                        path
                    };
                    (name, f.brep.clone())
                })
                .collect())
        }
    }
}
