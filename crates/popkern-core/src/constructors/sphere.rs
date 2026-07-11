//! Sphere constructor — proper 1-face BRep.

use std::f64::consts::{PI, TAU};

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

impl BRep {
    /// Build a sphere as a proper 1-face BRep.
    ///
    /// The sphere is centered at the origin with the given `radius`.
    /// The two vertices (poles) are at (0, 0, -radius) (south) and
    /// (0, 0, radius) (north). A single seam edge runs between the poles
    /// along the prime meridian (x > 0, y = 0).
    ///
    /// # Topology
    ///
    /// V=2, E=1, F=1, L=1, S=1.
    /// Euler check: V - E + F = 2 - 1 + 1 = 2.
    ///
    /// The single face has two half-edges (seam up and seam down) that are
    /// twins of each other, both on the same face. The loop traversal is:
    /// he_seam_up -> he_seam_down -> he_seam_up.
    ///
    /// # Geometry
    ///
    /// - Face: `SurfaceGeometry::Sphere`
    /// - Seam edge: the prime-meridian semicircle (x ≥ 0, y = 0) from
    ///   south to north pole, as an exact rational quadratic NURBS arc
    /// - Pcurves define the UV mapping on the rectangular domain
    ///   `[0, 2pi] x [0, pi]`
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` is not positive,
    /// or if validation fails.
    pub fn sphere(radius: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "sphere radius must be positive".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);

        // -- Vertices --
        let v_south = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, -radius),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v_north = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, radius),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // -- Solid / Shell --
        let solid_h = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell_h = brep.shells.insert(Shell {
            face: Default::default(),
            solid: solid_h,
        });
        brep.solids[solid_h].shell = shell_h;

        // -- Face & loop --
        let loop_face = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face: Default::default(),
        });

        let f_sphere = brep.faces.insert(Face {
            outer_loop: loop_face,
            inner_loops: Vec::new(),
            shell: shell_h,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });

        brep.loops[loop_face].face = f_sphere;
        brep.shells[shell_h].face = f_sphere;

        // -- Edge --
        let e_seam = brep.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        // -- Half-edges --
        // he_seam_up: south -> north (right seam, u=TAU)
        // he_seam_down: north -> south (left seam, u=0)
        let he_seam_up = brep.half_edges.insert(HalfEdge {
            origin: v_south,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_seam,
            loop_: loop_face,
            pcurve: None,
        });
        let he_seam_down = brep.half_edges.insert(HalfEdge {
            origin: v_north,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_seam,
            loop_: loop_face,
            pcurve: None,
        });

        // -- Wire twins --
        brep.half_edges[he_seam_up].twin = he_seam_down;
        brep.half_edges[he_seam_down].twin = he_seam_up;

        // -- Wire loop: next/prev --
        brep.half_edges[he_seam_up].next = he_seam_down;
        brep.half_edges[he_seam_up].prev = he_seam_down;
        brep.half_edges[he_seam_down].next = he_seam_up;
        brep.half_edges[he_seam_down].prev = he_seam_up;

        // -- Wire loop -> half_edge --
        brep.loops[loop_face].half_edge = he_seam_up;

        // -- Wire edge -> half_edge --
        brep.edges[e_seam].half_edge = he_seam_up;

        // -- Wire vertex -> outgoing --
        brep.vertices[v_south].outgoing = he_seam_up;
        brep.vertices[v_north].outgoing = he_seam_down;

        // =====================================================================
        // Geometry
        // =====================================================================

        // -- Surface geometry --
        brep.faces[f_sphere].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius,
            },
            orientation: true,
        });

        // -- Edge curve --
        // The prime-meridian semicircle from south to north pole
        // (x ≥ 0, y = 0), as an exact rational quadratic arc. The
        // historical pole-to-pole LINE was a chord through the
        // sphere's INTERIOR: every consumer trusting edge geometry
        // saw the chord (measured: the edge-face intersector reported
        // phantom seam×face hits at a cube's face centers — points
        // nowhere near the sphere — split the graph seam there, and
        // the boolean's sphere face lost its true seam borders).
        let w = std::f64::consts::FRAC_1_SQRT_2;
        brep.edges[e_seam].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(NurbsCurve {
                degree: 2,
                knots: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
                control_points: vec![
                    Point3::new(0.0, 0.0, -radius),
                    Point3::new(radius, 0.0, -radius),
                    Point3::new(radius, 0.0, 0.0),
                    Point3::new(radius, 0.0, radius),
                    Point3::new(0.0, 0.0, radius),
                ],
                weights: vec![1.0, w, 1.0, w, 1.0],
            }),
            t_start: 0.0,
            t_end: 1.0,
        });

        // -- Pcurves --
        // UV domain: [0, TAU] x [0, PI]
        // he_seam_up (south->north, right seam): (TAU, PI) -> (TAU, 0)
        // he_seam_down (north->south, left seam): (0, 0) -> (0, PI)
        let line_pcurve = |u0: f64, v0: f64, u1: f64, v1: f64| -> NurbsCurve {
            NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
                weights: vec![1.0, 1.0],
            }
        };

        brep.half_edges[he_seam_up].pcurve = Some(line_pcurve(TAU, PI, TAU, 0.0).into());
        brep.half_edges[he_seam_down].pcurve = Some(line_pcurve(0.0, 0.0, 0.0, PI).into());

        // =====================================================================
        // Validate
        // =====================================================================
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
