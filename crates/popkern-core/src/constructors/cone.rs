//! Cone constructor — proper 2-face BRep.
use popkern_math::PortableMath;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI, TAU};

use crate::error::CoreError;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

/// Build a degree-2 NURBS full-circle in the XY plane at height `z`.
fn xy_circle_nurbs(radius: f64, z: f64) -> NurbsCurve {
    let r = radius;
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0,
            0.0,
            0.0,
            FRAC_PI_2,
            FRAC_PI_2,
            PI,
            PI,
            PI + FRAC_PI_2,
            PI + FRAC_PI_2,
            TAU,
            TAU,
            TAU,
        ],
        control_points: vec![
            Point3::new(r, 0.0, z),
            Point3::new(r, r, z),
            Point3::new(0.0, r, z),
            Point3::new(-r, r, z),
            Point3::new(-r, 0.0, z),
            Point3::new(-r, -r, z),
            Point3::new(0.0, -r, z),
            Point3::new(r, -r, z),
            Point3::new(r, 0.0, z),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Build a degree-2 rational NURBS circle pcurve in the cap's UV plane.
///
/// The cap face's plane has `u_axis=(1,0,0), v_axis=(0,1,0)`, so UV maps
/// directly to XY. The pcurve is a circle of `radius` at the UV origin.
fn cap_circle_pcurve(radius: f64) -> NurbsCurve {
    let r = radius;
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, r, 0.0),
            Point3::new(0.0, r, 0.0),
            Point3::new(-r, r, 0.0),
            Point3::new(-r, 0.0, 0.0),
            Point3::new(-r, -r, 0.0),
            Point3::new(0.0, -r, 0.0),
            Point3::new(r, -r, 0.0),
            Point3::new(r, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

impl BRep {
    /// Build a cone as a proper 2-face BRep.
    ///
    /// The cone has its base centered at the origin in the XY plane and its
    /// apex at `(0, 0, height)`. The base circle has the given `radius`.
    ///
    /// # Topology
    ///
    /// V=2, E=2, F=2, L=2, S=1.
    /// Euler check: V - E + F = 2 - 2 + 2 = 2.
    ///
    /// Two faces: one lateral conical face and one planar base cap.
    /// The lateral face has a 3-half-edge loop (one base arc + two seam edges).
    /// The seam edge's two half-edges are both on the lateral face (they are
    /// twins of each other). The cap face has a single-half-edge self-loop.
    ///
    /// # Geometry
    ///
    /// - Lateral face: `SurfaceGeometry::Cone`
    /// - Base cap: `SurfaceGeometry::Plane`
    /// - Base circle edge: degree-2 NURBS full circle
    /// - Seam edge: `CurveGeometry::Line`
    /// - All half-edges carry pcurves
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` or `height` is not positive,
    /// or if validation fails.
    pub fn cone(radius: f64, height: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cone radius must be positive".into(),
            ));
        }
        if height <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cone height must be positive".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);
        let h = height;

        // -- Vertices --
        let v_base = brep.vertices.insert(Vertex {
            position: Point3::new(radius, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v_apex = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, h),
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

        // -- Faces & loops --
        let loop_lateral = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face: Default::default(),
        });
        let loop_cap = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face: Default::default(),
        });

        let f_lateral = brep.faces.insert(Face {
            outer_loop: loop_lateral,
            inner_loops: Vec::new(),
            shell: shell_h,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let f_cap = brep.faces.insert(Face {
            outer_loop: loop_cap,
            inner_loops: Vec::new(),
            shell: shell_h,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });

        brep.loops[loop_lateral].face = f_lateral;
        brep.loops[loop_cap].face = f_cap;
        brep.shells[shell_h].face = f_lateral;

        // -- Edges --
        let e_base_circle = brep.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        let e_seam = brep.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        // -- Half-edges --
        // Lateral face loop (3 half-edges):
        //   he_base_arc -> he_seam_up -> he_seam_down -> he_base_arc
        let he_base_arc = brep.half_edges.insert(HalfEdge {
            origin: v_base,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_base_circle,
            loop_: loop_lateral,
            pcurve: None,
        });
        let he_seam_up = brep.half_edges.insert(HalfEdge {
            origin: v_base,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_seam,
            loop_: loop_lateral,
            pcurve: None,
        });
        let he_seam_down = brep.half_edges.insert(HalfEdge {
            origin: v_apex,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_seam,
            loop_: loop_lateral,
            pcurve: None,
        });

        // Cap face loop (1 half-edge, self-loop):
        let he_base_cap = brep.half_edges.insert(HalfEdge {
            origin: v_base,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_base_circle,
            loop_: loop_cap,
            pcurve: None,
        });

        // -- Wire twins --
        brep.half_edges[he_base_arc].twin = he_base_cap;
        brep.half_edges[he_base_cap].twin = he_base_arc;
        brep.half_edges[he_seam_up].twin = he_seam_down;
        brep.half_edges[he_seam_down].twin = he_seam_up;

        // -- Wire lateral loop: next/prev --
        brep.half_edges[he_base_arc].next = he_seam_up;
        brep.half_edges[he_base_arc].prev = he_seam_down;
        brep.half_edges[he_seam_up].next = he_seam_down;
        brep.half_edges[he_seam_up].prev = he_base_arc;
        brep.half_edges[he_seam_down].next = he_base_arc;
        brep.half_edges[he_seam_down].prev = he_seam_up;

        // -- Wire cap self-loop --
        brep.half_edges[he_base_cap].next = he_base_cap;
        brep.half_edges[he_base_cap].prev = he_base_cap;

        // -- Wire loop -> half_edge --
        brep.loops[loop_lateral].half_edge = he_base_arc;
        brep.loops[loop_cap].half_edge = he_base_cap;

        // -- Wire edge -> half_edge --
        // CONTRACT (closed edges, see recompute_face_pcurves): the
        // canonical `half_edge` must traverse the stored curve FORWARD.
        // The base circle curve runs CCW, but the lateral chart's u
        // direction is CW here (axis (0,0,-1) puts perp2 at -ŷ), so the
        // lateral loop traverses the base circle BACKWARD — the cap
        // side is the forward (canonical) half.
        brep.edges[e_base_circle].half_edge = he_base_cap;
        brep.edges[e_seam].half_edge = he_seam_up;

        // -- Wire vertex -> outgoing --
        brep.vertices[v_base].outgoing = he_base_arc;
        brep.vertices[v_apex].outgoing = he_seam_down;

        // =====================================================================
        // Geometry
        // =====================================================================

        // -- Surface geometry --
        let half_angle = (radius / height).patan();
        brep.faces[f_lateral].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Cone {
                apex: Point3::new(0.0, 0.0, h),
                axis: Vec3::new(0.0, 0.0, -1.0),
                half_angle,
            },
            orientation: true,
        });
        brep.faces[f_cap].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, -1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        });

        // -- Edge curves --
        brep.edges[e_base_circle].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(xy_circle_nurbs(radius, 0.0)),
            t_start: 0.0,
            t_end: TAU,
        });
        brep.edges[e_seam].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(radius, 0.0, 0.0),
                direction: Vec3::new(-radius, 0.0, h),
            },
            t_start: 0.0,
            t_end: 1.0,
        });

        // -- Pcurves --
        // Lateral face pcurves (triangular UV domain):
        //   he_base_arc:   (0, height) -> (TAU, height)
        //   he_seam_up:    (TAU, height) -> (TAU, 0)   (apex at v=0)
        //   he_seam_down:  (0, 0) -> (0, height)       (apex at v=0)
        let line_pcurve = |u0: f64, v0: f64, u1: f64, v1: f64| -> NurbsCurve {
            NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
                weights: vec![1.0, 1.0],
            }
        };

        brep.half_edges[he_base_arc].pcurve = Some(line_pcurve(0.0, h, TAU, h).into());
        brep.half_edges[he_seam_up].pcurve = Some(line_pcurve(TAU, h, TAU, 0.0).into());
        brep.half_edges[he_seam_down].pcurve = Some(line_pcurve(0.0, 0.0, 0.0, h).into());

        // Cap face pcurve: degree-2 rational NURBS circle.
        brep.half_edges[he_base_cap].pcurve = Some(cap_circle_pcurve(radius).into());

        // =====================================================================
        // Validate
        // =====================================================================
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
#[path = "cone_tests.rs"]
mod tests;
