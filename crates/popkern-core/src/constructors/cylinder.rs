//! Cylinder constructor — proper 3-face BRep.

use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, PI, TAU};

use crate::error::CoreError;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

/// Build a NURBS full-circle in the XY plane at height `z`, centered at the
/// origin, starting at `(radius, 0, z)`.
///
/// Degree 2, 9 control points (standard 4-arc representation).
/// Knot domain: `[0, TAU]`.
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
    /// Build a cylinder as a proper 3-face BRep.
    ///
    /// The cylinder has its base centered at the origin in the XY plane and
    /// extends along the +Z axis by `height`. The radius is exact (not a
    /// polygonal approximation).
    ///
    /// # Topology
    ///
    /// V=2, E=3, F=3, L=3, S=1.
    /// Euler check: V - E + F = 2 - 3 + 3 = 2.
    ///
    /// Three faces: one lateral cylindrical face and two planar cap faces.
    /// The lateral face has a 4-half-edge loop (two arcs + two seam edges).
    /// Each cap face has a single-half-edge loop (one full circle).
    ///
    /// # Geometry
    ///
    /// - Lateral face: `SurfaceGeometry::Cylinder`
    /// - Cap faces: `SurfaceGeometry::Plane`
    /// - Circle edges: degree-2 NURBS full circles
    /// - Seam edge: `CurveGeometry::Line`
    /// - All half-edges carry pcurves
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` or `height` is not positive,
    /// or if validation fails.
    pub fn cylinder(radius: f64, height: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder radius must be positive".into(),
            ));
        }
        if height <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cylinder height must be positive".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);
        let h = height;

        // -- Vertices --
        let v_bot = brep.vertices.insert(Vertex {
            position: Point3::new(radius, 0.0, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v_top = brep.vertices.insert(Vertex {
            position: Point3::new(radius, 0.0, h),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // -- Solid / Shell (placeholders, wired below) --
        let solid_h = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell_h = brep.shells.insert(Shell {
            face: Default::default(),
            solid: solid_h,
        });
        brep.solids[solid_h].shell = shell_h;

        // -- Faces & loops (placeholders) --
        let loop_lateral = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face: Default::default(),
        });
        let loop_bot = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face: Default::default(),
        });
        let loop_top = brep.loops.insert(Loop {
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
        let f_bot = brep.faces.insert(Face {
            outer_loop: loop_bot,
            inner_loops: Vec::new(),
            shell: shell_h,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let f_top = brep.faces.insert(Face {
            outer_loop: loop_top,
            inner_loops: Vec::new(),
            shell: shell_h,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });

        brep.loops[loop_lateral].face = f_lateral;
        brep.loops[loop_bot].face = f_bot;
        brep.loops[loop_top].face = f_top;
        brep.shells[shell_h].face = f_lateral;

        // -- Edges --
        let e_bot_circle = brep.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        let e_top_circle = brep.insert_edge_tracked(Edge {
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
        // Lateral face loop (4 half-edges):
        //   he_bot_arc -> he_seam_up -> he_top_arc -> he_seam_down -> he_bot_arc
        let he_bot_arc = brep.half_edges.insert(HalfEdge {
            origin: v_bot,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_bot_circle,
            loop_: loop_lateral,
            pcurve: None,
        });
        let he_seam_up = brep.half_edges.insert(HalfEdge {
            origin: v_bot,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_seam,
            loop_: loop_lateral,
            pcurve: None,
        });
        let he_top_arc = brep.half_edges.insert(HalfEdge {
            origin: v_top,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_top_circle,
            loop_: loop_lateral,
            pcurve: None,
        });
        let he_seam_down = brep.half_edges.insert(HalfEdge {
            origin: v_top,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_seam,
            loop_: loop_lateral,
            pcurve: None,
        });

        // Cap face loops (1 half-edge each, self-loop):
        let he_bot_cap = brep.half_edges.insert(HalfEdge {
            origin: v_bot,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_bot_circle,
            loop_: loop_bot,
            pcurve: None,
        });
        let he_top_cap = brep.half_edges.insert(HalfEdge {
            origin: v_top,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: e_top_circle,
            loop_: loop_top,
            pcurve: None,
        });

        // -- Wire twins --
        brep.half_edges[he_bot_arc].twin = he_bot_cap;
        brep.half_edges[he_bot_cap].twin = he_bot_arc;
        brep.half_edges[he_top_arc].twin = he_top_cap;
        brep.half_edges[he_top_cap].twin = he_top_arc;
        brep.half_edges[he_seam_up].twin = he_seam_down;
        brep.half_edges[he_seam_down].twin = he_seam_up;

        // -- Wire lateral loop: next/prev --
        brep.half_edges[he_bot_arc].next = he_seam_up;
        brep.half_edges[he_bot_arc].prev = he_seam_down;
        brep.half_edges[he_seam_up].next = he_top_arc;
        brep.half_edges[he_seam_up].prev = he_bot_arc;
        brep.half_edges[he_top_arc].next = he_seam_down;
        brep.half_edges[he_top_arc].prev = he_seam_up;
        brep.half_edges[he_seam_down].next = he_bot_arc;
        brep.half_edges[he_seam_down].prev = he_top_arc;

        // -- Wire cap self-loops --
        brep.half_edges[he_bot_cap].next = he_bot_cap;
        brep.half_edges[he_bot_cap].prev = he_bot_cap;
        brep.half_edges[he_top_cap].next = he_top_cap;
        brep.half_edges[he_top_cap].prev = he_top_cap;

        // -- Wire loop -> half_edge --
        brep.loops[loop_lateral].half_edge = he_bot_arc;
        brep.loops[loop_bot].half_edge = he_bot_cap;
        brep.loops[loop_top].half_edge = he_top_cap;

        // -- Wire edge -> half_edge --
        // CONTRACT (closed edges, see recompute_face_pcurves): the
        // edge's canonical `half_edge` must traverse the stored curve
        // FORWARD. Both circle curves run CCW (angle increasing), which
        // the lateral loop traverses forward along the BOTTOM circle
        // but BACKWARD along the TOP circle — so the top circle's
        // canonical half is the cap side. Wiring the lateral side there
        // made every pcurve re-authoring (shell/offset assembly) wind
        // the lateral and cap loops against their orientation flags.
        brep.edges[e_bot_circle].half_edge = he_bot_arc;
        brep.edges[e_top_circle].half_edge = he_top_cap;
        brep.edges[e_seam].half_edge = he_seam_up;

        // -- Wire vertex -> outgoing --
        brep.vertices[v_bot].outgoing = he_bot_arc;
        brep.vertices[v_top].outgoing = he_top_arc;

        // =====================================================================
        // Geometry
        // =====================================================================

        // -- Surface geometry --
        brep.faces[f_lateral].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Cylinder {
                origin: Point3::new(0.0, 0.0, 0.0),
                axis: Vec3::new(0.0, 0.0, 1.0),
                radius,
            },
            orientation: true,
        });
        brep.faces[f_bot].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, -1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        });
        brep.faces[f_top].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, h),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        });

        // -- Edge curves --
        // Circle edges: NURBS full circles starting at (radius, 0, z).
        brep.edges[e_bot_circle].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(xy_circle_nurbs(radius, 0.0)),
            t_start: 0.0,
            t_end: TAU,
        });
        brep.edges[e_top_circle].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(xy_circle_nurbs(radius, h)),
            t_start: 0.0,
            t_end: TAU,
        });
        // Seam edge: line from V_bot to V_top.
        brep.edges[e_seam].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(radius, 0.0, 0.0),
                direction: Vec3::new(0.0, 0.0, h),
            },
            t_start: 0.0,
            t_end: 1.0,
        });

        // -- Pcurves --
        // Lateral face pcurves (UV domain [0,TAU] x [0,h]):
        //   he_bot_arc:    (0, 0) -> (TAU, 0)
        //   he_seam_up:    (TAU, 0) -> (TAU, h)
        //   he_top_arc:    (TAU, h) -> (0, h)
        //   he_seam_down:  (0, h) -> (0, 0)
        let line_pcurve = |u0: f64, v0: f64, u1: f64, v1: f64| -> NurbsCurve {
            NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
                weights: vec![1.0, 1.0],
            }
        };

        brep.half_edges[he_bot_arc].pcurve = Some(line_pcurve(0.0, 0.0, TAU, 0.0).into());
        brep.half_edges[he_seam_up].pcurve = Some(line_pcurve(TAU, 0.0, TAU, h).into());
        brep.half_edges[he_top_arc].pcurve = Some(line_pcurve(TAU, h, 0.0, h).into());
        brep.half_edges[he_seam_down].pcurve = Some(line_pcurve(0.0, h, 0.0, 0.0).into());

        // Cap face pcurves: degree-2 rational NURBS circles.
        brep.half_edges[he_bot_cap].pcurve = Some(cap_circle_pcurve(radius).into());
        brep.half_edges[he_top_cap].pcurve = Some(cap_circle_pcurve(radius).into());

        // =====================================================================
        // Validate
        // =====================================================================
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}
