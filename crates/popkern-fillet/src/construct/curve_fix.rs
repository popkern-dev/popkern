//! Edge-curve repair for cylindrical fillet ribbon faces.
//!
//! After topology operations and surface assignments, a cylinder
//! fillet face's outer loop can contain edges whose 3-D curves are
//! straight lines spanning the cylinder cross-section as *chords*
//! rather than arcs on the cylinder surface. OCCT-based slicers
//! (Bambu, FreeCAD, …) read the resulting STEP, find the boundary
//! lying off the cylindrical surface, and refuse to render the
//! face — so the user sees a sharp edge where there should be a
//! fillet ribbon.
//!
//! ## Construction-site complication
//!
//! The boundary vertex *positions* in the BRep also drift slightly
//! from the analytic cylinder surface (typically a few percent of
//! radius). Causes upstream:
//! - corner-blend insertion picks vertex positions that are a
//!   compromise between multiple incident surfaces (cylinder ×
//!   sphere/Nurbs cap × planar face);
//! - `surface_fix::correct_face_surfaces` re-fits surfaces but
//!   doesn't move vertices.
//!
//! Repairing the chord by replacing it with an *exact* analytic
//! arc on the cylinder fails: the arc's endpoints would lie at the
//! cylinder cross-section, not at the actual vertex positions, and
//! the resulting `ShellGap` validation error is bigger than the
//! original chord error.
//!
//! ## Strategy
//!
//! Build a **rational quadratic NURBS arc through 3 points** —
//! `(p_start, p_mid_on_cylinder, p_end)`. The middle point lives on
//! the analytic cylinder cross-section; the endpoints are the
//! actual face vertices. The unique circle through those 3 points
//! is *approximately* on the cylinder (close in the middle, exact
//! at the endpoints) and is representable as a degree-2 rational
//! NURBS. STEP exports this as a B_SPLINE_CURVE_WITH_KNOTS that
//! OCCT readers tessellate as a curved boundary, which is what the
//! slicer cares about.
//!
//! Edges whose endpoints span different cross-sections (helical
//! relationship on the cylinder) are skipped — they need a more
//! general fit and aren't produced by the present fillet pipeline.
//! Axial edges (same circumferential angle) keep their straight
//! line — the line lies on the cylinder.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Walk every face whose surface is cylindrical and replace
/// chord-line boundary edges with rational-quadratic arcs through
/// the existing vertex positions, bulging toward the cylinder
/// cross-section.
pub(crate) fn fix_cylinder_face_edges(brep: &mut BRep) {
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in face_handles {
        let Some((origin, axis, radius)) = read_cylinder(brep, fh) else {
            continue;
        };
        let Some(loop_h) = brep.faces.get(fh).map(|f| f.outer_loop) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(loop_h) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        if !brep.half_edges.contains_key(start_he) {
            continue;
        }
        let mut visited: Vec<EdgeHandle> = Vec::new();
        let mut cursor = start_he;
        let mut safety = 64;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            let he = &brep.half_edges[cursor];
            visited.push(he.edge);
            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
        for edge in visited {
            try_rewrite_chord_to_arc(brep, edge, origin, axis, radius);
        }
    }
}

fn read_cylinder(brep: &BRep, fh: FaceHandle) -> Option<(Point3, Vec3, f64)> {
    let face = brep.faces.get(fh)?;
    let sof = face.surface.as_ref()?;
    match &sof.geometry {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some((*origin, *axis, *radius)),
        _ => None,
    }
}

fn try_rewrite_chord_to_arc(
    brep: &mut BRep,
    edge: EdgeHandle,
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    // Only rewrite straight-line edges (Line geometry or degree-1
    // NURBS or no curve at all). Anything richer is left untouched.
    let is_chord_candidate = match &edge_data.curve {
        Some(CurveOnEdge {
            geometry: CurveGeometry::Line { .. },
            ..
        }) => true,
        Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(n),
            ..
        }) => n.degree == 1 && n.control_points.len() == 2,
        None => true,
        _ => return false,
    };
    if !is_chord_candidate {
        return false;
    }

    let he = edge_data.half_edge;
    let Some(he_data) = brep.half_edges.get(he) else {
        return false;
    };
    let twin = he_data.twin;
    let v_start = he_data.origin;
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return false;
    };
    let v_end = twin_data.origin;
    let Some(start_v) = brep.vertices.get(v_start) else {
        return false;
    };
    let Some(end_v) = brep.vertices.get(v_end) else {
        return false;
    };
    let p1 = start_v.position;
    let p2 = end_v.position;

    let chord = Vec3::new(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    let chord_len = chord.length();
    if chord_len < 1e-9 {
        return false;
    }

    let axis = unit(cyl_axis);
    let rel1 = Vec3::new(
        p1.x - cyl_origin.x,
        p1.y - cyl_origin.y,
        p1.z - cyl_origin.z,
    );
    let rel2 = Vec3::new(
        p2.x - cyl_origin.x,
        p2.y - cyl_origin.y,
        p2.z - cyl_origin.z,
    );
    let axial1 = rel1.dot(&axis);
    let axial2 = rel2.dot(&axis);

    // Detect axial / transverse / helical edges by the axial-vs-
    // angular split between the endpoints.
    let radial1 = sub_vec(rel1, scale_vec(axis, axial1));
    let radial2 = sub_vec(rel2, scale_vec(axis, axial2));
    let r1_len = radial1.length();
    let r2_len = radial2.length();
    if r1_len < cyl_radius * 1e-3 || r2_len < cyl_radius * 1e-3 {
        // Endpoints sit on the axis line itself — degenerate.
        return false;
    }

    // Angle subtended between the two radial directions, in [0, π].
    let cos_subtend = (radial1.dot(&radial2) / (r1_len * r2_len)).clamp(-1.0, 1.0);
    let subtend = cos_subtend.acos();

    // Axial line (same circumferential angle): the chord lies on
    // the cylinder surface — no rewrite needed.
    if subtend < 1e-4 {
        return false;
    }

    // Helical edges (significant axial-and-angular span) require a
    // helix fit which CurveGeometry::Arc can't represent. Skip.
    let axial_span = (axial1 - axial2).abs();
    let angular_span_arc_len = cyl_radius * subtend;
    if axial_span > angular_span_arc_len.max(cyl_radius) * 0.5 {
        return false;
    }

    // Build the cylinder-cross-section midpoint.
    //
    // 1. Pick the axial position to use for the cross-section:
    //    average of the two endpoint axials. The two endpoints may
    //    sit at slightly different axial positions; the analytic
    //    arc lives at the average plane.
    // 2. Find the angular bisector between the two radial dirs,
    //    in the plane perpendicular to the cylinder axis. The
    //    bisector points OUTWARD from the cylinder axis along the
    //    centre of the (shorter) chord direction.
    // 3. Sample a point on the analytic cylinder cross-section at
    //    the bisector angle: `mid_on_cyl = centre + r * bisector`.
    let axial_mid = 0.5 * (axial1 + axial2);
    let centre = Point3::new(
        cyl_origin.x + axial_mid * axis.x,
        cyl_origin.y + axial_mid * axis.y,
        cyl_origin.z + axial_mid * axis.z,
    );
    let bisector = bisector_dir(radial1, radial2, axis, &p1, &p2);
    let mid_on_cyl = Point3::new(
        centre.x + bisector.x * cyl_radius,
        centre.y + bisector.y * cyl_radius,
        centre.z + bisector.z * cyl_radius,
    );

    // 3-point circle through (p1, mid_on_cyl, p2). If the three
    // points are collinear, abort.
    let Some((c, n_arc, r_arc)) = three_point_circle(p1, mid_on_cyl, p2) else {
        return false;
    };

    // Tangent at p1 (perpendicular to radius, in arc plane) and at
    // p2. Their intersection gives the rational quadratic Bezier's
    // middle CP.
    let radial_p1 = Vec3::new(p1.x - c.x, p1.y - c.y, p1.z - c.z);
    let radial_p2 = Vec3::new(p2.x - c.x, p2.y - c.y, p2.z - c.z);
    let tangent_p1 = unit(n_arc.cross(&radial_p1));
    let tangent_p2 = unit(n_arc.cross(&radial_p2));
    let Some(cp_mid) = ray_intersection_2d(p1, tangent_p1, p2, tangent_p2, n_arc) else {
        return false;
    };

    // Subtended half-angle gives the middle weight.
    let cos_half_angle = (radial_p1.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0)))
        .dot(&radial_p2.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0)))
        .clamp(-1.0, 1.0)
        .acos()
        * 0.5;
    let weight_mid = cos_half_angle.cos();
    if !weight_mid.is_finite() || weight_mid < 1e-4 {
        return false;
    }

    let nurbs = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![p1, cp_mid, p2],
        weights: vec![1.0, weight_mid, 1.0],
    };

    // Sanity-check: the NURBS's parametric endpoints must coincide
    // with p1 and p2 to within tolerance, otherwise downstream
    // ShellGap validation will trip.
    use popkern_core::geometry::curve::CurveEval;
    let test_start = nurbs.point_at(0.0);
    let test_end = nurbs.point_at(1.0);
    let endpoint_err = test_start.distance_to(&p1).max(test_end.distance_to(&p2));
    if endpoint_err > (chord_len * 0.001).max(1e-9) {
        return false;
    }
    // Sanity-check: the arc midpoint should be near the cylinder
    // surface. A wildly diverging midpoint indicates we built the
    // arc the wrong way around the circle — skip rather than
    // ship a deformed boundary.
    let test_mid = nurbs.point_at(0.5);
    let mid_rel = Vec3::new(
        test_mid.x - cyl_origin.x,
        test_mid.y - cyl_origin.y,
        test_mid.z - cyl_origin.z,
    );
    let mid_axial = mid_rel.dot(&axis);
    let mid_radial = sub_vec(mid_rel, scale_vec(axis, mid_axial)).length();
    if (mid_radial - cyl_radius).abs() > cyl_radius * 0.5 {
        let _ = r_arc;
        return false;
    }

    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start: 0.0,
            t_end: 1.0,
        });
        return true;
    }
    let _ = r_arc;
    false
}

/// Bisector direction in the cylinder cross-section plane,
/// pointing toward the chord midpoint side. Used to pick which of
/// the two possible 3-point arcs through `(p1, p2)` to build.
fn bisector_dir(radial1: Vec3, radial2: Vec3, axis: Vec3, p1: &Point3, p2: &Point3) -> Vec3 {
    let r1 = radial1.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let r2 = radial2.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    let mut sum = Vec3::new(r1.x + r2.x, r1.y + r2.y, r1.z + r2.z);
    let len = sum.length();
    if len < 1e-9 {
        // The two radials point opposite ways — chord is a diameter.
        // Pick a perpendicular bisector by rotating r1 90° about
        // the cylinder axis toward the chord midpoint.
        let chord_mid = Vec3::new(
            0.5 * (p1.x + p2.x),
            0.5 * (p1.y + p2.y),
            0.5 * (p1.z + p2.z),
        );
        let candidate = axis.cross(&r1);
        let dot_with_mid = candidate.dot(&chord_mid);
        sum = if dot_with_mid >= 0.0 {
            candidate
        } else {
            Vec3::new(-candidate.x, -candidate.y, -candidate.z)
        };
    }
    unit(sum)
}

/// Compute the unique circumcircle through 3 non-collinear points.
fn three_point_circle(p1: Point3, p2: Point3, p3: Point3) -> Option<(Point3, Vec3, f64)> {
    let a = Vec3::new(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    let b = Vec3::new(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
    let aa = a.dot(&a);
    let bb = b.dot(&b);
    let ab = a.dot(&b);
    let denom = 2.0 * (aa * bb - ab * ab);
    if denom.abs() < 1e-18 {
        return None;
    }
    let alpha = bb * (aa - ab) / denom;
    let beta = aa * (bb - ab) / denom;
    let centre_minus_p1 = Vec3::new(
        a.x * alpha + b.x * beta,
        a.y * alpha + b.y * beta,
        a.z * alpha + b.z * beta,
    );
    let centre = Point3::new(
        p1.x + centre_minus_p1.x,
        p1.y + centre_minus_p1.y,
        p1.z + centre_minus_p1.z,
    );
    let radius = centre_minus_p1.length();
    let normal_raw = a.cross(&b);
    let normal = unit(normal_raw);
    Some((centre, normal, radius))
}

/// Intersection of two rays in a common plane (defined by `normal`).
/// Both rays start at given origins and travel along given (unit)
/// directions. Returns the intersection point or `None` if the rays
/// are parallel.
fn ray_intersection_2d(o1: Point3, d1: Vec3, o2: Point3, d2: Vec3, normal: Vec3) -> Option<Point3> {
    // 2D solve in the plane: parameterise line 1 as o1 + t·d1, line
    // 2 as o2 + s·d2. Project everything to the plane via the
    // normal cross trick.
    let cross_d = d1.cross(&d2);
    let cross_n = cross_d.dot(&normal);
    if cross_n.abs() < 1e-12 {
        return None;
    }
    let delta = Vec3::new(o2.x - o1.x, o2.y - o1.y, o2.z - o1.z);
    let t_num = delta.cross(&d2).dot(&normal);
    let t = t_num / cross_n;
    Some(Point3::new(
        o1.x + d1.x * t,
        o1.y + d1.y * t,
        o1.z + d1.z * t,
    ))
}

fn unit(v: Vec3) -> Vec3 {
    let len = v.length();
    if len < 1e-12 {
        return Vec3::new(0.0, 0.0, 1.0);
    }
    Vec3::new(v.x / len, v.y / len, v.z / len)
}

fn sub_vec(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn scale_vec(v: Vec3, s: f64) -> Vec3 {
    Vec3::new(v.x * s, v.y * s, v.z * s)
}
