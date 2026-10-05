//! Spindle-torus pinch caps (PK-327 design §5 Stage 4, §6 "spindle
//! torus after offset").
//!
//! A torus offset whose tube swallows the axis (`minor ≥ major`) is an
//! immersed spindle torus: the realized Minkowski boundary is the
//! APPLE — the surface band `v ∈ (−θ₀, θ₀)` with `cos θ₀ = −R/r'` —
//! and the lemon (the rest of the chart) lies strictly inside the
//! offset solid. The trim boundary is DEGENERATE: the circles
//! `v = ±θ₀` have zero radius (they are the two pinch points on the
//! axis at `±√(r'² − R²)`), so no transversal intersection curve
//! exists for the generic scan — this analytic arm serves the class
//! with exact geometry instead: delete every face touching the dead
//! window, then fan-cap the two exposed rings onto exact pinch-pole
//! vertices with exact meridian ellipse arcs (`u_basis = radial·r'`,
//! `v_basis = axis·r'` — the tube angle is the parameter).

use popkern_math::PortableMath;
use std::collections::{BTreeMap, BTreeSet};

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::body_distance::FaceTrimEntry;
use crate::error::{ShellError, TrimRegionReport};

/// One spindled surface group: faces of the result sharing a spindle
/// torus.
#[derive(Debug)]
pub(crate) struct SpindleGroup {
    /// Group faces.
    pub faces: Vec<FaceHandle>,
    /// Torus center.
    pub center: Point3,
    /// Unit axis.
    pub axis: Vec3,
    /// Major radius R.
    pub major: f64,
    /// Offset minor radius r' (≥ R).
    pub minor: f64,
}

/// Detect spindled surface groups (faces whose surface is a torus with
/// `minor ≥ major`), grouped by identical torus geometry.
pub(crate) fn detect_spindle_groups(body: &BRep, tolerance: &Tolerance) -> Vec<SpindleGroup> {
    let tol = tolerance.linear();
    let mut groups: Vec<SpindleGroup> = Vec::new();
    let mut face_handles: Vec<FaceHandle> = body.faces.keys().collect();
    face_handles.sort();
    for fh in face_handles {
        let Some(sof) = body.faces[fh].surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } = &sof.geometry
        else {
            continue;
        };
        if *minor_radius < *major_radius {
            continue;
        }
        let axis_n = match axis.normalized() {
            Some(a) => a,
            None => continue,
        };
        let existing = groups.iter_mut().find(|g| {
            (g.center - *center).length() <= tol
                && (g.axis.dot(&axis_n).abs() - 1.0).abs() <= 1e-9
                && (g.major - *major_radius).abs() <= tol
                && (g.minor - *minor_radius).abs() <= tol
        });
        match existing {
            Some(g) => g.faces.push(fh),
            None => groups.push(SpindleGroup {
                faces: vec![fh],
                center: *center,
                axis: axis_n,
                major: *major_radius,
                minor: *minor_radius,
            }),
        }
    }
    groups
}

/// Region reports for detected spindle groups (the Decline-policy
/// payload): witness = the pinch point, face pair = the group's first
/// face against itself (same-surface class).
pub(crate) fn spindle_region_reports(groups: &[SpindleGroup]) -> Vec<TrimRegionReport> {
    groups
        .iter()
        .map(|g| {
            let w = (g.minor * g.minor - g.major * g.major).sqrt();
            TrimRegionReport {
                face_a: g.faces[0],
                face_b: g.faces[0],
                witness: g.center + g.axis * w,
                curve_length: 0.0,
            }
        })
        .collect()
}

/// Serve one spindle group: delete the dead-window faces and fan-cap
/// the exposed rings onto exact pinch poles.
///
/// # Errors
///
/// [`ShellError::Topology`] when the group's faces do not isolate into
/// two cap rings (partial spindle patches — the fixture-less residual
/// class, recorded in the decline-class registry).
pub(crate) fn serve_spindle_group(
    body: &mut BRep,
    group: &SpindleGroup,
    tolerance: &Tolerance,
) -> Result<usize, ShellError> {
    let (rr, r_minor) = (group.major, group.minor);
    let theta0 = (-rr / r_minor).pacos();
    let dead_lo = theta0;
    let dead_hi = 2.0 * std::f64::consts::PI - theta0;
    let w = (r_minor * r_minor - rr * rr).sqrt();

    // Orientation flag to inherit onto the cap sectors.
    let orientation = group
        .faces
        .first()
        .and_then(|fh| body.faces[*fh].surface.as_ref())
        .map(|s| s.orientation)
        .unwrap_or(true);
    let torus_geometry = body.faces[group.faces[0]]
        .surface
        .as_ref()
        .map(|s| s.geometry.clone())
        .ok_or_else(|| ShellError::Topology("spindle group lost its surface".into()))?;

    // Faces whose v-window touches the open dead zone die.
    let margin = 1e-9;
    let mut dead: BTreeSet<FaceHandle> = BTreeSet::new();
    for &fh in &group.faces {
        let Some(entry) = FaceTrimEntry::build(body, fh) else {
            continue;
        };
        let Some((_, v0, _, v1)) = entry.uv_bbox() else {
            continue;
        };
        if v1 > dead_lo + margin && v0 < dead_hi - margin {
            dead.insert(fh);
        }
    }
    if dead.is_empty() {
        return Ok(0);
    }
    let removed_count = dead.len();

    let pre_boundary: BTreeSet<EdgeHandle> = body.boundary_edges.clone();
    body.remove_face_region_to_boundary(&dead)
        .map_err(|e| ShellError::Topology(format!("spindle region removal failed: {e}")))?;

    // The new boundary must fall into exactly two rings, one per pole.
    let new_boundary: Vec<EdgeHandle> = body
        .boundary_edges
        .iter()
        .copied()
        .filter(|eh| !pre_boundary.contains(eh))
        .collect();
    let rings = connected_rings(body, &new_boundary);
    if rings.len() != 2 {
        return Err(ShellError::Topology(format!(
            "spindle pinch expects two cap rings, found {} (partial spindle patch — \
             decline-class registry row)",
            rings.len()
        )));
    }

    let (perp1, perp2) = orthonormal_basis(group.axis);
    for ring in rings {
        // Ring latitude from a sample point.
        let sample = ring_sample_point(body, &ring)?;
        let rel = sample - group.center;
        let z_loc = rel.dot(&group.axis);
        let radial = rel - group.axis * z_loc;
        let rho = radial.length();
        let v_ring = normalize_angle(z_loc.patan2(rho - rr));
        // Pole side: caps on the upper half close at v = θ₀, lower at
        // v = 2π − θ₀.
        let (v_pole, pole) = if v_ring < std::f64::consts::PI {
            (theta0, group.center + group.axis * w)
        } else {
            (dead_hi, group.center - group.axis * w)
        };

        let ring_vec: Vec<EdgeHandle> = ring.to_vec();
        let fan = body
            .cap_boundary_ring_fan(&ring_vec, pole)
            .map_err(|e| ShellError::Topology(format!("spindle pinch fan failed: {e}")))?;

        // Author sector surfaces, spoke curves, and pcurves.
        author_fan_geometry(
            body,
            &fan,
            &torus_geometry,
            orientation,
            group,
            (perp1, perp2),
            v_ring,
            v_pole,
            tolerance,
        )?;
    }
    Ok(removed_count)
}

/// Group boundary edges into connected rings by shared vertices.
fn connected_rings(body: &BRep, edges: &[EdgeHandle]) -> Vec<Vec<EdgeHandle>> {
    let mut vertex_edges: BTreeMap<VertexHandle, Vec<EdgeHandle>> = BTreeMap::new();
    let ends = |eh: EdgeHandle| -> Option<(VertexHandle, VertexHandle)> {
        let he = body.half_edges.get(body.edges[eh].half_edge)?;
        let next = body.half_edges.get(he.next)?;
        Some((he.origin, next.origin))
    };
    for &eh in edges {
        if let Some((a, b)) = ends(eh) {
            vertex_edges.entry(a).or_default().push(eh);
            vertex_edges.entry(b).or_default().push(eh);
        }
    }
    let mut seen: BTreeSet<EdgeHandle> = BTreeSet::new();
    let mut rings = Vec::new();
    for &start in edges {
        if seen.contains(&start) {
            continue;
        }
        let mut ring = Vec::new();
        let mut stack = vec![start];
        while let Some(eh) = stack.pop() {
            if !seen.insert(eh) {
                continue;
            }
            ring.push(eh);
            if let Some((a, b)) = ends(eh) {
                for v in [a, b] {
                    for &other in vertex_edges.get(&v).into_iter().flatten() {
                        if !seen.contains(&other) {
                            stack.push(other);
                        }
                    }
                }
            }
        }
        rings.push(ring);
    }
    rings
}

/// A representative 3D point on the ring (an edge-curve midpoint, or a
/// vertex position).
fn ring_sample_point(body: &BRep, ring: &[EdgeHandle]) -> Result<Point3, ShellError> {
    for &eh in ring {
        if let Some(c) = body.edges[eh].curve.as_ref() {
            return Ok(c.geometry.point_at((c.t_start + c.t_end) * 0.5));
        }
        if let Some(he) = body.half_edges.get(body.edges[eh].half_edge) {
            return Ok(body.vertices[he.origin].position);
        }
    }
    Err(ShellError::Topology(
        "spindle pinch: ring has no sampleable edge".into(),
    ))
}

fn normalize_angle(a: f64) -> f64 {
    let tau = std::f64::consts::TAU;
    let mut v = a % tau;
    if v < 0.0 {
        v += tau;
    }
    v
}

/// Attach torus surfaces, exact meridian ellipse arcs, and linear
/// pcurves to a fresh pinch fan.
#[allow(clippy::too_many_arguments)]
fn author_fan_geometry(
    body: &mut BRep,
    fan: &popkern_core::euler::BoundaryRingFanResult,
    torus_geometry: &SurfaceGeometry,
    orientation: bool,
    group: &SpindleGroup,
    (perp1, perp2): (Vec3, Vec3),
    v_ring: f64,
    v_pole: f64,
    tolerance: &Tolerance,
) -> Result<(), ShellError> {
    let tol = tolerance.linear();
    // Azimuth of a ring vertex.
    let azimuth = |p: Point3| -> f64 {
        let rel = p - group.center;
        normalize_angle(rel.dot(&perp2).patan2(rel.dot(&perp1)))
    };

    // Spoke curves: exact meridian ellipse arcs, ranged ring → pole in
    // the tube-angle parameter.
    for &spoke in &fan.spoke_edges {
        let he = body.edges[spoke].half_edge;
        let he_data = &body.half_edges[he];
        let origin = he_data.origin;
        let next_origin = body.half_edges[he_data.next].origin;
        let (ring_vertex, forward) = if origin == fan.pole_vertex {
            (next_origin, false)
        } else {
            (origin, true)
        };
        let p_ring = body.vertices[ring_vertex].position;
        let phi = azimuth(p_ring);
        let radial_dir = perp1 * phi.pcos() + perp2 * phi.psin();
        let geometry = CurveGeometry::Ellipse {
            center: group.center + radial_dir * group.major,
            u_basis: radial_dir * group.minor,
            v_basis: group.axis * group.minor,
            start_angle: 0.0,
            end_angle: std::f64::consts::TAU,
        };
        // Verify the frame reproduces both endpoints before committing.
        let at = |t: f64| geometry.point_at(t);
        let pole_pos = body.vertices[fan.pole_vertex].position;
        if at(v_ring).distance_to(&p_ring) > tol * 100.0
            || at(v_pole).distance_to(&pole_pos) > tol * 100.0
        {
            return Err(ShellError::Topology(format!(
                "spindle spoke frame mismatch at φ={phi:.4} (ring residual {:.3e}, \
                 pole residual {:.3e})",
                at(v_ring).distance_to(&p_ring),
                at(v_pole).distance_to(&pole_pos)
            )));
        }
        let (t_start, t_end) = if forward {
            (v_ring, v_pole)
        } else {
            (v_pole, v_ring)
        };
        body.edges[spoke].curve = Some(CurveOnEdge {
            geometry,
            t_start,
            t_end,
        });
    }

    // Sector faces: torus surface + linear pcurves (rim v-line, spoke
    // u-lines). The pole side of the chart collapses (ρ = 0 at
    // v = v_pole) exactly like a sphere pole; the trim polygon is the
    // three pcurves.
    for (k, &face) in fan.sector_faces.iter().enumerate() {
        body.faces[face].surface = Some(SurfaceOnFace {
            geometry: torus_geometry.clone(),
            orientation,
        });
        // Walk the sector loop and author each half-edge's pcurve from
        // its endpoint geometry.
        let loop_h = body.faces[face].outer_loop;
        let start_he = body.loops[loop_h].half_edge;
        let mut cursor = start_he;
        let bound = body.half_edges.len().max(1);
        let mut phi_prev: Option<f64> = None;
        for _ in 0..=bound {
            let he = body.half_edges[cursor].clone();
            let origin_pos = body.vertices[he.origin].position;
            let dest = body.half_edges[he.next].origin;
            let dest_pos = body.vertices[dest].position;
            let uv_of = |p: Point3, is_pole: bool, phi_hint: Option<f64>| -> (f64, f64) {
                if is_pole {
                    // The pole keeps the neighbour's azimuth so the
                    // pcurve polygon stays chart-coherent.
                    (phi_hint.unwrap_or(0.0), v_pole)
                } else {
                    let mut phi = azimuth(p);
                    if let Some(prev) = phi_hint {
                        // Unwrap toward the previous azimuth.
                        while phi - prev > std::f64::consts::PI {
                            phi -= std::f64::consts::TAU;
                        }
                        while prev - phi > std::f64::consts::PI {
                            phi += std::f64::consts::TAU;
                        }
                    }
                    (phi, v_ring)
                }
            };
            let o_is_pole = he.origin == fan.pole_vertex;
            let d_is_pole = dest == fan.pole_vertex;
            let (u0, v0) = uv_of(origin_pos, o_is_pole, phi_prev);
            let (u1, v1) = uv_of(dest_pos, d_is_pole, Some(u0));
            body.half_edges[cursor].pcurve = Some(Pcurve::linear((u0, v0), (u1, v1), 0.0, 1.0));
            phi_prev = Some(u1);
            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
        let _ = k;
    }
    Ok(())
}
