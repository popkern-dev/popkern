//! Chamfer offset-row construction: interior-hint selection
//! (per-station planar Newell with measured sign calibration),
//! offset distances, and surface-projected offset rows.
//!
//! Split from `chamfer_curve.rs` to honour the structural ratchet.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::analysis::EdgeAnalysis;
use crate::types::ChamferKind;

/// Which side of ±(normal × tangent) is the face's material.
pub(super) enum InteriorHint {
    /// Constant direction — planar faces (whole-loop Newell, tuned).
    Constant(Vec3),
    /// Per-sample: toward this point (the face's boundary-sample
    /// centroid). A CLOSED rim's interior direction ROTATES around
    /// the ring — a cap's inward radial flips sign at the antipode —
    /// so no constant vector can sign-select every sample.
    Toward(Point3),
    /// Per-sample from a PLANAR face's material Newell: the inward
    /// direction at each station is newell × tangent — a constant
    /// `newell × tangent_at_midpoint` only sign-selects correctly
    /// while the tangent stays within 90° of the midpoint's, and an
    /// interrupted circular rim rotates its tangent ~330° (measured:
    /// the notched flange's rim chamfer offset half its stations
    /// OUTWARD, the conical-bevel fit read a 1.6 spread, and the
    /// swept fallback's terminal cut shipped a 0.29 curve lie).
    PlanarNewell(Vec3),
}

/// Compute offset distances on each face based on chamfer kind.
pub(super) fn compute_offset_distances(
    analysis: &EdgeAnalysis,
    kind: &ChamferKind,
) -> Result<(f64, f64), String> {
    match kind {
        ChamferKind::Equal(d) => Ok((*d, *d)),
        ChamferKind::TwoDistance {
            face_1: ref_face,
            dist_1,
            dist_2,
        } => {
            if *ref_face == analysis.face_1 {
                Ok((*dist_1, *dist_2))
            } else if *ref_face == analysis.face_2 {
                Ok((*dist_2, *dist_1))
            } else {
                Err("reference face does not match either adjacent face".into())
            }
        }
        ChamferKind::DistanceAngle {
            face: ref_face,
            distance,
            angle,
        } => {
            let other_dist = distance * angle.tan();
            if *ref_face == analysis.face_1 {
                Ok((*distance, other_dist))
            } else if *ref_face == analysis.face_2 {
                Ok((other_dist, *distance))
            } else {
                Err("reference face does not match either adjacent face".into())
            }
        }
    }
}

/// Compute offset points on a surface at a given distance from edge points.
///
/// For each edge point, computes the direction on the surface perpendicular
/// to the edge and offsets along that direction by the given distance.
/// Interior direction hint for `face` at the midpoint of `edge`:
/// outer-loop Newell normal × the face's own half-edge direction.
/// The Newell normal points outward by the loop-winding invariant,
/// so the cross product points into the face's material side.
pub(super) fn face_interior_hint(
    brep: &popkern_core::topology::BRep,
    edge: popkern_core::topology::EdgeHandle,
    face: popkern_core::topology::FaceHandle,
) -> Option<InteriorHint> {
    let face_data = brep.faces.get(face)?;
    // Newell normal over the outer loop.
    let start = brep.loops.get(face_data.outer_loop)?.half_edge;
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    let mut he = start;
    let mut he_on_edge = None;
    for _ in 0..brep.half_edges.len() + 1 {
        let d = brep.half_edges.get(he)?;
        let next = brep.half_edges.get(d.next)?;
        let p = brep.vertices.get(d.origin)?.position;
        let q = brep.vertices.get(next.origin)?.position;
        newell.x += (p.y - q.y) * (p.z + q.z);
        newell.y += (p.z - q.z) * (p.x + q.x);
        newell.z += (p.x - q.x) * (p.y + q.y);
        if d.edge == edge {
            he_on_edge = Some((p, q));
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    // The blended edge may be on an inner loop (hole rim chamfers).
    if he_on_edge.is_none() {
        for &lh in &face_data.inner_loops {
            let start = brep.loops.get(lh)?.half_edge;
            let mut he = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let d = brep.half_edges.get(he)?;
                let next = brep.half_edges.get(d.next)?;
                if d.edge == edge {
                    let p = brep.vertices.get(d.origin)?.position;
                    let q = brep.vertices.get(next.origin)?.position;
                    he_on_edge = Some((p, q));
                    break;
                }
                he = d.next;
                if he == start {
                    break;
                }
            }
            if he_on_edge.is_some() {
                break;
            }
        }
    }
    let (p, q) = he_on_edge?;
    // Midpoint and travel direction along this face's half-edge. On
    // a CLOSED edge (a rim: origin == dest) the vertex chord is
    // degenerate — sample the edge CURVE instead, flipped when this
    // face's half-edge is the reversed use.
    let (mid, he_dir) = if let Some(dir) = (q - p).normalized() {
        (p + (q - p) * 0.5, dir)
    } else {
        use popkern_core::geometry::curve::CurveEval;
        let edge_data = brep.edges.get(edge)?;
        let c = edge_data.curve.as_ref()?;
        let forward = {
            // The canonical half-edge is the forward use; find which
            // use belongs to this face.
            let canon = edge_data.half_edge;
            let canon_face = brep
                .half_edges
                .get(canon)
                .and_then(|h| brep.loops.get(h.loop_))
                .map(|l| l.face);
            canon_face == Some(face)
        };
        let t_mid = 0.5 * (c.t_start + c.t_end);
        let dt = (c.t_end - c.t_start) * 1e-4;
        let a = c.geometry.point_at(t_mid - dt);
        let b = c.geometry.point_at(t_mid + dt);
        let dir = (b - a).normalized()?;
        (
            c.geometry.point_at(t_mid),
            if forward { dir } else { dir * -1.0 },
        )
    };
    // The loop-winding invariant: the whole-loop Newell points INTO
    // the material on PLANAR faces (tuned across the chamfer suite)
    // — hint = newell × he_dir. On CURVED charts a whole-loop
    // Newell is meaningless (a periodic wall's two rims cancel or
    // flip it — measured: a cylinder rim chamfer walked the wall
    // trim UP off the solid, 1.5 above the cap), and normal-flag
    // conventions proved treacherous (both signs measured wrong on
    // different fixtures). Use MEASURED geometry instead: the
    // direction from the edge midpoint toward the face's
    // boundary-sample centroid, projected into the tangent plane.
    // The hint only selects the SIGN of ±(normal × tangent), so
    // centroid-side is decisive for any face whose bulk lies on its
    // material side — every chamfer support face qualifies. A
    // degenerate planar Newell (a cap bounded by one closed rim)
    // routes through the centroid path too.
    let surf = face_data.surface.as_ref()?;
    let planar_newell = matches!(surf.geometry, SurfaceGeometry::Plane { .. })
        .then(|| newell.normalized())
        .flatten();
    if let Some(n) = planar_newell {
        // Per-station inward = n_eff × curve_tangent. The SIGN of
        // n_eff is calibrated by MEASUREMENT at the midpoint station
        // — the loop-winding convention proved wrong on the notched
        // flange's interrupted rim (the constant-hint rule offset
        // half the stations outward, and the he-orientation rule
        // offset all of them outward): stepping from the edge
        // midpoint toward the face's boundary-sample centroid must
        // agree with the candidate inward direction.
        use popkern_core::geometry::curve::CurveEval;
        // Tangent AND position from the same curve-midpoint station —
        // the vertex-chord midpoint of a wide arc sits on the wrong
        // side of the face and inverts the calibration.
        let station = brep.edges.get(edge).and_then(|edge_data| {
            let c = edge_data.curve.as_ref()?;
            let t_mid = 0.5 * (c.t_start + c.t_end);
            let dt = ((c.t_end - c.t_start) * 1e-4).max(1e-9);
            let d =
                (c.geometry.point_at(t_mid + dt) - c.geometry.point_at(t_mid - dt)).normalized()?;
            Some((c.geometry.point_at(t_mid), d))
        });
        let samples = brep.face_boundary_sample_points(face, 8);
        if let (Some((mid_station, d)), false) = (station, samples.is_empty()) {
            let mut c = popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
            for p in &samples {
                c.x += p.x;
                c.y += p.y;
                c.z += p.z;
            }
            let k = samples.len() as f64;
            c.x /= k;
            c.y /= k;
            c.z /= k;
            let toward = c - mid_station;
            let candidate = n.cross(&d);
            let s_sign = if candidate.dot(&toward) >= 0.0 {
                1.0
            } else {
                -1.0
            };
            return Some(InteriorHint::PlanarNewell(n * s_sign));
        }
        return n.cross(&he_dir).normalized().map(InteriorHint::Constant);
    }
    let _ = (mid, he_dir, surf);
    let samples = brep.face_boundary_sample_points(face, 8);
    if samples.is_empty() {
        return None;
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for s in &samples {
        cx += s.x;
        cy += s.y;
        cz += s.z;
    }
    let n = samples.len() as f64;
    Some(InteriorHint::Toward(popkern_core::geometry::Point3::new(
        cx / n,
        cy / n,
        cz / n,
    )))
}

pub(super) fn compute_offset_points(
    edge_points: &[Point3],
    edge_tangents: &[Vec3],
    surface: &SurfaceGeometry,
    distance: f64,
    interior_hint: Option<&InteriorHint>,
    tolerance: &Tolerance,
) -> Result<Vec<Point3>, String> {
    if edge_points.len() < 2 {
        return Err("need at least 2 edge points".into());
    }

    let mut offset_points = Vec::with_capacity(edge_points.len());

    for (i, pt) in edge_points.iter().enumerate() {
        // Exact curve tangent; chord fallback only where the
        // derivative degenerates (see the sampling loop above).
        let exact = edge_tangents
            .get(i)
            .copied()
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        let tangent = if exact.length() > tolerance.numerical_floor() {
            exact
        } else if i == 0 {
            edge_points[1] - edge_points[0]
        } else if i == edge_points.len() - 1 {
            edge_points[i] - edge_points[i - 1]
        } else {
            edge_points[i + 1] - edge_points[i - 1]
        };
        let tangent_len = tangent.length();
        if tangent_len < tolerance.numerical_floor() {
            offset_points.push(*pt);
            continue;
        }
        let tangent_unit = Vec3::new(
            tangent.x / tangent_len,
            tangent.y / tangent_len,
            tangent.z / tangent_len,
        );

        // Find the surface normal at this point.
        let (u, v) = find_uv_for_point(surface, pt);
        let normal = surface.normal_at(u, v);
        let normal_len = normal.length();
        if normal_len < tolerance.numerical_floor() {
            offset_points.push(*pt);
            continue;
        }

        // Offset direction: perpendicular to the edge tangent, lying on the surface.
        // offset_dir = normalize(normal x tangent)
        let offset_raw = normal.cross(&tangent_unit);
        let offset_len = offset_raw.length();
        if offset_len < tolerance.numerical_floor() {
            offset_points.push(*pt);
            continue;
        }
        let mut offset_dir = Vec3::new(
            offset_raw.x / offset_len,
            offset_raw.y / offset_len,
            offset_raw.z / offset_len,
        );
        // Walk into the face's material side. The normal component
        // of a Toward vector is perpendicular to offset_dir, so the
        // raw dot sign-selects correctly without projection.
        let hint_vec = match interior_hint {
            Some(InteriorHint::Constant(v)) => Some(*v),
            Some(InteriorHint::Toward(c)) => Some(*c - *pt),
            Some(InteriorHint::PlanarNewell(n)) => Some(n.cross(&tangent_unit)),
            None => None,
        };
        if let Some(hint) = hint_vec
            && offset_dir.dot(&hint) < 0.0
        {
            offset_dir = offset_dir * -1.0;
        }

        let offset_pt = *pt + offset_dir * distance;

        // Project offset point back onto the surface for accuracy.
        let projected = project_onto_surface(surface, &offset_pt);
        offset_points.push(projected);
    }

    Ok(offset_points)
}

/// Project a 3D point onto a surface (closest point).
pub(super) fn project_onto_surface(surface: &SurfaceGeometry, point: &Point3) -> Point3 {
    let (u, v) = find_uv_for_point(surface, point);
    surface.point_at(u, v)
}

/// Find (u, v) parameters on a surface closest to a 3D point.
pub(super) fn find_uv_for_point(surface: &SurfaceGeometry, point: &Point3) -> (f64, f64) {
    if let Some(uv) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, point) {
        return uv;
    }

    // Fallback: grid search.
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    let n = 20;
    let mut best_u = (u_lo + u_hi) * 0.5;
    let mut best_v = (v_lo + v_hi) * 0.5;
    let mut best_dist = f64::INFINITY;

    for iu in 0..n {
        let u = u_lo + (u_hi - u_lo) * (iu as f64) / ((n - 1).max(1) as f64);
        for iv in 0..n {
            let v = v_lo + (v_hi - v_lo) * (iv as f64) / ((n - 1).max(1) as f64);
            let sp = surface.point_at(u, v);
            let d = sp.distance_to(point);
            if d < best_dist {
                best_dist = d;
                best_u = u;
                best_v = v;
            }
        }
    }

    (best_u, best_v)
}
