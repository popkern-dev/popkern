//! Seam-honest pcurve authoring for the primitive torus (split from
//! `torus.rs` for the file-size cap; PK-304).

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::topology::{BRep, HalfEdgeHandle};
use popkern_math::PortableMath;

/// Attach pcurves to all half-edges on the torus.
///
/// For each half-edge, determines the face it belongs to, retrieves the
/// torus surface geometry, and projects the 3D start/end positions to UV
/// via [`project_point_to_uv`]. Builds a degree-1 NURBS line pcurve in
/// `Point3(u, v, 0)` convention.
pub(super) fn attach_torus_pcurves(brep: &mut BRep) -> Result<(), CoreError> {
    let he_handles: Vec<HalfEdgeHandle> = brep.half_edges.keys().collect();
    for heh in he_handles {
        let he = &brep.half_edges[heh];
        let loop_handle = he.loop_;
        let face_handle = brep.loops[loop_handle].face;

        let surface = match &brep.faces[face_handle].surface {
            Some(sof) => sof.geometry.clone(),
            None => continue,
        };

        // Get the 3D start position (half-edge origin vertex).
        let start_pos = brep.vertices[brep.half_edges[heh].origin].position;
        // Get the 3D end position (next half-edge's origin vertex).
        let next_heh = brep.half_edges[heh].next;
        let end_pos = brep.vertices[brep.half_edges[next_heh].origin].position;

        // Determine half-edge direction relative to the edge curve.
        // If the edge has curve geometry, check whether the half-edge
        // runs forward or backward along it, and use the curve's
        // start/end points for accurate UV projection at seams.
        let (p_start, p_end) = if let Some(ref curve) = brep.edges[brep.half_edges[heh].edge].curve
        {
            let curve_start = curve.geometry.point_at(curve.t_start);
            let curve_end = curve.geometry.point_at(curve.t_end);
            // Compare against the curve's OWN endpoints: the previous
            // comparison against the next half-edge's origin was a
            // perfect tie for every TWIN (its far end IS the canonical
            // origin), so floating-point noise decided the direction
            // and a handful of pcurves were authored reversed (PK-304
            // — three faces' rings failed exact chaining and lost
            // cell resolution).
            let forward = start_pos.distance_to(&curve_start) < start_pos.distance_to(&curve_end);
            if forward {
                (
                    curve.geometry.point_at(curve.t_start),
                    curve.geometry.point_at(curve.t_end),
                )
            } else {
                (
                    curve.geometry.point_at(curve.t_end),
                    curve.geometry.point_at(curve.t_start),
                )
            }
        } else {
            (start_pos, end_pos)
        };

        let (u0, v0) = project_point_to_uv(&surface, &p_start).ok_or_else(|| {
            CoreError::PreconditionFailed("UV projection failed for pcurve start".into())
        })?;
        let (u1, v1) = project_point_to_uv(&surface, &p_end).ok_or_else(|| {
            CoreError::PreconditionFailed("UV projection failed for pcurve end".into())
        })?;

        // Seam-honest branch selection (PK-304): `project_point_to_uv`
        // normalizes each angle into [0, 2π) INDEPENDENTLY, so an edge
        // crossing a seam got a literal straight pcurve sweeping the
        // wrong way across the whole chart (measured: the seam-band
        // quadrant faces' boundary rings spanned 3/4 of the WRONG side
        // and 44% of the torus chart was covered by no face — every
        // seam-hugging boolean intersection arc was dropped). Anchor
        // both endpoints to the branch nearest the edge's 3D MIDPOINT
        // projection; the pcurve may then leave [0, 2π), which is the
        // per-half-edge chart convention (PK-286) consumers fold.
        let p_mid = Point3::new(
            0.5 * (p_start.x + p_end.x),
            0.5 * (p_start.y + p_end.y),
            0.5 * (p_start.z + p_end.z),
        );
        let mid_uv = brep.edges[brep.half_edges[heh].edge]
            .curve
            .as_ref()
            .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)))
            .and_then(|p| project_point_to_uv(&surface, &p))
            .or_else(|| project_point_to_uv(&surface, &p_mid));
        let (u0, v0, u1, v1) = if let Some((um, vm)) = mid_uv {
            let (pu, pv) = surface.periodicity();
            let near = |x: f64, anchor: f64, period: Option<f64>| match period {
                Some(p) if p.is_finite() && p > 0.0 => x - ((x - anchor) / p).round() * p,
                _ => x,
            };
            (
                near(u0, um, pu),
                near(v0, vm, pv),
                near(u1, um, pu),
                near(v1, vm, pv),
            )
        } else {
            (u0, v0, u1, v1)
        };

        let pcurve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
            weights: vec![1.0, 1.0],
        };

        brep.half_edges[heh].pcurve = Some(pcurve.into());
    }
    Ok(())
}

/// Author exact 3D edge curves for the primitive torus (PK-304): every
/// edge is a circular arc — rim edges (constant tube angle) lie on
/// MAJOR circles about the axis, side edges (constant major angle) on
/// MINOR circles around the tube. The edges carried NO curves before
/// this (measured fallout: the split stage's T-junction pre-split had
/// no windows to split, the open-chain boundary test fell back to
/// vertex chords a sagitta off the rims, and edge-adjacent repairs had
/// nothing to reconcile against). Arc parameter t IS the angle in
/// `circle_basis(normal)`, and the CANONICAL half-edge's origin sits
/// at `t_start` (the kernel direction convention, PK-303).
pub(super) fn attach_torus_edge_curves(brep: &mut BRep, major_radius: f64) {
    use crate::geometry::CurveOnEdge;
    use crate::geometry::Vec3;
    use crate::geometry::analytic_curves::circle_basis;
    use crate::geometry::curve::CurveGeometry;

    let eps = 1.0e-9;
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    for eh in edge_handles {
        let he = brep.edges[eh].half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let a = brep.vertices[he_data.origin].position;
        let b = brep.vertices[brep.half_edges[twin].origin].position;
        let ra = (a.x * a.x + a.y * a.y).sqrt();
        let rb = (b.x * b.x + b.y * b.y).sqrt();

        let (center, normal, radius) = if (a.z - b.z).abs() < eps && (ra - rb).abs() < eps {
            // Rim edge: major circle at this height.
            (Point3::new(0.0, 0.0, a.z), Vec3::new(0.0, 0.0, 1.0), ra)
        } else {
            // Side edge: minor circle in this major-angle plane.
            // Snap near-zero components before atan2: a vertex at the
            // u=π plane carries y = ±1e-16 noise and the angle (hence
            // the circle frame) flipped sign between edges of one
            // face, breaking the pcurve ring's exact-multiple
            // chaining (faces at that sector lost cell resolution).
            let snap = |x: f64| if x.abs() < 1.0e-12 { 0.0 } else { x };
            let u = snap(a.y).patan2(snap(a.x));
            let center = Point3::new(major_radius * u.pcos(), major_radius * u.psin(), 0.0);
            let normal = Vec3::new(-u.psin(), u.pcos(), 0.0);
            let radius = {
                let d = a - center;
                d.length()
            };
            (center, normal, radius)
        };
        let (bu, bv) = circle_basis(normal);
        let angle_of = |p: &Point3| -> f64 {
            let d = *p - center;
            d.dot(&bv).patan2(d.dot(&bu))
        };
        let theta_a = angle_of(&a);
        let mut theta_b = angle_of(&b);
        // Shortest-way branch: primitive edges span a quarter turn.
        theta_b -= ((theta_b - theta_a) / std::f64::consts::TAU).round() * std::f64::consts::TAU;
        if (theta_b - theta_a).abs() < eps {
            continue; // degenerate (shouldn't happen on the primitive)
        }
        brep.edges[eh].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle: theta_a.min(theta_b),
                end_angle: theta_a.max(theta_b),
            },
            t_start: theta_a,
            t_end: theta_b,
        });
    }
}
