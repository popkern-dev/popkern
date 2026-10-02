//! The variational-fit cap's sector faces (PK-574): the carved region
//! split into the fit's sectors by spokes from the centre, each face on
//! its sector surface and each rim edge on its sector's rim iso-curve.
//! Split from `setback_fit_fan` for the file budget; it shares that
//! module's helpers.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::variational_cap::VariationalCap;
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, NurbsSurface, SurfaceOnFace};

use super::*;

/// Split the blend face into the fit's sectors; returns the sector
/// faces in sector order.
pub(super) fn build_sector_faces(
    brep: &mut BRep,
    blend_face: FaceHandle,
    cap: &VariationalCap,
    sides: &[popkern_core::geometry::variational_cap::VariationalCapSide],
) -> Result<Vec<FaceHandle>, String> {
    let n = cap.sectors.len();
    // 1. Spoke vertices on the sides, moved exactly onto the fit's m_s.
    let mut m_vertices = Vec::with_capacity(n);
    for &m in &cap.spoke_points {
        let v = split_blend_boundary_at(brep, blend_face, m)?;
        brep.vertices[v].position = m;
        m_vertices.push(v);
    }
    // 2. Spokes: o by mev from m_0, the rest closed onto o by mef.
    let o = cap.centre;
    let mut o_vertex: Option<VertexHandle> = None;
    for (s, &m_v) in m_vertices.iter().enumerate() {
        let spoke = spoke_curve(&cap.sectors[s]);
        match o_vertex {
            None => {
                // The half-edge leaving m_0 on the region's own loop: a
                // corner patch's seam vertex borders two cap regions.
                let he = face_loop_half_edges(brep, blend_face)?
                    .into_iter()
                    .find(|&h| brep.half_edges[h].origin == m_v)
                    .ok_or("spoke foot not on the region's loop")?;
                let mev = brep.mev(he, o).map_err(|e| format!("mev spoke 0: {e}"))?;
                brep.edges[mev.edge].curve = Some(spoke);
                o_vertex = Some(mev.vertex);
            }
            Some(o_v) => {
                let (he_m, he_o) = half_edge_pair_in_common_loop(brep, m_v, o_v)?;
                let mef = brep
                    .mef(he_m, he_o)
                    .map_err(|e| format!("mef spoke {s}: {e}"))?;
                brep.edges[mef.edge].curve = Some(spoke);
            }
        }
    }
    let o_vertex = o_vertex.ok_or("no spoke vertex")?;
    brep.vertices[o_vertex].position = o;
    // 3. Sector faces: the faces around o, matched to sectors by their
    //    corners (o, m_s, c_s, m_{s+1}).
    let around: Vec<FaceHandle> = {
        let mut faces: Vec<FaceHandle> = brep
            .half_edges
            .values()
            .filter(|h| h.origin == o_vertex)
            .map(|h| brep.loops[h.loop_].face)
            .collect();
        faces.sort();
        faces.dedup();
        faces
    };
    if around.len() != n {
        return Err(format!(
            "{} faces around the centre, expected {n}",
            around.len()
        ));
    }
    let mut assigned: Vec<Option<FaceHandle>> = vec![None; n];
    for &face in &around {
        let verts: Vec<Point3> = face_loop_half_edges(brep, face)?
            .iter()
            .map(|&he| brep.vertices[brep.half_edges[he].origin].position)
            .collect();
        let score = |s: usize| {
            let sec = &cap.sectors[s];
            [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                .iter()
                .map(|&(u, v)| {
                    let p = sec.point_at(u, v);
                    verts
                        .iter()
                        .map(|q| q.distance_to(&p))
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(0.0, f64::max)
        };
        let (best, d) = (0..n)
            .map(|s| (s, score(s)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .ok_or("no sector")?;
        if d > 1e-6 || assigned[best].is_some() {
            return Err(format!(
                "face {face:?} matches no free sector (best {best} at {d:.3e})"
            ));
        }
        assigned[best] = Some(face);
    }
    let faces: Vec<FaceHandle> = assigned
        .into_iter()
        .collect::<Option<_>>()
        .ok_or("unassigned sector")?;
    // 4. Surfaces, oriented by the neighbour's outward normal at m_s.
    for (s, &face) in faces.iter().enumerate() {
        let sec = &cap.sectors[s];
        let side = &sides[s];
        let p = sec.point_at(1.0, 0.5);
        let n_cap = sec.partial_u(1.0, 0.5).cross(&sec.partial_v(1.0, 0.5));
        let (u, v) =
            popkern_core::geometry::uv_projection::project_point_to_uv(&side.neighbour, &p)
                .ok_or("sector rim does not project onto its neighbour")?;
        let mut n_nb = side.neighbour.normal_at(u, v);
        if side.neighbour_reversed {
            n_nb = n_nb * -1.0;
        }
        brep.faces[face].surface = Some(SurfaceOnFace {
            geometry: cap.faces[s].clone(),
            orientation: n_cap.dot(&n_nb) > 0.0,
        });
    }
    // 5. Rim edges: each non-spoke cap edge takes its sector's rim
    //    iso-curve over its own parameter window; rim vertices snap onto
    //    it.
    for (s, &face) in faces.iter().enumerate() {
        let sec = &cap.sectors[s];
        let rims = [sec.iso_curve_at_u(1.0), sec.iso_curve_at_v(1.0)];
        for he in face_loop_half_edges(brep, face)? {
            let edge = brep.half_edges[he].edge;
            let (a, b) = (
                brep.half_edges[he].origin,
                brep.half_edges[brep.half_edges[he].twin].origin,
            );
            if a == o_vertex || b == o_vertex {
                continue;
            }
            let (pa, pb) = (brep.vertices[a].position, brep.vertices[b].position);
            // A rim piece lies on one iso-curve (its corner end sits on
            // both): take the rim that carries both ends.
            let (rim, ta, tb, d) = rim_carrying(&rims, pa, pb);
            if d > 1e-3 {
                return Err(format!(
                    "rim edge {d:.3e} off sector {s}'s rims (ends {pa:?}, {pb:?})"
                ));
            }
            brep.vertices[a].position = rims[rim].point_at(ta);
            brep.vertices[b].position = rims[rim].point_at(tb);
            brep.edges[edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(rims[rim].clone()),
                t_start: ta.min(tb),
                t_end: ta.max(tb),
            });
        }
    }
    Ok(faces)
}

/// Spoke `s` as its sector's `v = 0` iso-curve (o at u = 0, m_s at 1).
fn spoke_curve(sector: &NurbsSurface) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(sector.iso_curve_at_v(0.0)),
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// The rim iso-curve carrying both ends of a rim edge: the one whose
/// worse end distance is smallest, with both end parameters and that
/// distance.
fn rim_carrying(
    rims: &[popkern_core::geometry::NurbsCurve; 2],
    pa: Point3,
    pb: Point3,
) -> (usize, f64, f64, f64) {
    let fit = |k: usize| {
        let (ta, tb) = (param_on(&rims[k], pa), param_on(&rims[k], pb));
        let d = rims[k]
            .point_at(ta)
            .distance_to(&pa)
            .max(rims[k].point_at(tb).distance_to(&pb));
        (k, ta, tb, d)
    };
    let (first, second) = (fit(0), fit(1));
    if first.3 <= second.3 { first } else { second }
}
