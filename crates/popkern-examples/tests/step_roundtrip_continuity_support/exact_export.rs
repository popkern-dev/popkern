//! Gates for the exact-Gregory export opt-in round trip (PK-556):
//! representation identity, sampled seam parity, and the OCC oracle.

use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle};

/// The OCC oracle verdict (`occ_check.py`: reader checks, BRepCheck,
/// no open free bounds) on `step`; `None` when the venv is absent.
pub(crate) fn occ_oracle(step: &std::path::Path) -> Option<bool> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let check = std::path::Path::new("../../tools/occ-oracle/occ_check.py");
    if !python.exists() || !check.exists() {
        return None;
    }
    let output = std::process::Command::new(&python)
        .arg(check)
        .arg(step)
        .output()
        .ok()?;
    Some(String::from_utf8_lossy(&output.stdout).contains("\"passed\": true"))
}

/// The imported face whose NURBS surface carries `exact`'s control
/// net (points and weights within `tol`, relative to the weight scale).
pub(crate) fn carries_net(imported: &[BRep], exact: &NurbsSurface, tol: f64) -> bool {
    let w_scale = exact.weights.iter().copied().fold(0.0_f64, f64::max);
    imported.iter().flat_map(|b| b.faces.values()).any(|face| {
        let Some(SurfaceGeometry::Nurbs(n)) = face.surface.as_ref().map(|s| &s.geometry) else {
            return false;
        };
        n.degree_u == exact.degree_u
            && n.degree_v == exact.degree_v
            && n.control_points.len() == exact.control_points.len()
            && n.weights
                .iter()
                .zip(&exact.weights)
                .all(|(a, b)| (a - b).abs() <= tol * w_scale)
            && n.control_points
                .iter()
                .zip(&exact.control_points)
                .all(|(a, b)| a.distance_to(b) <= tol)
    })
}

/// Sampled seam continuity (`measure_cap_seam_continuity`, 200
/// samples) of `edge`'s two faces: `(g1 normal sine, g2 residual)`.
pub(crate) fn sampled_seam(brep: &BRep, edge: EdgeHandle) -> Option<(f64, f64)> {
    let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    let faces = [
        brep.loops.get(he.loop_)?.face,
        brep.loops.get(twin.loop_)?.face,
    ];
    let m = popkern_blend::cap_seam_continuity::measure_cap_seam_continuity(brep, &faces, 200);
    m.measured()
        .then_some((m.worst_g1_normal, m.worst_g2_residual))
}
