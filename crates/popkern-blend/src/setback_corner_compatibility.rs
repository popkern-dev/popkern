//! Setback-corner compatibility (PK-575): whether the two surfaces a
//! vertex-blend cap meets at a boundary corner agree well enough for
//! one cap to be G1 or G2 to both there.
//!
//! A cap boundary turns from neighbour face A to neighbour face B at a
//! corner vertex `c`. A single C² cap has ONE normal and ONE second
//! fundamental form at `c`, so it can be G1 to both neighbours only if
//! their normals agree at `c`, and G2 to both only if their second
//! fundamental forms agree too. When they don't, the residual at `c` is
//! the INPUT's, and no cap construction can remove it. This record makes
//! that a named, structured fact on the result rather than something
//! absorbed into a rim tolerance.
//!
//! Measured on the apex family (flagship, pentagon5/hexagon6 blunt,
//! hexagon6_mixed, pentagon5_variable): every setback corner agrees to
//! |ΔII| <= 2e-9, since both strips are flat where the rim meets them.
//! The large rim residuals PK-575 retired were construction-limited.

use popkern_core::geometry::certify::{Iv, SurfaceJetSource};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend_continuity_evidence::ContinuityGates;

/// One cap-boundary corner's input compatibility.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SetbackCornerCompatibility {
    /// The corner vertex in the result body.
    pub vertex: VertexHandle,
    /// Its position.
    pub point: Point3,
    /// The two neighbour faces the cap boundary turns between.
    pub neighbours: (FaceHandle, FaceHandle),
    /// Sine of the angle between the neighbours' normals at the corner.
    pub normal_sine: f64,
    /// Largest normal-curvature difference between the neighbours over
    /// every tangent direction at the corner (spectral norm of ΔII).
    pub form_jump: f64,
    /// Highest continuity order a single cap can meet at this corner
    /// against the gates: 2 (compatible to second order), 1 (normals
    /// agree, curvatures do not), 0 (the neighbours crease here).
    pub compatible_order: u8,
}

/// Exact normal and second-form operator of `surface` at the point
/// nearest `p`: `(n, [II(e1,e1), II(e1,e2), II(e2,e2)], (e1, e2))`, in
/// the frame `frame` when given, else `e1 = Ŝ_u`, `e2 = n × e1`.
type SecondForm = (Vec3, [f64; 3], (Vec3, Vec3));

fn second_form_at(
    surface: &SurfaceGeometry,
    p: &Point3,
    frame: Option<(Vec3, Vec3)>,
) -> Option<SecondForm> {
    let (u, v) = project_point_to_uv(surface, p)?;
    let mut source = SurfaceJetSource::new(surface).ok()?;
    let (iu, iv) = (Iv::point(u), Iv::point(v));
    let span = source.spans_touching(iu, iv).ok()?.first().copied()?;
    let jet = source.jet_in(iu, iv, span).ok()?;
    let (su, sv) = (jet.partial(1, 0).mid(), jet.partial(0, 1).mid());
    let (suu, suv, svv) = (
        jet.partial(2, 0).mid(),
        jet.partial(1, 1).mid(),
        jet.partial(0, 2).mid(),
    );
    let n = su.cross(&sv).normalized()?;
    let (e1, e2) = match frame {
        Some(f) => f,
        None => {
            let e1 = su.normalized()?;
            (e1, n.cross(&e1))
        }
    };
    let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = e * g - f_ * f_;
    if det <= 1e-24 * e * g {
        return None;
    }
    let chart = |d: Vec3| {
        let (r1, r2) = (d.dot(&su), d.dot(&sv));
        ((r1 * g - r2 * f_) / det, (r2 * e - r1 * f_) / det)
    };
    let (l, m, nn) = (suu.dot(&n), suv.dot(&n), svv.dot(&n));
    let ii = |a: (f64, f64), b: (f64, f64)| {
        a.0 * l * b.0 + (a.0 * m * b.1 + a.1 * m * b.0) + a.1 * nn * b.1
    };
    let (w1, w2) = (chart(e1), chart(e2));
    Some((n, [ii(w1, w1), ii(w1, w2), ii(w2, w2)], (e1, e2)))
}

/// Spectral norm of the symmetric 2×2 form `[[a, b], [b, c]]`.
fn spectral_norm(a: f64, b: f64, c: f64) -> f64 {
    0.5 * (a + c).abs() + (0.25 * (a - c) * (a - c) + b * b).sqrt()
}

/// Every corner of the cap faces' outer boundary where the boundary
/// turns between two different non-cap neighbours, with the
/// neighbours' compatibility there against `gates` (G1 on the normal
/// sine, G2 on the second-form jump).
pub fn setback_corner_compatibility(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    gates: &ContinuityGates,
) -> Vec<SetbackCornerCompatibility> {
    let caps: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    // Cap-boundary seams: (non-cap neighbour, endpoint vertices).
    let mut rims: Vec<(FaceHandle, [VertexHandle; 2])> = Vec::new();
    for edge in brep.edges.values() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
            continue;
        };
        let (fa, fb) = (la.face, lb.face);
        let neighbour = match (caps.contains(&fa), caps.contains(&fb)) {
            (true, false) => fb,
            (false, true) => fa,
            _ => continue,
        };
        rims.push((neighbour, [he.origin, tw.origin]));
    }
    let mut out = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for (i, a) in rims.iter().enumerate() {
        for b in &rims[i + 1..] {
            if a.0 == b.0 {
                continue;
            }
            let Some(&vertex) = a.1.iter().find(|v| b.1.contains(v)) else {
                continue;
            };
            if !seen.insert((vertex, a.0.min(b.0), a.0.max(b.0))) {
                continue;
            }
            let Some(point) = brep.vertices.get(vertex).map(|v| v.position) else {
                continue;
            };
            let surface = |f: FaceHandle| {
                brep.faces
                    .get(f)
                    .and_then(|x| x.surface.as_ref())
                    .map(|s| &s.geometry)
            };
            let (Some(sa), Some(sb)) = (surface(a.0), surface(b.0)) else {
                continue;
            };
            let Some((na, ii_a, frame)) = second_form_at(sa, &point, None) else {
                continue;
            };
            let Some((nb, ii_b, _)) = second_form_at(sb, &point, Some(frame)) else {
                continue;
            };
            let sign = if na.dot(&nb) < 0.0 { -1.0 } else { 1.0 };
            let normal_sine = na.cross(&nb).length();
            let form_jump = spectral_norm(
                ii_a[0] - sign * ii_b[0],
                ii_a[1] - sign * ii_b[1],
                ii_a[2] - sign * ii_b[2],
            );
            let compatible_order = if normal_sine > gates.g1_normal {
                0
            } else if form_jump > gates.g2_curvature {
                1
            } else {
                2
            };
            out.push(SetbackCornerCompatibility {
                vertex,
                point,
                neighbours: (a.0.min(b.0), a.0.max(b.0)),
                normal_sine,
                form_jump,
                compatible_order,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spectral_norm_reads_the_worst_direction() {
        assert!((spectral_norm(0.0, 1.0, 0.0) - 1.0).abs() < 1e-15);
        assert!((spectral_norm(2.0, 0.0, -3.0) - 3.0).abs() < 1e-15);
        assert!((spectral_norm(1.0, 0.0, 1.0) - 1.0).abs() < 1e-15);
    }
}
