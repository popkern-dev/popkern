//! Cone-band dress-up gates: bodies whose side band is an analytic
//! cone authored by a rim chamfer or by drafting a cylinder wall.
//! Both classes historically shipped cracked-seam bands that meshed
//! open while the BRep validated.

use crate::{Dir, EdgeQuery, FaceQuery, Plane, Profile, SurfaceKind};

/// Vertical-ray parity containment against the default-quality mesh,
/// with the ray nudged off exact seam-station vertex columns (an
/// on-vertex ray double-counts its crossing).
fn contains(part: &crate::Part, q: (f64, f64, f64)) -> bool {
    let m = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let q = (q.0 + 1.03e-7, q.1 + 1.71e-7, q.2);
    let mut crossings = 0;
    for t in &m.triangles {
        let (a, b, c) = (
            m.vertices[t[0]].position,
            m.vertices[t[1]].position,
            m.vertices[t[2]].position,
        );
        let (e1x, e1y) = (b.x - a.x, b.y - a.y);
        let (e2x, e2y) = (c.x - a.x, c.y - a.y);
        let det = e1x * e2y - e1y * e2x;
        if det.abs() < 1e-14 {
            continue;
        }
        let (rx, ry) = (q.0 - a.x, q.1 - a.y);
        let u = (rx * e2y - ry * e2x) / det;
        let v = (e1x * ry - e1y * rx) / det;
        if u < -1e-12 || v < -1e-12 || u + v > 1.0 + 1e-12 {
            continue;
        }
        let z = a.z + u * (b.z - a.z) + v * (c.z - a.z);
        if z > q.2 {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

/// PIN → GATE (countersink-tool campaign): rim chamfers with LARGE
/// legs stay conical and watertight. The offset walk used chord
/// tangents (20 samples ⇒ ~1.45e-2·distance row spread on the cap,
/// where back-projection cannot absorb it); past ~2 mm the spread
/// crossed the 1e-3·r conical-recognition gate and the bevel shipped
/// as a cracked-seam NURBS ring (~2.8k open mesh edges, and the tool
/// then tripped the boolean's sheet-body pre-check). Closed by
/// sampling the EXACT curve tangent (derivative_at) with the chord
/// as degenerate-only fallback.
#[test]
fn large_leg_rim_chamfer_stays_conical_and_watertight() {
    for leg in [3.0, 8.0] {
        let cyl = Profile::on(Plane::XY)
            .circle(44.0)
            .extrude(20.0)
            .expect("cyl");
        let out = cyl
            .chamfer(leg, move |e: EdgeQuery| e.circular().at_height(20.0))
            .expect("chamfer");
        let opens = out.open_edge_count().expect("mesh");
        assert_eq!(opens, 0, "leg {leg}: cracked bevel band");
        assert_eq!(out.face_count(), 4, "leg {leg}");
        // Material probes: rim corner gone, cone flank present.
        assert!(!contains(&out, (43.8, 0.0, 19.8)));
        assert!(contains(&out, (44.0 - leg - 0.2, 0.0, 19.8)));
        assert!(contains(&out, (43.8, 0.0, 20.0 - leg - 0.2)));
    }
}

/// PIN → GATE (draft cone support): drafting a cylinder side wall
/// yields a WATERTIGHT exact frustum at any angle, both signs. Four
/// stacked defects previously left it open (638 mesh opens at every
/// angle while the BRep "validated"): the cone apex was offset from
/// the axis (placed above the neutral line's WALL point), the rim
/// circles and seam line kept their cylinder-era curves (the edge
/// recompute was Plane×Plane-only), rim vertices stayed at the old
/// radius, and the full-span cylinder-chart rim pcurves survived the
/// degenerate-span keep gate so fix_pcurve_endpoints collapsed them
/// to zero span (now cleared by the image-staleness check).
#[test]
fn drafted_cylinder_frustum_is_watertight_and_conical() {
    for angle_deg in [2.0_f64, 8.5, 45.0, -8.5] {
        let out = Profile::on(Plane::XY)
            .circle(38.0)
            .extrude(30.0)
            .expect("cyl")
            .draft(angle_deg.to_radians(), Dir::Z, |f: FaceQuery| {
                f.kind(SurfaceKind::Cylinder)
            })
            .expect("draft");
        let opens = out.open_edge_count().expect("mesh");
        assert_eq!(opens, 0, "{angle_deg}deg: cracked cone band");
        assert_eq!(out.face_count(), 3, "{angle_deg}deg");
        // Cone flank position: at the top (z 30) the radius is
        // 38 ∓ 30·tan|angle| (positive angle narrows toward +Z).
        let r_top = 38.0 - 30.0 * angle_deg.to_radians().tan();
        assert!(
            contains(&out, (r_top - 0.3, 0.0, 29.8)),
            "{angle_deg}deg: material missing inside the top rim"
        );
        assert!(
            !contains(&out, (r_top + 0.3, 0.0, 29.8)),
            "{angle_deg}deg: material outside the top rim"
        );
        // The neutral (bottom) circle is unchanged. Probe just above
        // the base with margins that clear the wall at any taper
        // (cone radius at z 0.02 lies within 38 ∓ 0.02·tan 45°).
        assert!(contains(&out, (37.5, 0.0, 0.02)));
        assert!(!contains(&out, (38.5, 0.0, 0.02)));
    }
}
