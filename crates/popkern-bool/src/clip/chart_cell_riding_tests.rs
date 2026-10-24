use super::super::{cell_windows, face_chart_cell};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A record along the meridian patch edge at u = 0.7 of a torus turned
/// 0.7 about its axis, its u wobbling +-3e-14 as a marched record's
/// does: exactly one of the edge's two faces claims it, whole.
#[test]
fn a_record_riding_a_patch_edge_is_claimed_by_one_face_whole() {
    let tol = Tolerance::default();
    let torus = BRep::torus(1.0, 0.3, tol)
        .and_then(|t| t.transform(&Transform3::rotation_z(0.7)))
        .expect("torus");
    let n = 65;
    let control_points: Vec<Point3> = (0..n)
        .map(|i| {
            let wobble = if i % 2 == 0 { 3e-14 } else { -3e-14 };
            Point3::new(0.7 + wobble, 1.7 + 1.3 * i as f64 / (n - 1) as f64, 0.0)
        })
        .collect();
    let mut knots = vec![0.0];
    knots.extend((0..n).map(|i| i as f64 / (n - 1) as f64));
    knots.push(1.0);
    let record = NurbsCurve {
        degree: 1,
        knots,
        control_points,
        weights: vec![1.0; n],
    };
    let mut claims = Vec::new();
    for (face, _) in &torus.faces {
        let Some(cell) = face_chart_cell(&torus, face, &tol) else {
            continue;
        };
        let windows = cell_windows(&cell, &record);
        if !windows.is_empty() {
            claims.push((cell.u_range, windows));
        }
    }
    assert_eq!(claims.len(), 1, "claims: {claims:?}");
    let (u_range, windows) = &claims[0];
    assert!(
        (u_range.0 - 0.7).abs() < 1e-12,
        "the low-u side claims: {u_range:?}"
    );
    assert_eq!(windows, &vec![(0.0, 1.0)]);
}
