use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::modify_chamfer_inplace::modify_chamfer_distance;
use crate::types::{ChamferKind, ChamferOptions, ChamferSpec};

/// Chamfer one edge of a unit cube and return (result, chamfer face).
fn chamfered_cube(dist: f64) -> (BRep, popkern_core::topology::FaceHandle) {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let out = crate::chamfer_edges(
        &cube,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(dist),
        }],
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap()
    .brep;
    // The bevel is the one face whose plane normal is diagonal
    // (two equal-magnitude components) — curvature-analysis
    // identification alone classifies EVERY plane flanked by
    // non-tangent planes as a chamfer, so it cannot select the strip
    // among the cube's original faces.
    let face = out
        .faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                if let popkern_core::geometry::SurfaceGeometry::Plane { normal, .. } = &sof.geometry
                {
                    let mut c = [normal.x.abs(), normal.y.abs(), normal.z.abs()];
                    c.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                    c[0] < 1.0e-9 && (c[1] - c[2]).abs() < 1.0e-9 && c[1] > 0.1
                } else {
                    false
                }
            })
        })
        .map(|(f, _)| f)
        .expect("chamfered cube carries a diagonal bevel face");
    (out, face)
}

/// In-place distance change matches a fresh chamfer at the target
/// distance: same validity, and every strip corner within tolerance
/// of the fresh build's strip corners (the Hausdorff-style
/// equivalence the fillet in-place suite uses).
#[test]
fn equal_distance_change_matches_fresh_chamfer() {
    let (brep, face) = chamfered_cube(0.1);
    let modified =
        modify_chamfer_distance(&brep, face, ChamferKind::Equal(0.2), Tolerance::default())
            .expect("in-place chamfer distance change serves");
    modified.validate().expect("modified body validates");

    let (fresh, fresh_face) = chamfered_cube(0.2);
    let corner_positions = |b: &BRep, f: popkern_core::topology::FaceHandle| -> Vec<[f64; 3]> {
        let mut pts: Vec<[f64; 3]> = b
            .face_outer_loop_positions(f)
            .into_iter()
            .map(|p| [p.x, p.y, p.z])
            .collect();
        pts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        pts
    };
    let got = corner_positions(&modified, face);
    let want = corner_positions(&fresh, fresh_face);
    assert_eq!(got.len(), want.len(), "same strip corner count");
    for (g, w) in got.iter().zip(want.iter()) {
        for k in 0..3 {
            assert!(
                (g[k] - w[k]).abs() < 1.0e-9,
                "in-place corner {g:?} deviates from fresh-chamfer corner {w:?}"
            );
        }
    }
    // Topology untouched: same face/edge/vertex counts as before the
    // modification (the whole point of the in-place path).
    assert_eq!(modified.faces.len(), brep.faces.len());
    assert_eq!(modified.edges.len(), brep.edges.len());
    assert_eq!(modified.vertices.len(), brep.vertices.len());
}

/// Shrinking works too, and round-trips back to the original corners.
#[test]
fn distance_round_trip_restores_corners() {
    let (brep, face) = chamfered_cube(0.15);
    let bigger =
        modify_chamfer_distance(&brep, face, ChamferKind::Equal(0.25), Tolerance::default())
            .expect("grow serves");
    let back = modify_chamfer_distance(
        &bigger,
        face,
        ChamferKind::Equal(0.15),
        Tolerance::default(),
    )
    .expect("shrink back serves");
    let orig: Vec<_> = brep.face_outer_loop_positions(face);
    let restored: Vec<_> = back.face_outer_loop_positions(face);
    for (a, b) in orig.iter().zip(restored.iter()) {
        assert!(
            a.distance_to(b) < 1.0e-9,
            "round-trip corner drifted: {a:?} vs {b:?}"
        );
    }
}

/// A fillet face declines with the structured kind reason (the
/// chamfer path must not silently reshape fillets).
#[test]
fn fillet_face_declines_with_kind_reason() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let out = crate::fillet_edges(
        &cube,
        &[crate::FilletSpec {
            edge,
            radius: crate::FilletRadius::Constant(0.1),
        }],
        &crate::FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap()
    .brep;
    let fillet_face = out
        .faces
        .keys()
        .find(|&f| {
            crate::identify::identify_blend(&out, f)
                .is_some_and(|b| b.kind == crate::identify::BlendKind::Fillet)
        })
        .expect("fillet face identifiable");
    let err = modify_chamfer_distance(
        &out,
        fillet_face,
        ChamferKind::Equal(0.2),
        Tolerance::default(),
    )
    .expect_err("fillet face must decline");
    assert!(
        err.to_string().contains("Chamfer"),
        "decline names the accepted kind: {err}"
    );
}

/// DistanceAngle declines with the named remove+recreate route.
#[test]
fn distance_angle_declines_named() {
    let (brep, face) = chamfered_cube(0.1);
    let err = modify_chamfer_distance(
        &brep,
        face,
        ChamferKind::DistanceAngle {
            face: brep.faces.keys().next().unwrap(),
            distance: 0.2,
            angle: 0.5,
        },
        Tolerance::default(),
    )
    .expect_err("DistanceAngle declines");
    assert!(err.to_string().contains("chamfer_edges"), "{err}");
}
