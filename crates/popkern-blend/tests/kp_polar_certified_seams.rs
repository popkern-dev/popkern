//! PK-499 certified seam bounds on the KP polar setback fan (the
//! pentagon5_uniform fixture of `kp_polar_commercial_gates`).
//!
//! The sampled gates bound these seams at 5e-2. The certifier bounds
//! them over the whole edge:
//! - plain bi-nonic cap↔cap seams: G1 below 1e-6, G2 below 1e-3
//!   (measured 3.4e-7 and 4.0e-4);
//! - cylinder strip↔cap seams: G1 below 1e-5 (measured 3.2e-6), and a
//!   certified transverse curvature jump of ~0.52 — the G1 fan's real
//!   cap↔strip G2 residual, now a certified number rather than a
//!   sample;
//! - seams on the dual-jet corner pieces are the recorded frontier
//!   (blended corners the enclosures decline at, and a degenerate
//!   chart along the corner spokes): counted, not certified here.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::edge_seam_certificate::{certify_edge_seam, seam_certify_options};
use popkern_core::topology::{BRep, FaceHandle};

mod kp_gate_helpers;
use kp_gate_helpers::{fillet_corner_full, incident_edges, pick_vertex_with_valence};

#[derive(Debug, Clone, Copy, PartialEq)]
enum Side {
    PlainCap,
    DualJetCap,
    Strip,
    Other,
}

fn side_of(brep: &BRep, face: FaceHandle) -> Side {
    match brep.faces[face].surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(_))) => Side::PlainCap,
        Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(_))) => {
            Side::DualJetCap
        }
        Some(SurfaceGeometry::Cylinder { .. }) => Side::Strip,
        _ => Side::Other,
    }
}

#[test]
fn pentagon5_fan_seams_certify_over_the_whole_edge() {
    let pyr = BRep::n_gonal_pyramid(5, 1.0, 0.8, Tolerance::default()).expect("pyramid");
    let vertex = pick_vertex_with_valence(&pyr, 5).expect("apex");
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner_full(&pyr, &edges, &[0.05; 5], true, true);
    let options = seam_certify_options(&brep);
    let (mut cap_cap, mut strip_cap, mut frontier) = (0, 0, 0);
    for (edge, data) in brep.edges.iter() {
        let he = &brep.half_edges[data.half_edge];
        let (fa, fb) = (
            brep.loops[he.loop_].face,
            brep.loops[brep.half_edges[he.twin].loop_].face,
        );
        if fa == fb || !(faces.contains(&fa) || faces.contains(&fb)) {
            continue;
        }
        let sides = [side_of(&brep, fa), side_of(&brep, fb)];
        if sides.contains(&Side::DualJetCap) {
            frontier += 1;
            continue;
        }
        let cert = certify_edge_seam(&brep, edge, &options);
        match sides {
            [Side::PlainCap, Side::PlainCap] => {
                cap_cap += 1;
                assert!(cert.is_certified(), "cap-cap {edge:?}: {cert:?}");
                assert!(
                    cert.g1_upper < 1e-6 && cert.g2_upper < 1e-3,
                    "cap-cap {edge:?}: {cert:?}"
                );
            }
            [Side::PlainCap, Side::Strip] | [Side::Strip, Side::PlainCap] => {
                strip_cap += 1;
                assert!(cert.is_certified(), "strip-cap {edge:?}: {cert:?}");
                assert!(cert.g1_upper < 1e-5, "strip-cap {edge:?}: {cert:?}");
                // The G1 fan's real cap↔strip curvature jump, certified
                // from both sides.
                assert!(
                    cert.lower[2] > 0.4 && cert.g2_upper < 0.6,
                    "strip-cap {edge:?}: {cert:?}"
                );
            }
            other => panic!("unexpected fan seam kinds {other:?} on {edge:?}"),
        }
    }
    assert_eq!(cap_cap, 20, "plain cap-cap seams");
    assert_eq!(strip_cap, 10, "plain cap-strip seams");
    assert_eq!(
        frontier, 20,
        "dual-jet corner seams (the recorded frontier)"
    );
}
