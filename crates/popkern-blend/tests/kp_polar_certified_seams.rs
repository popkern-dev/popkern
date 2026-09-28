//! PK-499 certified seam bounds on the KP polar setback fan (the
//! pentagon5_uniform fixture of `kp_polar_commercial_gates`).
//!
//! The sampled gates bound these seams at 5e-2. The certifier bounds
//! them over the whole edge:
//! - plain bi-nonic cap↔cap seams: G1 below 1e-6, G2 below 1e-3
//!   (subdivision seams ~2e-8 / 3.7e-4; the five spokes 2.9e-7 /
//!   7.5e-4 since PK-582's spoke re-reconcile, which restores the
//!   spoke G1 relation after the rim passes: they read 1.2e-5 / 3.1e-3
//!   from PK-575 to PK-582);
//! - cylinder strip↔cap seams: G1 below 1e-5 (measured 7e-7), and a
//!   certified transverse curvature jump of 0.093 (0.52 before PK-575).
//!   This is the G1 fan's real cap↔strip G2 residual, a certified number
//!   rather than a sample;
//! - dual-jet corner pieces (PK-541), certified through the patch's
//!   polynomial boundary surrogate:
//!   - corner-piece↔cap spokes: G1 below 1e-6 and G2 below 1e-3;
//!   - corner-piece↔strip rims: G1 below 1e-5. This is a G1 request,
//!     so the corner claims are not G2-reconciled, and the corner piece
//!     turns umbilic approaching its setback corner. Its cross
//!     curvature climbs from 1.5 at a tenth of the seam to the strip's
//!     1/r = 20 at the corner, where the cylinder reads 0. That jump is
//!     certified from both sides. The G2 request serves G2 there
//!     (`setback_rim_g2_evidence_gate`).

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

/// A cap face's ROLE: the setback corner piece is the one whose
/// boundary meets TWO strips (its rim sides meet at the setback corner).
/// Every other piece is a plain cap piece, whatever its patch variant.
/// Since PK-575 the pieces at spoke junctions are dual-jet patches too.
fn side_of(brep: &BRep, face: FaceHandle) -> Side {
    let is_strip = |f: FaceHandle| {
        matches!(
            brep.faces[f].surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        )
    };
    match brep.faces[face].surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(_)))
        | Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(_))) => {
            let strips = brep
                .edges
                .values()
                .filter(|e| {
                    let he = &brep.half_edges[e.half_edge];
                    let (fa, fb) = (
                        brep.loops[he.loop_].face,
                        brep.loops[brep.half_edges[he.twin].loop_].face,
                    );
                    (fa == face && is_strip(fb)) || (fb == face && is_strip(fa))
                })
                .count();
            if strips >= 2 {
                Side::DualJetCap
            } else {
                Side::PlainCap
            }
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
    let (mut cap_cap, mut strip_cap, mut corner_spoke, mut corner_rim) = (0, 0, 0, 0);
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
        let cert = certify_edge_seam(&brep, edge, &options);
        match sides {
            [Side::PlainCap, Side::PlainCap] => {
                cap_cap += 1;
                assert!(cert.is_certified(), "cap-cap {edge:?}: {cert:?}");
                // The five SPOKE seams carry the rim's junction
                // translation near their rim ends (the rim's transverse
                // curvature at a spoke junction is the spoke curve's own
                // normal curvature). PK-582 re-fits their cross rows to
                // the pre-pass G1 relation: certified g1 2.1-2.9e-7,
                // g2 6.9-7.5e-4 (1.2e-5 / 3.1e-3 without the re-fit).
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
                // from both sides. The rim reconciles' cluster fields
                // brought it from ~0.52 to 0.093 (PK-575).
                assert!(
                    cert.lower[2] > 0.05 && cert.g2_upper < 0.15,
                    "strip-cap {edge:?}: {cert:?}"
                );
            }
            [Side::DualJetCap, Side::PlainCap] | [Side::PlainCap, Side::DualJetCap] => {
                corner_spoke += 1;
                assert!(cert.is_certified(), "corner spoke {edge:?}: {cert:?}");
                assert!(
                    cert.g1_upper < 1e-6 && cert.g2_upper < 1e-3,
                    "corner spoke {edge:?}: {cert:?}"
                );
            }
            [Side::DualJetCap, Side::Strip] | [Side::Strip, Side::DualJetCap] => {
                corner_rim += 1;
                assert!(cert.is_certified(), "corner rim {edge:?}: {cert:?}");
                assert!(cert.g1_upper < 1e-5, "corner rim {edge:?}: {cert:?}");
                // The umbilic setback corner of a G1 corner piece: the
                // 1/r jump, attained within 1e-4 of a seam end.
                let at = cert.worst_at[2];
                assert!(
                    cert.lower[2] > 19.9 && cert.g2_upper < 20.2,
                    "corner rim {edge:?}: {cert:?}"
                );
                assert!(
                    at.1 <= 1e-4 || at.0 >= 1.0 - 1e-4,
                    "corner rim {edge:?}: {cert:?}"
                );
            }
            other => panic!("unexpected fan seam kinds {other:?} on {edge:?}"),
        }
    }
    assert_eq!(cap_cap, 20, "plain cap-cap seams");
    assert_eq!(strip_cap, 10, "plain cap-strip seams");
    assert_eq!(corner_spoke, 10, "dual-jet corner spokes");
    assert_eq!(corner_rim, 10, "dual-jet corner rims");
}
