//! Quick validity check: does the committed cube_extreme binonic
//! output self-intersect? The promotion pipeline has a
//! self-intersection screen (at 1e-6 tolerance), but the twist
//! optimizer doesn't — so it's possible the optimizer returned
//! a self-intersecting patch that the screen accepted on a
//! different tolerance. Sanity-check.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{collect_patches, fillet_corner, incident_edges};

#[test]
fn cube_extreme_binonic_is_not_self_intersecting() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, &[0.02, 0.10, 0.20], true);
    let patches = collect_patches(&brep, &faces);
    for (i, p) in patches.iter().enumerate() {
        match p {
            // PK-615: the corner-cell tangent fit serves dual-jet sectors
            // here when it is fairer than the plain fit (26.8 against
            // 576 / r_min); they take the same check.
            GregoryPatch::RectangularBiNonic(_) | GregoryPatch::RectangularBiNonicDualJet(_) => {
                let self_int_1e6 = p.self_intersects_with_tolerance(1e-6);
                let self_int_1e4 = p.self_intersects_with_tolerance(1e-4);
                println!(
                    "patch {i}: binonic — self_int@1e-6={self_int_1e6}, self_int@1e-4={self_int_1e4}"
                );
                assert!(
                    !self_int_1e6,
                    "patch {i} self-intersects at 1e-6 tolerance — optimizer output is invalid"
                );
            }
            GregoryPatch::RectangularTwoStripe(_) => {
                println!("patch {i}: two-stripe (guard reverted)");
            }
            other => panic!("unexpected variant: {other:?}"),
        }
    }
}
