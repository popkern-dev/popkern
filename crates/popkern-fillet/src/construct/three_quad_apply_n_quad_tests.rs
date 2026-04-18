//! §5.3 item #4 — N ≥ 4 corner dispatch fixture tests.
//!
//! Uses `BRep::square_pyramid` (which has a 4-valence apex with 4
//! planar triangular side faces) to validate the N-quad Gregory
//! pipeline end-to-end beyond the synthetic topology pokes in
//! `poke_hex_tests.rs`. A successful dispatch should:
//!
//! 1. Run `fillet_edges` without error on the 4 slant edges.
//! 2. Produce 4 Gregory sub-patches at the apex (one per sub-quad).
//! 3. Yield a valid BRep post-fillet.

use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CornerMethod, fillet_edges};

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// §5.3 item #4 — fillet all 4 slant edges of a square pyramid and
/// verify the apex routes through the N=4 three-quad (now "N-quad")
/// dispatch, producing 4 Gregory sub-patches.
#[test]
fn item4_square_pyramid_apex_produces_four_gregory_subpatches() {
    let tol = Tolerance::default();
    let pyramid = BRep::square_pyramid(1.0, 0.8, tol).expect("build square pyramid");

    // Identify the apex: 4 incident edges, none on the base plane.
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_h, v)| v.position.z > 0.1)
        .map(|(h, _)| h)
        .expect("apex (z > 0) must exist on the pyramid");

    // Collect the 4 slant edges emanating from the apex.
    let slant_edges: Vec<_> = pyramid
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    let slant_edges: Vec<_> = {
        // Dedupe (each edge visited twice — once per twin).
        let mut seen = std::collections::BTreeSet::new();
        slant_edges
            .into_iter()
            .filter(|e| seen.insert(*e))
            .collect()
    };
    assert_eq!(
        slant_edges.len(),
        4,
        "square pyramid apex should have 4 incident slant edges"
    );

    // Mixed radii: uniform radii would let the rolling-ball sphere
    // fit first and bypass the Gregory N-quad path. Varying radii
    // forces the sphere fit to fail and routes to Gregory.
    let radii = [0.07_f64, 0.09, 0.11, 0.13];
    let specs: Vec<FilletSpec> = slant_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();

    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&pyramid, &specs, &opts, tol, None)
        .expect("fillet pipeline must succeed on pyramid apex");

    // The apex should be replaced by 4 corner patches (one per quad).
    let corner_patches: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == apex)
        .collect();
    assert_eq!(
        corner_patches.len(),
        4,
        "N=4 dispatch should yield 4 corner patches at the apex; \
         got {}",
        corner_patches.len()
    );

    for (idx, patch) in corner_patches.iter().enumerate() {
        assert_eq!(patch.method, CornerMethod::Gregory);
        assert_eq!(patch.n_sides, 4);
        let face = result.brep.faces.get(patch.face).expect("sub-face");
        let sof = face.surface.as_ref().expect("sub-face has surface");
        let is_gregory_rect = matches!(
            &sof.geometry,
            SurfaceGeometry::Gregory(
                GregoryPatch::Rectangular(_) | GregoryPatch::RectangularDeg5(_)
            )
        );
        assert!(
            is_gregory_rect,
            "apex sub-patch {idx} must be a rectangular Gregory surface; \
             got {:?}",
            std::mem::discriminant(&sof.geometry)
        );
    }

    // Euler invariant check: N-gon poke adds 1 vertex, N edges,
    // N-1 faces. For N=4: +1V, +4E, +3F compared to the baseline
    // (polar-NURBS) pipeline.
    let baseline_opts = FilletOptions {
        use_rational_gregory: false,
        ..FilletOptions::default()
    };
    let baseline = fillet_edges(&pyramid, &specs, &baseline_opts, tol, None)
        .expect("baseline polar-NURBS pipeline must succeed");
    assert_eq!(
        result.brep.vertices.len(),
        baseline.brep.vertices.len() + 1,
        "N=4 dispatch: +1 central vertex"
    );
    assert_eq!(
        result.brep.edges.len(),
        baseline.brep.edges.len() + 4,
        "N=4 dispatch: +4 spoke edges"
    );
    assert_eq!(
        result.brep.faces.len(),
        baseline.brep.faces.len() + 3,
        "N=4 dispatch: +3 sub-faces (1 → 4)"
    );
}
