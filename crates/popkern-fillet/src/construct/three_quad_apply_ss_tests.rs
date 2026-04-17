//! §5.3 R8 C8e — Shirman-Séquin central-vertex G2 test.
//!
//! Companion to `three_quad_apply_g2_tests.rs` (the §2.4 boundary
//! scorecard). Split out to keep the G2 test file under the 500-line
//! structural cap.

use crate::types::{FilletOptions, FilletRadius, FilletSpec};

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// R8 C8e — Shirman-Séquin central-vertex G2 vertex-enclosure
/// property. Runs the mixed-radii Gregory dispatch and verifies
/// that the three sub-patches' mixed second partials `τ_k` at the
/// shared central vertex agree in the **tangent-plane projection**
/// (the SS "relaxed G2" condition).
///
/// Before SS lands (C4b.2b common-plane projection only), the
/// tangent-plane components of `τ_k` differ across sub-patches,
/// reflecting the independent per-patch parallelogram offsets.
/// SS forces them to share a common tangent-plane component —
/// visible as `max_pairwise_diff < 1e-9` after the solve.
#[test]
fn r8_c8e_central_vertex_shirman_sequin_tangent_agreement() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let tol = Tolerance::default();
    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let pipeline = crate::fillet_edges(&cube, &specs, &opts, tol, None)
        .expect("fillet pipeline must succeed on mixed-radii cube corner");
    let corner_faces: Vec<_> = pipeline
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .map(|p| p.face)
        .collect();
    assert_eq!(
        corner_faces.len(),
        3,
        "mixed-radii dispatch must produce 3 Gregory sub-patches"
    );

    // Extract each sub-patch's Chiyokura-Kimura mixed second partial
    // `τ` at its corner-3 (the shared central vertex). For a cubic
    // Bezier at corner (1, 0): τ = 9·(P_20 + P_31 − P_30 − P_21).
    let mut taus: Vec<Vec3> = Vec::with_capacity(3);
    let mut central_point: Option<Point3> = None;
    for &fh in &corner_faces {
        let face = pipeline.brep.faces.get(fh).unwrap();
        let sof = face.surface.as_ref().unwrap();
        let rect = match &sof.geometry {
            SurfaceGeometry::Gregory(GregoryPatch::Rectangular(r)) => r,
            _ => panic!("expected RectangularGregory sub-patch"),
        };
        let p20 = rect.boundary_cps[2][0];
        let p31 = rect.boundary_cps[3][1];
        let p30 = rect.boundary_cps[3][0];
        let p21 = rect.inner_cps[3].0;
        let tau = Vec3::new(
            9.0 * (p20.x + p31.x - p30.x - p21.x),
            9.0 * (p20.y + p31.y - p30.y - p21.y),
            9.0 * (p20.z + p31.z - p30.z - p21.z),
        );
        taus.push(tau);
        central_point = Some(match central_point {
            None => p30,
            Some(prev) => {
                assert!(
                    (prev.x - p30.x).abs() < 1e-9
                        && (prev.y - p30.y).abs() < 1e-9
                        && (prev.z - p30.z).abs() < 1e-9,
                    "all 3 sub-patches must share the same central vertex"
                );
                prev
            }
        });
    }

    // Reconstruct the 6 hex vertices from the sub-patch corners and
    // recompute the hex-plane normal exactly the way the SS code does
    // — this keeps the test's tangent basis aligned with SS's, so the
    // tangent-plane agreement we're checking is the one the solve
    // guarantees (not a different "surface-normal-averaged" basis).
    let rects: Vec<_> = corner_faces
        .iter()
        .map(|&fh| {
            let face = pipeline.brep.faces.get(fh).unwrap();
            let sof = face.surface.as_ref().unwrap();
            match &sof.geometry {
                SurfaceGeometry::Gregory(GregoryPatch::Rectangular(r)) => r.as_ref().clone(),
                _ => panic!("expected RectangularGregory sub-patch"),
            }
        })
        .collect();
    let hex: [Point3; 6] = [
        rects[0].boundary_cps[0][0],
        rects[0].boundary_cps[0][3],
        rects[1].boundary_cps[0][0],
        rects[1].boundary_cps[0][3],
        rects[2].boundary_cps[0][0],
        rects[2].boundary_cps[0][3],
    ];
    let c_pt = central_point.unwrap();
    let mut n_sum = Vec3::new(0.0, 0.0, 0.0);
    for k in 0..6 {
        let a = Vec3::new(hex[k].x - c_pt.x, hex[k].y - c_pt.y, hex[k].z - c_pt.z);
        let b = Vec3::new(
            hex[(k + 1) % 6].x - c_pt.x,
            hex[(k + 1) % 6].y - c_pt.y,
            hex[(k + 1) % 6].z - c_pt.z,
        );
        let cr = a.cross(&b);
        n_sum = Vec3::new(n_sum.x + cr.x, n_sum.y + cr.y, n_sum.z + cr.z);
    }
    let n_c = n_sum
        .normalized()
        .expect("hex-plane normal should be non-degenerate for this fixture");

    // Tangent-plane component of each τ: τ − (τ·n_c) n_c.
    let project_to_plane = |v: &Vec3| -> Vec3 {
        let d = v.dot(&n_c);
        Vec3::new(v.x - d * n_c.x, v.y - d * n_c.y, v.z - d * n_c.z)
    };
    let tau_tan: Vec<Vec3> = taus.iter().map(project_to_plane).collect();

    // SS property: all 3 tangent-plane components agree.
    let mut max_diff = 0.0_f64;
    for i in 0..3 {
        for j in (i + 1)..3 {
            let d = Vec3::new(
                tau_tan[i].x - tau_tan[j].x,
                tau_tan[i].y - tau_tan[j].y,
                tau_tan[i].z - tau_tan[j].z,
            );
            max_diff = max_diff.max(d.length());
        }
    }
    let tau_scale = tau_tan.iter().map(|v| v.length()).fold(0.0_f64, f64::max);

    eprintln!(
        "R8 C8e Shirman-Séquin check (mixed radii {radii:?}):\n  \
         |τ_tan|_max = {tau_scale:.4e}\n  \
         max pairwise |τ_tan_k − τ_tan_k'| = {max_diff:.4e}"
    );

    assert!(
        max_diff < 1e-9,
        "R8 C8e: Shirman-Séquin tangent-plane agreement violated \
         (max pairwise diff {max_diff:.4e}, τ_tan scale {tau_scale:.4e})"
    );
}
