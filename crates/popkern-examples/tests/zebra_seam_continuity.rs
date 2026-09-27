//! PK-503: the surface-evaluated zebra view reads seam continuity
//! class from the image. Across a plane→round seam the reflection
//! coordinate `s(x)` along an image row:
//! - breaks (value jump) at a G0 crease — the unrounded cube edge;
//! - kinks (slope jump, no value jump) at a G1-only seam — a circular
//!   round meets its planes with a curvature step 0 → 1/r;
//! - neither at a G2 seam — the curvature-continuous round.
//!
//! The surface view evaluates each face's own normal at the point the
//! pixel's ray hits, so these readings belong to the geometry, not to
//! the display mesh.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::zebra_attribution::summarize;
use popkern_viz::camera::Camera;
use popkern_viz::zebra::{ZebraNormals, ZebraSeamCrossing, rasterize_zebra_frame};

/// The cube's top edge along y at x = 1, z = 1.
fn top_edge(b: &BRep) -> popkern_core::topology::EdgeHandle {
    b.edges
        .iter()
        .find(|(_, e)| {
            let he = &b.half_edges[e.half_edge];
            let p = b.vertices[he.origin].position;
            let q = b.vertices[b.half_edges[he.twin].origin].position;
            (p.x - 1.0).abs() < 1e-9
                && (q.x - 1.0).abs() < 1e-9
                && (p.z - 1.0).abs() < 1e-9
                && (q.z - 1.0).abs() < 1e-9
        })
        .map(|(h, _)| h)
        .expect("top +x edge")
}

fn body(continuity: Option<ContinuityLevel>) -> BRep {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let Some(continuity) = continuity else {
        return cube;
    };
    let edge = top_edge(&cube);
    let options = FilletOptions {
        continuity,
        ..FilletOptions::default()
    };
    fillet_edges(
        &cube,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.3),
        }],
        &options,
        tol,
        None,
    )
    .expect("fillet")
    .brep
}

/// Looking down the rounded edge's bisector with y up: the edge runs
/// vertically through the image and every row crosses top → round →
/// side.
fn camera() -> Camera {
    Camera {
        eye: glam::Vec3::new(3.2, 0.5, 3.0),
        target: glam::Vec3::new(0.85, 0.5, 0.85),
        up: glam::Vec3::Y,
        fov_y: 0.35,
        aspect: 1.0,
        near: 0.01,
        far: 100.0,
    }
}

/// Worst-seam median value and slope jumps (per-seam statistics: the
/// round's two rims are different populations), and the crossing
/// count.
fn seam_reading(brep: &BRep, normals: ZebraNormals) -> (f64, f64, usize) {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    let frame = rasterize_zebra_frame(brep, &mesh, &camera(), 320, 320);
    let crossings: Vec<ZebraSeamCrossing> = frame.seam_crossings(normals);
    assert!(
        crossings.len() > 100,
        "rows must cross seams ({} crossings)",
        crossings.len()
    );
    let s = summarize(&crossings);
    assert!(s.seams >= 1);
    (s.value_p50, s.slope_p50, crossings.len())
}

#[test]
fn surface_zebra_reads_crease_kink_and_smooth_seams() {
    let crease = seam_reading(&body(None), ZebraNormals::Surface);
    let g1 = seam_reading(&body(Some(ContinuityLevel::G1)), ZebraNormals::Surface);
    let g2 = seam_reading(&body(Some(ContinuityLevel::G2)), ZebraNormals::Surface);
    eprintln!("crease {crease:?}\ng1 {g1:?}\ng2 {g2:?}");
    // G0 crease: the stripes break.
    assert!(crease.0 > 0.2, "crease must break the stripes: {crease:?}");
    // G1-only: continuous stripes that kink.
    assert!(g1.0 < 0.1 * crease.0, "G1 seam must not break: {g1:?}");
    // Measured 2026-09-27 (worst-seam medians, per pixel): G1 5.6e-3,
    // G2 1.2e-4 — a 48x separation.
    assert!(
        g1.1 > 20.0 * g2.1,
        "G1 seam must kink where G2 does not: g1 {g1:?} vs g2 {g2:?}"
    );
    // G2: neither.
    assert!(g2.0 < 0.1 * crease.0, "G2 seam must not break: {g2:?}");
}

/// The display mesh's view of the same G2 seam: its slope reading is
/// the tessellation's, recorded so the attribution has both numbers.
#[test]
fn mesh_zebra_on_the_g2_seam_is_recorded_against_the_surface_view() {
    let brep = body(Some(ContinuityLevel::G2));
    let surface = seam_reading(&brep, ZebraNormals::Surface);
    let mesh = seam_reading(&brep, ZebraNormals::Mesh);
    eprintln!("g2 seam: surface {surface:?} mesh {mesh:?}");
    assert_eq!(surface.2, mesh.2, "same frame, same crossings");
    // Measured 2026-09-27: surface 1.2e-4, mesh 1.8e-3 per pixel. The
    // excess kink in the mesh view is the display tessellation's.
    assert!(
        mesh.1 > 5.0 * surface.1,
        "surface {surface:?} mesh {mesh:?}"
    );
}
