//! Calibration tests for point↔mesh queries: the mesh-based SDF must
//! reproduce closed-form signed distances on analytic bodies to within
//! tessellation chord error. This is the trust anchor for the offset
//! SDF oracle — if these drift, every oracle verdict is suspect.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::PointMeshIndex;
use crate::TessellationOptions;

fn index_of(brep: &BRep) -> PointMeshIndex {
    let mesh = crate::tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    PointMeshIndex::build(&mesh)
}

/// Analytic SDF of the axis-aligned box `[lo, hi]`.
fn box_sdf(p: Point3, lo: Point3, hi: Point3) -> f64 {
    let cx = (lo.x - p.x).max(p.x - hi.x);
    let cy = (lo.y - p.y).max(p.y - hi.y);
    let cz = (lo.z - p.z).max(p.z - hi.z);
    let ox = cx.max(0.0);
    let oy = cy.max(0.0);
    let oz = cz.max(0.0);
    let outside = (ox * ox + oy * oy + oz * oz).sqrt();
    let inside = cx.max(cy).max(cz).min(0.0);
    outside + inside
}

fn unit_cube() -> BRep {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    BRep::prism(
        &profile,
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .expect("cube prism")
}

#[test]
fn sphere_signed_distance_matches_analytic() {
    let brep = BRep::sphere(1.0, Tolerance::default()).expect("sphere");
    let index = index_of(&brep);
    // Chord error bound for the default tessellation.
    let band = 0.01;
    let samples = [
        (Point3::new(0.0, 0.0, 0.0), -1.0),
        (Point3::new(0.5, 0.0, 0.0), -0.5),
        (Point3::new(0.0, 1.5, 0.0), 0.5),
        (Point3::new(1.2, 1.2, 1.2), (3.0_f64).sqrt() * 1.2 - 1.0),
        (Point3::new(0.0, 0.0, -2.0), 1.0),
    ];
    for (p, expected) in samples {
        let got = index.signed_distance(p);
        assert!(
            (got - expected).abs() <= band,
            "sphere sdf at {p:?}: got {got}, expected {expected}"
        );
    }
}

#[test]
fn cube_signed_distance_matches_analytic() {
    let index = index_of(&unit_cube());
    let (lo, hi) = (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
    // Planar faces tessellate exactly; only vertex weld tolerance
    // separates mesh from analytic distances.
    let band = 1e-6;
    let mut worst = 0.0_f64;
    let n = 11;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let p = Point3::new(
                    -0.5 + 2.0 * i as f64 / (n - 1) as f64,
                    -0.5 + 2.0 * j as f64 / (n - 1) as f64,
                    -0.5 + 2.0 * k as f64 / (n - 1) as f64,
                );
                let got = index.signed_distance(p);
                let expected = box_sdf(p, lo, hi);
                worst = worst.max((got - expected).abs());
            }
        }
    }
    assert!(worst <= band, "worst cube sdf deviation {worst}");
}

#[test]
fn containment_parity_on_cylinder() {
    let brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    let index = index_of(&brep);
    assert!(index.contains(Point3::new(0.0, 0.0, 1.0)));
    assert!(index.contains(Point3::new(0.9, 0.0, 0.1)));
    assert!(!index.contains(Point3::new(1.1, 0.0, 1.0)));
    assert!(!index.contains(Point3::new(0.0, 0.0, 2.1)));
    assert!(!index.contains(Point3::new(0.0, 0.0, -0.1)));
}

#[test]
fn empty_mesh_reports_infinite_distance() {
    let index = PointMeshIndex::build(&crate::TriMesh::default());
    assert!(index.distance(Point3::new(0.0, 0.0, 0.0)).is_infinite());
    assert!(!index.contains(Point3::new(0.0, 0.0, 0.0)));
}
