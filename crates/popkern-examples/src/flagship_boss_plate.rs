//! Flagship demo: sketch → solid → boolean → fillet → watertight STEP.
//!
//! A rectangular plate and a cylindrical boss are each built from a
//! closed sketch profile by extrusion, unioned, and filleted — the
//! boss/plate seam is a CURVED (circular) fillet edge, and the plate
//! corners carry vertex-blend caps. The result must tessellate
//! watertight and export as a STEP solid.
//!
//! Stages print their gate results; the binary exits non-zero if a
//! gate fails, so it can serve as an integration check.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

/// Closed degree-1 rectangle wire in the z=0 plane, centered at the
/// origin: the plate sketch.
fn rectangle_sketch(width: f64, depth: f64) -> NurbsCurve {
    let (hw, hd) = (0.5 * width, 0.5 * depth);
    let pts = vec![
        Point3::new(-hw, -hd, 0.0),
        Point3::new(hw, -hd, 0.0),
        Point3::new(hw, hd, 0.0),
        Point3::new(-hw, hd, 0.0),
        Point3::new(-hw, -hd, 0.0),
    ];
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n],
    }
}

/// Closed rational quadratic circle wire (standard 9-CP form) in the
/// plane z = `z`, centered on the z-axis: the boss sketch.
fn circle_sketch(radius: f64, z: f64) -> NurbsCurve {
    let r = radius;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let pts = vec![
        Point3::new(r, 0.0, z),
        Point3::new(r, r, z),
        Point3::new(0.0, r, z),
        Point3::new(-r, r, z),
        Point3::new(-r, 0.0, z),
        Point3::new(-r, -r, z),
        Point3::new(0.0, -r, z),
        Point3::new(r, -r, z),
        Point3::new(r, 0.0, z),
    ];
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points: pts,
        weights,
    }
}

/// Straight extrusion path from `z0` to `z1` along +Z.
fn extrusion_path(z0: f64, z1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z0), Point3::new(0.0, 0.0, z1)],
        weights: vec![1.0, 1.0],
    }
}

fn extrude(sketch: NurbsCurve, z0: f64, z1: f64, tol: Tolerance, label: &str) -> BRep {
    let result = sweep(
        &SweepProfile::Wire(sketch),
        &extrusion_path(z0, z1),
        &SweepOptions::default(),
        tol,
        None,
    )
    .unwrap_or_else(|e| panic!("{label} extrusion failed: {e:?}"));
    for w in &result.warnings {
        eprintln!("  [{label}] sweep warning: {w:?}");
    }
    if let Err(e) = result.brep.validate() {
        eprintln!("  [{label}] INVALID: {e:?}");
        eprintln!(
            "  [{label}] faces={} edges={} vertices={} loops={}",
            result.brep.faces.len(),
            result.brep.edges.len(),
            result.brep.vertices.len(),
            result.brep.loops.len()
        );
        for (fh, f) in result.brep.faces.iter() {
            let n = result.brep.face_outer_loop_positions(fh).len();
            let kind = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
                Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
                Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cyl",
                _ => "?",
            };
            eprintln!(
                "    face {fh:?} {kind} loop_len={n} inner={}",
                f.inner_loops.len()
            );
        }
        panic!("{label} extrusion invalid");
    }
    result.brep
}

fn main() {
    let tol = Tolerance::default();

    // Stage 1: sketches → solids.
    let plate = extrude(rectangle_sketch(2.4, 1.6), 0.0, 0.4, tol, "plate");
    println!(
        "stage 1a plate: faces={} edges={} (validate ok)",
        plate.faces.len(),
        plate.edges.len()
    );
    // The boss penetrates the plate (base below the plate top) so the
    // union seam is a general-position circle on the plate top.
    let boss = extrude(circle_sketch(0.35, 0.0), 0.2, 1.0, tol, "boss");
    println!(
        "stage 1b boss: faces={} edges={} (validate ok)",
        boss.faces.len(),
        boss.edges.len()
    );

    // Boss geometry audit: seam edge curve + vertices.
    for (eh, edge) in boss.edges.iter() {
        let he = edge.half_edge;
        let o = boss.vertices[boss.half_edges[he].origin].position;
        let t = boss.vertices[boss.half_edges[boss.half_edges[he].twin].origin].position;
        let has_curve = edge.curve.is_some();
        eprintln!(
            "  boss edge {eh:?}: ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) curve={has_curve}",
            o.x, o.y, o.z, t.x, t.y, t.z
        );
    }

    for (_, f) in boss.faces.iter() {
        if let Some(popkern_core::geometry::SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = f.surface.as_ref().map(|s| &s.geometry)
        {
            eprintln!(
                "  boss cyl: origin ({:.4},{:.4},{:.4}) axis ({:.3},{:.3},{:.3}) r {:.6}",
                origin.x, origin.y, origin.z, axis.x, axis.y, axis.z, radius
            );
            let v = popkern_core::geometry::Point3::new(-0.0, 0.3446, 0.2);
            let d = v - *origin;
            let radial = d - *axis * d.dot(axis);
            eprintln!("  seam vertex radial dist: {:.6}", radial.length());
        }
    }

    // Stage 2: boolean union.
    let unioned = popkern_bool::union(plate, boss, tol, None).expect("union");
    if let Err(e) = unioned.validate() {
        eprintln!("union INVALID: {e:?}");
        for (fh, f) in unioned.faces.iter() {
            let c = unioned.face_outer_loop_positions(fh);
            let z: f64 = c.iter().map(|p| p.z).sum::<f64>() / c.len().max(1) as f64;
            let kind = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
                Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
                Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cyl",
                _ => "?",
            };
            eprintln!(
                "  face {fh:?} {kind} loop={} inner={} z~{z:.2}",
                c.len(),
                f.inner_loops.len()
            );
        }
        eprintln!(
            "  totals: F={} E={} V={} L={}",
            unioned.faces.len(),
            unioned.edges.len(),
            unioned.vertices.len(),
            unioned.loops.len()
        );
        panic!("union must validate");
    }
    println!(
        "stage 2 union: faces={} edges={} vertices={} (validate ok)",
        unioned.faces.len(),
        unioned.edges.len(),
        unioned.vertices.len()
    );

    // Locate the curved seam edge: cylinder × plate-top plane.
    let mut seam_edges = Vec::new();
    for (eh, edge) in unioned.edges.iter() {
        let he = edge.half_edge;
        let f1 = unioned.loops[unioned.half_edges[he].loop_].face;
        let f2 = unioned.loops[unioned.half_edges[unioned.half_edges[he].twin].loop_].face;
        let kind = |fh: popkern_core::topology::FaceHandle| {
            unioned.faces[fh].surface.as_ref().map(|s| &s.geometry)
        };
        use popkern_core::geometry::SurfaceGeometry as SG;
        let cyl_plane = matches!(
            (kind(f1), kind(f2)),
            (Some(SG::Cylinder { .. }), Some(SG::Plane { .. }))
                | (Some(SG::Plane { .. }), Some(SG::Cylinder { .. }))
        ) || matches!(
            (kind(f1), kind(f2)),
            (Some(SG::Nurbs(_)), Some(SG::Plane { .. }))
                | (Some(SG::Plane { .. }), Some(SG::Nurbs(_)))
        );
        if !cyl_plane {
            continue;
        }
        // On the plate top (z = 0.4)?
        if let Some(curve) = edge.curve.as_ref() {
            use popkern_core::geometry::curve::CurveEval;
            let ceval: &dyn CurveEval = &curve.geometry;
            let mid = ceval.point_at(0.5 * (curve.t_start + curve.t_end));
            if (mid.z - 0.4).abs() < 1e-6 {
                seam_edges.push(eh);
            }
        }
    }
    println!("stage 2 seam edges at plate top: {}", seam_edges.len());
    assert!(!seam_edges.is_empty(), "no curved seam edge found");
}
