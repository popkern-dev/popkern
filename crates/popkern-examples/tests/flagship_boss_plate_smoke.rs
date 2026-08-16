//! Flagship boss-plate smoke gate (PK-294).
//!
//! The flagship example rotted invisibly for weeks — the sweep stayed
//! green while `cargo run --bin flagship_boss_plate` panicked (an
//! empty-curve split via the PnS poke normalization) and, once that
//! was fixed, tessellated 98 open weld edges (a periodic-loop winding
//! tear on the mitred seam sliver). Examples are not executed by the
//! test sweep, so this test distils the example's serving pipeline —
//! plate ∪ boss, the 9-edge fillet chain (1 circular seam + 8 plate
//! edges), tessellation — and gates the same outcomes the example
//! asserts: validation and watertight welding.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_sweep::{SweepOptions, SweepProfile, sweep};
use std::collections::HashMap;

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

fn circle_sketch(radius: f64, z: f64) -> NurbsCurve {
    let r = radius;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, z),
            Point3::new(r, r, z),
            Point3::new(0.0, r, z),
            Point3::new(-r, r, z),
            Point3::new(-r, 0.0, z),
            Point3::new(-r, -r, z),
            Point3::new(0.0, -r, z),
            Point3::new(r, -r, z),
            Point3::new(r, 0.0, z),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

fn extrude(sketch: NurbsCurve, z0: f64, z1: f64, tol: Tolerance) -> BRep {
    let path = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z0), Point3::new(0.0, 0.0, z1)],
        weights: vec![1.0, 1.0],
    };
    let result = sweep(
        &SweepProfile::Wire(sketch),
        &path,
        &SweepOptions::default(),
        tol,
        None,
    )
    .expect("extrusion");
    result.brep.validate().expect("extrusion valid");
    result.brep
}

#[test]
fn flagship_boss_plate_pipeline_serves_watertight() {
    let tol = Tolerance::default();
    let plate = extrude(rectangle_sketch(2.4, 1.6), 0.0, 0.4, tol);
    let boss = extrude(circle_sketch(0.35, 0.0), 0.2, 1.0, tol);
    let unioned = popkern_bool::union(plate, boss, tol, None).expect("union");
    unioned.validate().expect("union valid");

    // The seam edge (boss cylinder × plate top) + the 8 plate edges.
    let mut specs = Vec::new();
    for (eh, d) in unioned.edges.iter() {
        let he = &unioned.half_edges[d.half_edge];
        let f1 = unioned.loops[he.loop_].face;
        let f2 = unioned.loops[unioned.half_edges[he.twin].loop_].face;
        if f1 == f2 {
            continue; // cylinder seam spine — not filletable
        }
        let a = unioned.vertices[he.origin].position;
        let b = unioned.vertices[unioned.half_edges[he.twin].origin].position;
        // The flagship chain: the circular boss seam + the plate's
        // top and bottom rim edges (z = 0.4 and z = 0). The boss top
        // rim (z = 1.0) and the vertical plate edges are excluded.
        let top = (a.z - 0.4).abs() < 1e-9 && (b.z - 0.4).abs() < 1e-9;
        let bottom = a.z.abs() < 1e-9 && b.z.abs() < 1e-9;
        if !top && !bottom {
            continue;
        }
        specs.push(FilletSpec {
            edge: eh,
            radius: FilletRadius::Constant(0.08),
        });
    }
    assert_eq!(
        specs.len(),
        9,
        "flagship chain is 1 seam + 8 plate-top edges, found {}",
        specs.len()
    );

    let result = fillet_edges(&unioned, &specs, &FilletOptions::default(), tol, None)
        .expect("flagship fillet must complete");
    result.brep.validate().expect("filleted flagship validates");

    // Tessellate and weld-count exactly as the example does.
    let mesh =
        popkern_mesh::tessellate(&result.brep, &popkern_mesh::TessellationOptions::default())
            .expect("tessellation");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut canon: Vec<Point3> = Vec::new();
    let id = |p: Point3, key: &mut HashMap<(i64, i64, i64), usize>, canon: &mut Vec<Point3>| {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        *key.entry(k).or_insert_with(|| {
            canon.push(p);
            canon.len() - 1
        })
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position, &mut key, &mut canon),
            id(mesh.vertices[t[1]].position, &mut key, &mut canon),
            id(mesh.vertices[t[2]].position, &mut key, &mut canon),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let open = counts.values().filter(|&&c| c == 1).count();
    assert_eq!(open, 0, "flagship must tessellate watertight, {open} open");
}
