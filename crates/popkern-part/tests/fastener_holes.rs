//! PK-460 receiving-hole serves: tapped holes (straight and ANGLED),
//! clearance holes with countersink/counterbore, external threads on
//! a shank, and the screw-fits-its-hole disjointness gate.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_part::Part;
use popkern_part::fastener::{
    ClearanceHoleOptions, HoleDepth, TapOptions, ThreadSpec, ThreadedRod,
};
use popkern_part::prelude::*;

fn block(w: f64, d: f64, h: f64) -> Part {
    Profile::on(Plane::XY).rect(w, d).extrude(h).expect("block")
}

/// Ray-parity point-in-mesh: count +x crossings of the mesh triangles.
fn inside_mesh(mesh: &popkern_mesh::TriMesh, p: Point3) -> bool {
    let mut crossings = 0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let b = mesh.vertices[t[1]].position;
        let c = mesh.vertices[t[2]].position;
        // Ray p + s·(1, 3e-7, 7e-7) (jittered to dodge edge hits).
        let d = Vec3::new(1.0, 3.0e-7, 7.0e-7);
        let e1 = b - a;
        let e2 = c - a;
        let pv = d.cross(&e2);
        let det = e1.dot(&pv);
        if det.abs() < 1e-14 {
            continue;
        }
        let inv = 1.0 / det;
        let tv = p - a;
        let u = tv.dot(&pv) * inv;
        if !(0.0..=1.0).contains(&u) {
            continue;
        }
        let qv = tv.cross(&e1);
        let v = d.dot(&qv) * inv;
        if v < 0.0 || u + v > 1.0 {
            continue;
        }
        let s = e2.dot(&qv) * inv;
        if s > 0.0 {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

#[test]
fn tapped_through_hole_serves_and_rod_fits_with_clearance() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let plate = block(14.0, 14.0, 3.5);
    let tapped = plate
        .tapped_hole(
            &spec,
            Point3::new(0.0, 0.0, 3.5),
            Vec3::new(0.0, 0.0, -1.0),
            HoleDepth::Through,
            &TapOptions {
                clearance: 0.2,
                entry_chamfer: true,
            },
        )
        .expect("tapped hole");
    assert_eq!(tapped.open_edge_count().expect("mesh"), 0, "watertight");
    // The screw-fits gate: a nominal rod on the same axis must be
    // fully disjoint from the tapped plate (positive modeled
    // clearance everywhere). The rod threads run off both ends so
    // the whole engaged stretch is thread-true.
    let rod = ThreadedRod::new(spec, 5.25)
        .tip_chamfer(false)
        .build(Tolerance::default())
        .expect("rod")
        .translate(0.0, 0.0, -1.0)
        .expect("place rod");
    let plate_mesh = tapped
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("plate mesh");
    let rod_mesh = rod
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("rod mesh");
    let mut inside = 0;
    for v in &rod_mesh.vertices {
        if inside_mesh(&plate_mesh, v.position) {
            inside += 1;
        }
    }
    assert_eq!(
        inside, 0,
        "{inside} rod surface vertices penetrate the tapped plate \
         (fit clearance violated)"
    );
}

#[test]
#[ignore = "blocked on the PK-460 walk-class campaign: stranded/loop-gap/LH/two-start merge defects (fork in flight)"]
fn angled_tapped_hole_serves_watertight() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let plate = block(20.0, 20.0, 8.0);
    // A non-axis-aligned fitting: 30° off vertical, tilted in x and y.
    let direction = Vec3::new(0.35, 0.2, -1.0);
    let tapped = plate
        .tapped_hole(
            &spec,
            Point3::new(1.0, -0.5, 8.0),
            direction,
            HoleDepth::Blind(4.5),
            &TapOptions {
                clearance: 0.15,
                entry_chamfer: true,
            },
        )
        .expect("angled tapped hole");
    assert_eq!(tapped.open_edge_count().expect("mesh"), 0, "watertight");
}

#[test]
fn clearance_holes_with_countersink_and_counterbore_serve() {
    let spec = ThreadSpec::parse("M5").expect("M5");
    let plate = block(30.0, 16.0, 6.0);
    let with_csk = plate
        .clearance_hole(
            &spec,
            Point3::new(-8.0, 0.0, 6.0),
            Vec3::new(0.0, 0.0, -1.0),
            HoleDepth::Through,
            &ClearanceHoleOptions {
                countersink: Some(0.2),
                ..ClearanceHoleOptions::default()
            },
        )
        .expect("countersunk clearance hole");
    let both = with_csk
        .clearance_hole(
            &spec,
            Point3::new(8.0, 0.0, 6.0),
            Vec3::new(0.0, 0.0, -1.0),
            HoleDepth::Through,
            &ClearanceHoleOptions {
                counterbore: Some(0.25),
                ..ClearanceHoleOptions::default()
            },
        )
        .expect("counterbored clearance hole");
    assert_eq!(both.open_edge_count().expect("mesh"), 0, "watertight");
}

#[test]
#[ignore = "blocked on the PK-460 walk-class campaign: stranded/loop-gap/LH/two-start merge defects (fork in flight)"]
fn external_thread_on_existing_shank_serves() {
    let spec = ThreadSpec::parse("M8x1").expect("M8 fine");
    let shank = Profile::on(Plane::XY)
        .circle(4.0)
        .extrude(10.0)
        .expect("shank");
    let threaded = shank
        .external_thread(
            &spec,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            6.0,
            10.0,
        )
        .expect("external thread");
    assert_eq!(threaded.open_edge_count().expect("mesh"), 0, "watertight");
}
