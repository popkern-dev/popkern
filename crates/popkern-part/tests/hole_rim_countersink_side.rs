//! Countersinks on hole rims at every side of a plate's middle.
//!
//! A chamfer's offset on a planar face must walk into that face's
//! material. The sign was calibrated toward the face's boundary-sample
//! centroid, which for a hole rim is the plate's middle: true only
//! when the rim's midpoint station (angle π, its −x side) faces it.
//! On the project box lid the two x < 0 holes were offset INTO the
//! hole, the bevel cone was built inverted (apex above the face,
//! opening down) and capped the hole to a pinhole while every
//! validator passed: the chamfer ADDED material. Holes all around the
//! middle put the midpoint station on both sides of the centroid.
//!
//! Oracle: each rim chamfer of `d` removes the ring between a 45°
//! frustum (r → r + d over depth d) and the bore, exactly.
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_part::{EdgeQuery, Part, Plane, Profile};

const R: f64 = 0.04;
const D: f64 = 0.03;
const T: f64 = 0.08;

/// Hole centres all around the plate's middle, the corners included.
const CENTRES: [(f64, f64); 8] = [
    (-0.61, -0.31),
    (0.61, -0.31),
    (0.61, 0.31),
    (-0.61, 0.31),
    (-0.61, 0.0),
    (0.61, 0.0),
    (0.0, -0.31),
    (0.0, 0.31),
];

fn drilled_plate() -> Part {
    let mut plate = Profile::on(Plane::XY)
        .rect(1.6, 1.0)
        .extrude(T)
        .expect("plate");
    for (x, y) in CENTRES {
        let bore = Profile::on(Plane::XY)
            .circle(R)
            .extrude(T + 0.2)
            .expect("bore")
            .translate(x, y, -0.1)
            .expect("bore placement");
        plate = plate.subtract(&bore).expect("drill");
    }
    plate
}

fn signed_mesh_volume(part: &Part) -> f64 {
    // Fine chords: the rims are 0.04 circles, and the default
    // chord height loses ~5% of the countersink rings to faceting.
    let options = popkern_mesh::TessellationOptions {
        chord_height: 1e-6,
        max_subdivisions: 512,
        ..Default::default()
    };
    let mesh = part.mesh(&options).expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

#[test]
fn every_hole_rim_countersinks_into_the_plate_whatever_its_side() {
    let plate = drilled_plate();
    let sunk = plate
        .chamfer(D, |e: EdgeQuery| e.circular().at_height(T))
        .expect("countersinks");

    // Every bevel is a cone opening UP out of the plate: apex below
    // the top face, on the hole's axis.
    let mut apexes = Vec::new();
    for face in sunk.brep().faces.values() {
        if let Some(SurfaceGeometry::Cone { apex, axis, .. }) =
            face.surface.as_ref().map(|s| &s.geometry)
        {
            apexes.push((*apex, *axis));
        }
    }
    assert_eq!(apexes.len(), CENTRES.len(), "one cone per rim");
    for (apex, axis) in &apexes {
        assert!(
            axis.z > 0.999,
            "cone at ({:.2}, {:.2}) opens along {axis:?}, not up out of the plate",
            apex.x,
            apex.y
        );
        assert!(
            (apex.z - (T - D - R)).abs() < 1e-9,
            "cone at ({:.2}, {:.2}): apex z {} (want {})",
            apex.x,
            apex.y,
            apex.z,
            T - D - R
        );
    }

    // Each countersink removes the frustum ring around its bore.
    let pi = std::f64::consts::PI;
    let (r0, r1) = (R, R + D);
    let ring = pi * D / 3.0 * (r0 * r0 + r0 * r1 + r1 * r1) - pi * r0 * r0 * D;
    let removed = signed_mesh_volume(&plate) - signed_mesh_volume(&sunk);
    let want = ring * CENTRES.len() as f64;
    assert!(
        (removed - want).abs() < 1e-3 * want,
        "countersinks removed {removed:.4e}, the frustum rings {want:.4e}"
    );
    assert_eq!(sunk.open_edge_count().expect("mesh"), 0, "watertight");
}
