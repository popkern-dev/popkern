//! PK-632: the intersection of two toothed prisms whose caps are
//! COPLANAR (two mating gears turned off their mesh phase) serves,
//! closed, with the face count and volume OpenCASCADE's
//! `BRepAlgoAPI_Common` gives for the same solids (recorded below,
//! 2026-10-03). Three clip defects made it decline or serve open:
//!
//! - the clip's boundary-hugging guard treated a TRANSVERSAL tooth-gap
//!   interval near a flank as a boundary-running trace and fell back to
//!   the crossing hull, which bridged two teeth through the gap;
//! - the far-tier endpoint re-trim moved an endpoint that was ALREADY
//!   an exact junction onto the next junction along its curve;
//! - the lateral slop band snapped an exact plane × plane line's end
//!   onto a neighbouring chord vertex 1e-4 away.

use popkern_part::prelude::*;

fn inv(a: f64) -> f64 {
    a.tan() - a
}
fn roll(rb: f64, r: f64) -> f64 {
    let x = rb / r;
    if x >= 1.0 { 0.0 } else { x.acos() }
}

#[allow(clippy::too_many_arguments)]
fn outline(
    z: usize,
    rb: f64,
    ra: f64,
    rf: f64,
    cx: f64,
    cy: f64,
    widen: f64,
    phase: f64,
    pressure: f64,
    steps: usize,
    teeth: usize,
) -> Vec<(f64, f64)> {
    let pi = std::f64::consts::PI;
    let inv_p = inv(pressure);
    let half_p = pi / (2.0 * z as f64);
    let pitch = 2.0 * pi / z as f64;
    let r_lo = rb.max(rf);
    let half_lo = half_p + inv_p - inv(roll(rb, r_lo)) + widen;
    let has_valley = r_lo - rf >= 1e-12;
    let mut pts = Vec::new();
    for k in 0..teeth {
        let rot = phase + k as f64 * pitch;
        for s in 0..=steps {
            let r = r_lo + (ra - r_lo) * s as f64 / steps as f64;
            let half = half_p + inv_p - inv(roll(rb, r)) + widen;
            pts.push((cx + r * (rot - half).cos(), cy + r * (rot - half).sin()));
        }
        for s in 0..=steps {
            let r = ra - (ra - r_lo) * s as f64 / steps as f64;
            let half = half_p + inv_p - inv(roll(rb, r)) + widen;
            pts.push((cx + r * (rot + half).cos(), cy + r * (rot + half).sin()));
        }
        if has_valley {
            let vf = rot + half_lo;
            let vt = rot + pitch - half_lo;
            for s in 0..=3 {
                let th = vf + (vt - vf) * s as f64 / 3.0;
                pts.push((cx + rf * th.cos(), cy + rf * th.sin()));
            }
        }
    }
    pts
}

fn gears(steps: usize, misphase_deg: f64) -> (Part, Part) {
    let pi = std::f64::consts::PI;
    let pressure = 30f64.to_radians();
    let z = 9usize;
    let module = 1.5 / 27.0;
    let rp = module * z as f64 / 2.0;
    let rb = rp * pressure.cos();
    let ra = rp + module;
    let rf = rp - 1.25 * module;
    let sun = outline(z, rb, ra, rf, 0.0, 0.0, -0.004, 0.0, pressure, steps, z);
    let psi = pi + pi / z as f64 + misphase_deg.to_radians();
    let planet = outline(
        z,
        rb,
        ra,
        rf,
        2.0 * rp,
        0.0,
        -0.004,
        psi,
        pressure,
        steps,
        z,
    );
    (
        Profile::on(Plane::XY)
            .polygon(&sun)
            .extrude(0.4)
            .expect("sun"),
        Profile::on(Plane::XY)
            .polygon(&planet)
            .extrude(0.4)
            .expect("planet"),
    )
}

fn mesh_volume(p: &Part) -> f64 {
    let m = p
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    m.triangles
        .iter()
        .map(|t| {
            let (a, b, c) = (
                m.vertices[t[0]].position,
                m.vertices[t[1]].position,
                m.vertices[t[2]].position,
            );
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum()
}

/// (chords per flank, misphase °, OpenCASCADE faces, OpenCASCADE volume).
const CASES: &[(usize, f64, usize, f64)] = &[
    (1, 5.0, 5, 3.33714061205813e-05),
    (1, 10.0, 8, 0.000627113164049551),
    (1, 15.0, 8, 0.0012797851487246282),
    (1, 20.0, 8, 0.00148289550577439),
    (2, 2.0, 5, 1.2223118293839281e-05),
    (2, 5.0, 7, 0.0002878655285171347),
    (2, 10.0, 10, 0.0010662674609961064),
    (2, 20.0, 12, 0.0020560016112605474),
    (6, 2.0, 9, 7.054803860103956e-05),
    (6, 5.0, 14, 0.0003930647136888164),
    (6, 10.0, 18, 0.001206948028766886),
    (6, 15.0, 20, 0.0019552397148431483),
    (6, 20.0, 20, 0.0021937927705436133),
];

#[test]
fn coplanar_capped_toothed_prisms_intersect_like_opencascade() {
    let mut failures = Vec::new();
    for &(steps, misphase, faces, volume) in CASES {
        let (sun, planet) = gears(steps, misphase);
        match sun.intersect(&planet) {
            Ok(common) => {
                let v = mesh_volume(&common);
                let open = common.open_edge_count().expect("mesh");
                if common.face_count() != faces || open != 0 || ((v - volume) / volume).abs() > 1e-9
                {
                    failures.push(format!(
                        "{steps} chords, {misphase}°: faces {} (OCC {faces}), open {open}, volume {v:e} (OCC {volume:e})",
                        common.face_count()
                    ));
                }
            }
            Err(e) => failures.push(format!("{steps} chords, {misphase}°: declined: {e}")),
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
