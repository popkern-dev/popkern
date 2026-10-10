//! PK-752: a disc (axis y) unioned with a boss whose axis runs through
//! the disc's axis, the boss entering from every azimuth so that its
//! intersection curve crosses the disc's seam at every relative angle,
//! including the curve's apex exactly on the seam. The scenes are
//! congruent (the boss turns about the disc's axis), so every azimuth
//! serves on its charts — no graph-merge fallback, no face walked
//! off-chart — watertight, with one volume.

use popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events;
use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const RF: f64 = 50.0;
const RB: f64 = 42.5;
const YB: f64 = -45.0;

/// The disc: radius RF about the y axis, y from -13 to -1, its seam at
/// azimuth `seam_deg` about the y axis (measured from +x toward +z).
fn disc(seam_deg: f64) -> BRep {
    // BRep::cylinder runs along +z with its seam at +x: turn z to -y
    // (rotation about x by +90°), then spin the seam about y.
    let rod = BRep::cylinder(RF, 12.0, Tolerance::default()).expect("disc");
    let t = Transform3::translation(0.0, -1.0, 0.0)
        .compose(&Transform3::rotation_y(-seam_deg.to_radians()))
        .compose(&Transform3::rotation_x(std::f64::consts::FRAC_PI_2));
    rod.transform(&t).expect("placed")
}

/// The boss: radius RB, axis through (0, YB, 0), entering from azimuth
/// `deg` and ending on the disc's axis.
fn boss(deg: f64) -> BRep {
    let rod = BRep::cylinder(RB, 90.0, Tolerance::default()).expect("boss");
    // Along +z from z = 0..90; turn +z to point from azimuth `deg`
    // toward the axis, then lift to y = YB.
    let t = Transform3::translation(0.0, YB, 0.0)
        .compose(&Transform3::rotation_y((-deg - 90.0).to_radians()))
        .compose(&Transform3::translation(0.0, 0.0, -90.0));
    rod.transform(&t).expect("placed")
}

fn mesh_volume(body: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.vertices[i].position);
            let (a, b, c) = (
                Vec3::new(a.x, a.y, a.z),
                Vec3::new(b.x, b.y, b.z),
                Vec3::new(c.x, c.y, c.z),
            );
            a.dot(&b.cross(&c)) / 6.0
        })
        .sum()
}

#[test]
fn the_boss_unions_from_every_azimuth_across_the_disc_seam() {
    // The boss enters from azimuth 180° (the -x side) in the repro; the
    // disc's seam sits at 180° too, then sweeps a full turn relative to it.
    let probe = boss(180.0);
    let entry =
        probe
            .vertices
            .values()
            .map(|v| v.position)
            .fold(
                Point3::new(0.0, 0.0, 0.0),
                |a, p| if p.x < a.x { p } else { a },
            );
    assert!(entry.x < -80.0, "boss enters from -x: {entry:?}");
    let mut volumes = Vec::new();
    let mut failures = Vec::new();
    // Every 30° from the curve's apex (180°: the apex exactly on the
    // seam) through a full turn.
    for k in 0..12 {
        let seam = 180.0 + 30.0 * k as f64;
        let _ = take_merge_fallback_events();
        let _ = take_chart_walk_fallback_events();
        let body = match boolean(
            BoolOperation::Union,
            disc(seam),
            boss(180.0),
            Tolerance::default(),
            None,
        ) {
            Ok(b) => b,
            Err(e) => {
                failures.push(format!("seam {seam}°: {e}"));
                continue;
            }
        };
        if let Err(e) = body.validate() {
            failures.push(format!("seam {seam}°: invalid {e}"));
            continue;
        }
        let merges = take_merge_fallback_events();
        let walks = take_chart_walk_fallback_events();
        if !merges.is_empty() || !walks.is_empty() {
            failures.push(format!("seam {seam}°: fallbacks {merges:?} {walks:?}"));
        }
        // Tessellation dominates the cost: volumes every 90° (the apex
        // on the seam among them) carry the congruence check.
        if k % 3 == 0 {
            volumes.push((seam, mesh_volume(&body)));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
    let (lo, hi) = volumes
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), (_, v)| {
            (lo.min(*v), hi.max(*v))
        });
    // Congruent scenes: one volume to the mesher's chord tolerance.
    assert!((hi - lo) / hi < 1e-3, "volumes {volumes:?}");
}
