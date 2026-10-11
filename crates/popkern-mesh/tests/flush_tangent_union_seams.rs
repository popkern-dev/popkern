//! PK-750: a closed-rim rod (`BRep::cylinder`: one wall face, one seam)
//! unioned with a box flush with its caps and tangent to it along the
//! box's own vertical edges, the rod's seam turned to every class of
//! placement: on either tangent line (the box's cut line takes the seam
//! over), on the far side, inside the box (the seam vertex lies in the
//! coplanar overlap), and at a generic azimuth. Every placement serves
//! the same half stadium on its charts — no graph-merge fallback, no
//! face walked off-chart — watertight, with its volume.
//!
//! Each placement has failed differently (PK-777 added the last two): rays crossing the cap plane
//! just outside the rim counted as cap hits (a below-rim probe read
//! Inside the rod, the cap's outward normal calibrated inverted, and
//! both copies of the coplanar overlap were dropped); the redundant-
//! edge pass compared orientation FLAGS, not oriented normals, and left
//! an anti-parallel-normal cap split; a cut line lying on the seam left
//! the chart without a seam; with the seam vertex inside the overlap
//! the cap's rim piece and the box's overlap record on the same arc
//! became two edges (the assembler compared PARAMETER midpoints, and
//! an angle arc and a rational quarter-circle span disagree off the
//! quarter grid) and the union declined non-manifold at every seam
//! angle but 0°; and a closed rim's closest-point clamp read a record
//! near the seam off the wall's boundary, so the wall's chart declined
//! at −60° and −75°. With the seam vertex inside the overlap the shared
//! cap arc is a two-edge chain, and the chain merge still compared
//! orientation flags: the caps stayed split (8 faces) at every such
//! placement (PK-778).

use popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events;
use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const R: f64 = 4.5;
const H: f64 = 16.0;

fn seamed_rod(seam_deg: f64) -> BRep {
    let rod = BRep::cylinder(R, H, Tolerance::default()).expect("rod");
    let t = Transform3::translation(-6.0, 0.0, 0.0)
        .compose(&Transform3::rotation_z(seam_deg.to_radians()));
    let seam = t.apply_point(&Point3::new(R, 0.0, 0.0));
    let (s, c) = seam_deg.to_radians().sin_cos();
    let want = Point3::new(-6.0 + R * c, R * s, 0.0);
    assert!(seam.distance_to(&want) < 1e-12, "seam at {seam:?}");
    rod.transform(&t).expect("placed")
}

/// The 12 × 2R box between x = ±6, flush with the rod's caps.
fn slab() -> BRep {
    let profile = [
        Point3::new(-6.0, -R, 0.0),
        Point3::new(6.0, -R, 0.0),
        Point3::new(6.0, R, 0.0),
        Point3::new(-6.0, R, 0.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), H, Tolerance::default()).expect("slab")
}

fn mesh_volume_and_open_edges(body: &BRep) -> (f64, usize) {
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    let volume = mesh
        .triangles
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
        .sum();
    // Welded edge multiplicity: every undirected edge used exactly twice.
    let key = |i: usize| {
        let p = mesh.vertices[i].position;
        (
            (p.x * 1e6).round() as i64,
            (p.y * 1e6).round() as i64,
            (p.z * 1e6).round() as i64,
        )
    };
    let mut uses: std::collections::HashMap<_, usize> = std::collections::HashMap::new();
    for t in &mesh.triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (ka, kb) = (key(a), key(b));
            let e = if ka < kb { (ka, kb) } else { (kb, ka) };
            *uses.entry(e).or_default() += 1;
        }
    }
    (volume, uses.values().filter(|&&n| n != 2).count())
}

#[test]
fn a_seamed_rod_unions_with_the_flush_tangent_box_wherever_its_seam_lies() {
    let want = (12.0 * 2.0 * R + 0.5 * std::f64::consts::PI * R * R) * H;
    let mut failures = Vec::new();
    // Every 15° of seam azimuth: on either tangent line (±90°), inside
    // the coplanar overlap, and on the far side.
    for k in -6..=18 {
        let seam_deg = 15.0 * k as f64;
        let _ = take_merge_fallback_events();
        let _ = take_chart_walk_fallback_events();
        let served = boolean(
            BoolOperation::Union,
            seamed_rod(seam_deg),
            slab(),
            Tolerance::default(),
            None,
        );
        let body = match served {
            Ok(b) => b,
            Err(e) => {
                failures.push(format!("{seam_deg}°: {e}"));
                continue;
            }
        };
        if let Err(e) = body.validate() {
            failures.push(format!("{seam_deg}°: invalid {e}"));
            continue;
        }
        let merges = take_merge_fallback_events();
        let walks = take_chart_walk_fallback_events();
        if !merges.is_empty() || !walks.is_empty() {
            failures.push(format!("{seam_deg}°: fallbacks {merges:?} {walks:?}"));
        }
        // Top, bottom, two flat sides, the far box side, the half wall,
        // wherever the seam lies: a seam vertex inside the overlap splits
        // the shared cap arc into a two-edge chain, which merges like a
        // single edge (PK-778).
        let faces = 6;
        if body.faces.len() != faces {
            failures.push(format!(
                "{seam_deg}°: {} faces, want {faces}",
                body.faces.len()
            ));
        }
        let (volume, open) = mesh_volume_and_open_edges(&body);
        // The mesh chords the curved wall: ~1e-4 relative.
        let rel = (volume - want).abs() / want;
        if open != 0 || rel > 1e-3 {
            failures.push(format!(
                "{seam_deg}°: {open} open mesh edges, volume {volume} vs {want} (rel {rel:.1e})"
            ));
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}
