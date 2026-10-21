//! PK-505: the certified chordal bound over the fillet mileage corpus
//! (design `docs/architecture/certified-chordal-bound.md` § Measured):
//! the certificate's cost against tessellation, and which faces the
//! default mesher meets its chord height on, proven.

use popkern_mesh::{ChordalStatus, TessellationOptions, chordal_certificate, tessellate};
use popkern_part::prelude::*;

#[path = "fillet_mileage_corpus_support/parts.rs"]
mod parts;
use parts::*;

type CorpusBuild = fn() -> Result<Part, PartError>;

fn corpus() -> Vec<(&'static str, CorpusBuild)> {
    vec![
        ("mounting_bracket", build_mounting_bracket),
        ("ribbed_plate", build_ribbed_plate),
        ("bossed_housing", build_bossed_housing),
        ("turned_bushing", build_turned_bushing),
        ("pocketed_tray", build_pocketed_tray),
        ("slotted_block", build_slotted_block),
        ("countersunk_plate", build_countersunk_plate),
        ("oblique_countersunk_plate", build_oblique_countersunk_plate),
        ("drafted_cup", build_drafted_cup),
        ("cross_drilled_boss", build_cross_drilled_boss),
        ("coaxial_bushing", build_coaxial_bushing),
        ("tilted_hole_plate", build_tilted_hole_plate),
        ("offset_bracket", build_offset_bracket),
        ("rim_chamfered_disc", build_rim_chamfered_disc),
        ("star_socket_plate", build_star_socket_plate),
        ("mixed_corner_block", build_mixed_corner_block),
        ("chamfered_peg", build_chamfered_peg),
        ("skewed_dress_bar", build_skewed_dress_bar),
        ("counterbored_plate", build_counterbored_plate),
        ("angled_rib_plate", build_angled_rib_plate),
        ("double_boss_plate", build_double_boss_plate),
        ("notched_flange", build_notched_flange),
        ("pentagon_tower", build_pentagon_tower),
    ]
}

#[test]
fn the_corpus_meshes_certify_or_name_their_faces() {
    let options = TessellationOptions::default();
    let (mut mesh_time, mut cert_time) = (0.0, 0.0);
    let (mut faces, mut certified, mut exceeds, mut inconclusive) = (0, 0, 0, 0);
    let mut loose = 0;
    for (name, build) in corpus() {
        let Ok(part) = build() else {
            eprintln!("[pk505] {name}: did not build");
            continue;
        };
        let brep = part.brep();
        let t = std::time::Instant::now();
        let Ok(mesh) = tessellate(brep, &options) else {
            eprintln!("[pk505] {name}: did not mesh");
            continue;
        };
        let meshed = t.elapsed().as_secs_f64();
        let t = std::time::Instant::now();
        let cert = chordal_certificate(brep, &mesh, options.chord_height);
        let proved = t.elapsed().as_secs_f64();
        mesh_time += meshed;
        cert_time += proved;
        for f in &cert.faces {
            // The proven bound covers the sampled witness. An excess is
            // the mesh's when its witness is past the tolerance too, and
            // otherwise the bound's (`loose`).
            assert!(f.surface_bound >= f.witness, "{name}: {f:?}");
            if f.status == ChordalStatus::Exceeds && f.witness <= options.chord_height {
                loose += 1;
            }
        }
        let count =
            |f: &dyn Fn(&ChordalStatus) -> bool| cert.faces.iter().filter(|x| f(&x.status)).count();
        let c = count(&|s| *s == ChordalStatus::Certified);
        let e = count(&|s| *s == ChordalStatus::Exceeds);
        let i = count(&|s| matches!(s, ChordalStatus::Inconclusive(_)));
        faces += cert.faces.len();
        certified += c;
        exceeds += e;
        inconclusive += i;
        eprintln!(
            "[pk505] {name}: {} tris, mesh {meshed:.3}s cert {proved:.3}s, faces {c} ok / {e} exceed / {i} inconclusive, worst {:.3e}",
            mesh.triangles.len(),
            cert.worst()
        );
        for f in cert
            .faces
            .iter()
            .filter(|f| f.status != ChordalStatus::Certified)
        {
            let kind = brep
                .faces
                .get(f.face)
                .and_then(|x| x.surface.as_ref())
                .map(|s| s.geometry.kind_name());
            eprintln!(
                "[pk505]   {kind:?} {:?}: surface {:.3e} (witness {:.3e}) boundary {:.3e}",
                f.status, f.surface_bound, f.witness, f.boundary_bound
            );
        }
    }
    // Measured 2026-10-21 (PK-910): all 271 faces certified. The last
    // five, pentagon_tower Gregory caps beside activated corners (two
    // triangles over the whole patch, 5.4e-3 and 4.3e-2 against 2.2e-3
    // witnesses), certify through refined pairings (4.7e-3..5.0e-3) once
    // box ends within rounding of the corner snap to it (the two at
    // 4.3e-2). PK-909's difference hull and bent corner pieces certified
    // the ten before. A ratchet: the counts may only improve.
    let witnessed = exceeds - loose;
    assert!(
        certified >= 271 && exceeds == 0 && witnessed == 0 && inconclusive == 0,
        "{certified} certified, {exceeds} exceed ({witnessed} witnessed), {inconclusive} inconclusive"
    );
    eprintln!(
        "[pk505] corpus: {faces} faces, {certified} certified, {exceeds} exceed ({loose} unwitnessed), {inconclusive} inconclusive; mesh {mesh_time:.2}s, certificate {cert_time:.2}s"
    );
}

#[test]
fn narrow_fillet_strips_meet_the_chord_height() {
    // The PK-522 witness (2026-09-24): a unit cube with r = 0.05 fillet
    // strips meshed one chord across a quarter arc, 1.27e-2 against a
    // 0.005 chord height. The certificate now proves where it stands.
    let part = Profile::on(Plane::XY)
        .rect(1.0, 1.0)
        .extrude(1.0)
        .and_then(|b| b.fillet(0.05, |e| e.straight()))
        .expect("filleted cube");
    let options = TessellationOptions::default();
    let mesh = tessellate(part.brep(), &options).expect("mesh");
    let cert = chordal_certificate(part.brep(), &mesh, options.chord_height);
    eprintln!(
        "[pk505] narrow strips: {} faces, worst {:.3e}, exceed {}, inconclusive {}",
        cert.faces.len(),
        cert.worst(),
        cert.exceeding().len(),
        cert.inconclusive().len()
    );
    assert!(cert.inconclusive().is_empty(), "{:?}", cert.inconclusive());
    assert!(cert.certified(), "{:?}", cert.exceeding());
}

#[test]
fn a_box_ending_within_rounding_of_an_activated_corner_encloses_like_one_on_it() {
    // PK-910: pentagon_tower's caps mesh a vertex at u = 1 − 1.06e-13
    // beside their activated corner (1, 1). Ending there, a box missed
    // the exact corner the zero tracking needs and its hull spread to
    // 4.2e-2 at any size (against 4.9e-4 on the corner), so the caps'
    // corner pieces never certified. Box ends that near a knot snap to
    // it.
    use popkern_core::geometry::SurfaceEval;
    use popkern_core::geometry::certify::{Iv, IvVec3, SurfaceJetSource};
    let part = build_pentagon_tower().expect("pentagon tower");
    let mut caps = 0;
    for (_, face) in part.brep().faces.iter() {
        let Some(surface) = face.surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        if surface.kind_name() != "Gregory" {
            continue;
        }
        let mut source = SurfaceJetSource::new(surface).expect("jet source");
        let corner = IvVec3::from_point(surface.point_at(1.0, 1.0));
        let mut spread = |hi: f64| {
            let hull = source
                .position_hull(Iv::new(0.97, hi), Iv::new(0.97, 1.0))
                .expect("hull");
            (hull - corner).norm().hi
        };
        let (on, near) = (spread(1.0), spread(1.0 - 1.06e-13));
        assert!(
            near <= on * 1.001 + 1e-12,
            "{near} against {on} on the corner"
        );
        caps += 1;
    }
    assert!(caps >= 5, "{caps} Gregory faces");
}
