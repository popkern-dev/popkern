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
            // The proven bound covers the sampled witness; an excess is
            // the mesh's (its witness is past the tolerance too).
            assert!(f.surface_bound >= f.witness, "{name}: {f:?}");
            if f.status == ChordalStatus::Exceeds && f.surface_bound > options.chord_height {
                assert!(
                    f.witness > options.chord_height,
                    "{name}: a loose bound, not a miss: {f:?}"
                );
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
    // Measured 2026-10-20 (PK-505): 271 faces, 243 certified, 13 proven
    // past the chord height (PK-903: cones, tori, Gregory corner caps),
    // 15 inconclusive (PK-904: Gregory corners). A ratchet: the counts
    // may only improve.
    assert!(
        certified >= 243 && exceeds <= 13 && inconclusive <= 15,
        "{certified} certified, {exceeds} exceed, {inconclusive} inconclusive"
    );
    eprintln!(
        "[pk505] corpus: {faces} faces, {certified} certified, {exceeds} exceed, {inconclusive} inconclusive; mesh {mesh_time:.2}s, certificate {cert_time:.2}s"
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
