//! PK-509(c): negative controls for the replay determinism check. The
//! workspace whole-body fingerprint (FNV-1a over the serialized body)
//! is the evidence behind every "replay is bit-identical" claim, so it
//! must see the smallest possible geometry change: a one-ulp move of a
//! single stored NURBS control point, weight or vertex coordinate.
//! Positive control: two replays of the reference document agree.

mod support;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// FNV-1a over the serialized body — the workspace-standard
/// whole-body fingerprint (identity pin idiom, PK-276 row 457).
fn body_fingerprint(brep: &BRep) -> u64 {
    let bytes = serde_json::to_vec(brep).expect("output BRep must serialize");
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for b in bytes {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// The next representable `f64` above `x`.
fn one_ulp_up(x: f64) -> f64 {
    f64::from_bits(x.to_bits() + 1)
}

fn replayed_reference_body() -> BRep {
    let mut tree = support::author_reference_document();
    tree.replay(Tolerance::default()).expect("reference replay");
    tree.current_brep().expect("reference body").clone()
}

#[test]
fn the_replay_fingerprint_sees_a_one_ulp_geometry_change() {
    let body = replayed_reference_body();
    let baseline = body_fingerprint(&body);
    assert_eq!(
        baseline,
        body_fingerprint(&replayed_reference_body()),
        "positive control: repeated replay is fingerprint-identical"
    );

    let nurbs_face = body
        .faces
        .keys()
        .find(|&f| {
            matches!(
                body.faces[f].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_))
            )
        })
        .expect("the G2 rounds are NURBS surfaces");
    let mutate_nurbs =
        |mutate: &dyn Fn(&mut popkern_core::geometry::nurbs_surface::NurbsSurface)| {
            let mut corrupted = body.clone();
            if let Some(SurfaceGeometry::Nurbs(surface)) = corrupted.faces[nurbs_face]
                .surface
                .as_mut()
                .map(|s| &mut s.geometry)
            {
                mutate(surface);
            }
            body_fingerprint(&corrupted)
        };
    let control_point = mutate_nurbs(&|s| {
        let mid = s.control_points.len() / 2;
        let p = &mut s.control_points[mid];
        p.x = one_ulp_up(p.x);
    });
    let weight = mutate_nurbs(&|s| {
        let mid = s.weights.len() / 2;
        let w = &mut s.weights[mid];
        *w = one_ulp_up(*w);
    });
    let vertex = {
        let mut corrupted = body.clone();
        let v = corrupted.vertices.keys().next().expect("vertices");
        let p = &mut corrupted.vertices[v].position;
        p.z = one_ulp_up(p.z);
        body_fingerprint(&corrupted)
    };
    for (class, fingerprint) in [
        ("control_point", control_point),
        ("weight", weight),
        ("vertex", vertex),
    ] {
        assert_ne!(
            fingerprint, baseline,
            "a one-ulp {class} change must change the replay fingerprint"
        );
    }
}
