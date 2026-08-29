//! PK-244 consumer-level gates for orientation-flag honesty.
//!
//! The per-face audit (`BRep::audit_winding_honesty`) measures the
//! invariant directly; these tests verify it through the consumers
//! that break when a flag lies: mesher triangle winding (a dishonest
//! face's triangles come out mirrored, so the divergence-theorem sum
//! over raw mesh winding diverges from the flag-folded quadrature
//! volume) and the corpus serving boundary.
//!
//! Fixture: the pyramid4 apex fillet — its N=4 partition Gregory
//! caps shipped flag-dishonest (loop CW in UV with `orientation:
//! true`) before the serving-boundary reconcile.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn pyramid4_apex_fillet() -> BRep {
    let brep = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let apex = brep
        .vertices
        .keys()
        .find(|&v| {
            brep.half_edges
                .values()
                .filter(|he| he.origin == v)
                .map(|he| he.edge)
                .collect::<std::collections::BTreeSet<_>>()
                .len()
                == 4
        })
        .unwrap();
    let edges: Vec<_> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.15),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: true,
        use_biseptic_gregory: true,
        use_two_stripe_gregory: true,
        use_kp_binonic_gregory: true,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &opts, Tolerance::default(), None)
        .unwrap()
        .brep
}

/// The audit itself must be clean at the fillet serving boundary.
#[test]
fn pyramid4_apex_serves_winding_honest_flags() {
    let out = pyramid4_apex_fillet();
    let dishonest = out.dishonest_winding_faces();
    assert!(
        dishonest.is_empty(),
        "pyramid4 apex fillet serves {} flag-dishonest face(s): {dishonest:?}",
        dishonest.len()
    );
}

/// Mesher consumer: the divergence-theorem volume summed over RAW
/// mesh triangle winding must agree with the flag-folded measure
/// volume. A flag-dishonest face mirrors its triangles, so the two
/// diverge by twice that face's (position-dependent) contribution —
/// pre-fix the pyramid4 caps were exactly such faces.
#[test]
fn pyramid4_mesh_winding_agrees_with_measured_volume() {
    let out = pyramid4_apex_fillet();
    let tol = Tolerance::default();
    let measured = popkern_measure::body_volume(&out, &tol).expect("volume measurable");

    let mesh = popkern_mesh::tessellate(&out, &popkern_mesh::TessellationOptions::default())
        .expect("tessellates");
    let mesh_signed: f64 = mesh
        .triangles
        .iter()
        .map(|tri| {
            let p0 = mesh.vertices[tri[0]].position;
            let p1 = mesh.vertices[tri[1]].position;
            let p2 = mesh.vertices[tri[2]].position;
            let v0 = popkern_core::geometry::Vec3::new(p0.x, p0.y, p0.z);
            let v1 = popkern_core::geometry::Vec3::new(p1.x, p1.y, p1.z);
            let v2 = popkern_core::geometry::Vec3::new(p2.x, p2.y, p2.z);
            v0.dot(&v1.cross(&v2)) / 6.0
        })
        .sum();

    let rel = (mesh_signed.abs() - measured).abs() / measured.max(1.0e-12);
    assert!(
        rel < 0.02,
        "mesh-winding volume {:.6} vs measured volume {:.6} (rel {:.3}): \
         a divergence here means some face's triangles are wound against \
         the shell (flag-dishonest at the mesher)",
        mesh_signed.abs(),
        measured,
        rel
    );
}

#[path = "fillet_mileage_corpus_support/parts.rs"]
mod parts;

/// Serving-boundary audit across the mileage corpus: every part that
/// builds must serve zero flag-dishonest faces. This reaches the
/// boolean/sweep/draft/fillet author mix of real part recipes, not
/// just the fillet fixture classes.
#[test]
fn mileage_corpus_serves_winding_honest_flags() {
    use popkern_part::prelude::*;
    type PartRecipe = (&'static str, fn() -> Result<Part, PartError>);
    let corpus: Vec<PartRecipe> = vec![
        ("mounting_bracket", parts::build_mounting_bracket),
        ("ribbed_plate", parts::build_ribbed_plate),
        ("bossed_housing", parts::build_bossed_housing),
        ("turned_bushing", parts::build_turned_bushing),
        ("pocketed_tray", parts::build_pocketed_tray),
        ("slotted_block", parts::build_slotted_block),
        ("countersunk_plate", parts::build_countersunk_plate),
        (
            "oblique_countersunk_plate",
            parts::build_oblique_countersunk_plate,
        ),
        ("drafted_cup", parts::build_drafted_cup),
        ("cross_drilled_boss", parts::build_cross_drilled_boss),
        ("coaxial_bushing", parts::build_coaxial_bushing),
        ("tilted_hole_plate", parts::build_tilted_hole_plate),
        ("offset_bracket", parts::build_offset_bracket),
        ("rim_chamfered_disc", parts::build_rim_chamfered_disc),
        ("star_socket_plate", parts::build_star_socket_plate),
        ("mixed_corner_block", parts::build_mixed_corner_block),
        ("chamfered_peg", parts::build_chamfered_peg),
        ("skewed_dress_bar", parts::build_skewed_dress_bar),
        ("counterbored_plate", parts::build_counterbored_plate),
        ("angled_rib_plate", parts::build_angled_rib_plate),
        ("double_boss_plate", parts::build_double_boss_plate),
        ("notched_flange", parts::build_notched_flange),
        ("pentagon_tower", parts::build_pentagon_tower),
    ];
    let mut offenders: Vec<String> = Vec::new();
    for (name, build) in corpus {
        let Ok(part) = build() else {
            // Build declines are the corpus baseline test's concern.
            continue;
        };
        let brep = part.brep();
        let dishonest: Vec<String> = brep
            .dishonest_winding_faces()
            .into_iter()
            .map(|f| {
                let kind = brep
                    .faces
                    .get(f)
                    .and_then(|d| d.surface.as_ref())
                    .map(|s| match &s.geometry {
                        popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                        popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                        popkern_core::geometry::SurfaceGeometry::Cone { .. } => "Cone",
                        popkern_core::geometry::SurfaceGeometry::Sphere { .. } => "Sphere",
                        popkern_core::geometry::SurfaceGeometry::Torus { .. } => "Torus",
                        popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
                        popkern_core::geometry::SurfaceGeometry::Gregory(_) => "Gregory",
                        _ => "Other",
                    })
                    .unwrap_or("<none>");
                format!("{f:?}[{kind}]")
            })
            .collect();
        eprintln!(
            "[flag-audit] {name}: {} dishonest {dishonest:?}",
            dishonest.len()
        );
        if !dishonest.is_empty() {
            offenders.push(format!("{name}: {}", dishonest.join(", ")));
        }
    }
    assert!(
        offenders.is_empty(),
        "corpus parts serve flag-dishonest faces:\n  {}",
        offenders.join("\n  ")
    );
}
