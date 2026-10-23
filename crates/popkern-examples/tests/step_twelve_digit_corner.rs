//! PK-915: why OpenCASCADE's rewrite of the notched saddle's STEP export
//! drops two corner-layer seams from certified G2 to G1 that popkern's
//! own reader keeps at G2.
//!
//! Attribution (2026-10-23): OCC writes reals to 12 significant digits
//! (popkern 16). Rounding popkern's own export to 12 digits reproduces
//! OCC's readings on every seam (the two drop seams: 5.61e-2 and
//! 6.69e-2 against OCC's 5.67e-2 and 6.67e-2, gate 5e-2), so the
//! rewrite's pcurve forms do not matter. The drop is geometric: on the
//! 12-digit file the certified LOWER bound of the curvature jump exceeds
//! or meets the gate within 1.2e-4 of the seams' cap-corner end (4.96e-2
//! and 5.14e-2), where the export's conditioned corner block amplifies a
//! 1e-12 coordinate change into a curvature change. Natively the lower
//! bound there is 0. A 200-sample continuity measure misses that first
//! 1e-4 of the seam and reads 1e-6. Which side of the gate the rounded
//! corner lands on is phase luck in the corner control points, which is
//! why a 1e-3 relative cap change (PK-600) flipped it. These seams join
//! the interchange class of docs/decisions/
//! 2026-09-26-interchange-continuity-at-setback-corners.md.
//!
//! This test keeps the attribution without the OCC venv: the 12-digit
//! leg must reproduce the OCC leg's saddle drops, and the two extra
//! drops must be certified geometric (lower bound over 0.9x the gate).

use popkern_blend::ContinuityGates;
use popkern_blend::blend_continuity_certification::evidence_certify_options;
use popkern_core::geometry::Point3;
use popkern_core::geometry::certify::SeamCertifyOptions;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::edge_seam_certificate::certify_edge_seam;
use popkern_core::topology::{BRep, EdgeHandle};

#[path = "step_roundtrip_continuity_support/kp_caps.rs"]
#[allow(dead_code)]
mod kp_caps;

/// Every real in the DATA section rewritten to `digits` significant
/// digits, as a reader that writes fewer digits would.
fn round_reals(step: &str, digits: usize) -> String {
    let Some(start) = step.find("DATA;") else {
        return step.to_string();
    };
    let (head, data) = step.split_at(start);
    let bytes = data.as_bytes();
    let mut out = String::with_capacity(data.len());
    let mut i = 0;
    let numeric =
        |c: u8| c.is_ascii_digit() || c == b'.' || c == b'E' || c == b'e' || c == b'-' || c == b'+';
    while i < bytes.len() {
        let c = bytes[i];
        let prev = if i == 0 { b' ' } else { bytes[i - 1] };
        let starts = (c.is_ascii_digit() || c == b'-' || c == b'.')
            && !(prev.is_ascii_alphanumeric() || prev == b'#' || prev == b'_' || prev == b'.');
        if starts {
            let mut j = i + 1;
            while j < bytes.len() && numeric(bytes[j]) {
                j += 1;
            }
            let token = &data[i..j];
            if token.contains('.')
                && let Ok(x) = token.parse::<f64>()
            {
                out.push_str(&format!("{:.*E}", digits - 1, x));
            } else {
                out.push_str(token);
            }
            i = j;
        } else {
            out.push(c as char);
            i += 1;
        }
    }
    format!("{head}{out}")
}

fn export(brep: &BRep) -> String {
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: "saddle".into(),
            id: "SADDLE".into(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("step_twelve_digit_corner");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("saddle.step");
    popkern_io::export_step_file(&doc, &path, &popkern_io::StepExportOptions::default())
        .expect("export");
    std::fs::read_to_string(&path).expect("read back")
}

fn import(text: &str) -> Vec<BRep> {
    let result =
        popkern_io::import_step(text.as_bytes(), &popkern_io::StepImportOptions::default())
            .expect("import");
    result
        .document
        .products
        .into_iter()
        .flat_map(|p| p.shapes)
        .flat_map(|s| s.solids.clone())
        .collect()
}

/// The seam of `bodies` with `edge`'s endpoints and midpoint.
fn matched(
    brep: &BRep,
    edge: EdgeHandle,
    bodies: &[BRep],
    tol: f64,
) -> Option<(usize, EdgeHandle)> {
    let ends = |b: &BRep, e: EdgeHandle| -> Option<([Point3; 2], Point3)> {
        let d = b.edges.get(e)?;
        let c = d.curve.as_ref()?;
        let h = b.half_edges.get(d.half_edge)?;
        let t = b.half_edges.get(h.twin)?;
        Some((
            [
                b.vertices.get(h.origin)?.position,
                b.vertices.get(t.origin)?.position,
            ],
            c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
        ))
    };
    let (key, mid) = ends(brep, edge)?;
    let near = |a: &Point3, b: &Point3| a.distance_to(b) <= tol;
    bodies.iter().enumerate().find_map(|(i, b)| {
        b.edges
            .keys()
            .find(|&e| {
                ends(b, e).is_some_and(|(k, m)| {
                    ((near(&key[0], &k[0]) && near(&key[1], &k[1]))
                        || (near(&key[0], &k[1]) && near(&key[1], &k[0])))
                        && near(&mid, &m)
                })
            })
            .map(|e| (i, e))
    })
}

#[test]
fn a_twelve_digit_export_crosses_the_saddle_corner_seams_at_g1() {
    let gates = ContinuityGates::commercial();
    let options = |b: &BRep| SeamCertifyOptions {
        budget: 1024,
        ..evidence_certify_options(b, &gates)
    };
    let body = kp_caps::notched_saddle().expect("saddle");
    let text = export(&body);
    let full = import(&text);
    let twelve = import(&round_reals(&text, 12));
    let g2 = |b: &BRep, e: EdgeHandle| certify_edge_seam(b, e, &options(b));
    let (mut both, mut twelve_only) = (0, Vec::new());
    for edge in body.edges.keys() {
        let native = g2(&body, edge);
        if !(native.bounds_cover_seam() && native.g2_upper <= gates.g2_curvature) {
            continue;
        }
        let read = |bodies: &[BRep]| {
            matched(&body, edge, bodies, gates.g0_gap).map(|(i, e)| g2(&bodies[i], e))
        };
        let (Some(f), Some(t)) = (read(&full), read(&twelve)) else {
            panic!("native G2 seam {edge:?} unmatched after import");
        };
        let drops =
            |c: &popkern_core::geometry::certify::SeamCertificate| c.g2_upper > gates.g2_curvature;
        match (drops(&f), drops(&t)) {
            (true, true) => both += 1,
            (false, true) => twelve_only.push(t),
            (true, false) => panic!("{edge:?}: drops at 16 digits but not at 12"),
            (false, false) => {}
        }
    }
    // The popkern leg's three corner-layer seams (PK-575), and the two
    // the 12-digit file adds, as OCC's leg measures.
    assert_eq!(both, 3, "seams dropping at both precisions");
    assert_eq!(twelve_only.len(), 2, "seams dropping at 12 digits only");
    for c in &twelve_only {
        // Geometric, not certifier slack: a certified true value of the
        // jump sits at the gate, at the seam's cap-corner end.
        assert!(
            c.lower[2] > 0.9 * gates.g2_curvature,
            "lower bound {:.3e} at {:?}",
            c.lower[2],
            c.worst_at[2]
        );
        let (a, b) = c.worst_at[2];
        assert!(
            a.min(1.0 - b) < 1e-3,
            "worst interval {:?} away from both ends",
            (a, b)
        );
    }
}
