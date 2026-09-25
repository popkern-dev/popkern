//! PK-509: negative controls for the STEP round-trip audit. A
//! known-good conic-fillet plate is exported, re-imported with healing
//! off, and audited against the native body by four independent gates:
//!
//! - `diagnostics`: our own file must import without a single repair,
//!   decline or deviation-envelope warning;
//! - `validate`: `validate_all` is clean;
//! - `closure`: one solid, with the native face/edge/vertex counts;
//! - `geometry`: native surface samples over each face's trim box lie on
//!   the matched re-imported surface.
//!
//! The positive control passes every gate. One targeted corruption per
//! class of the Part 21 text must then be rejected by at least one
//! gate, and the per-gate rejection counts are asserted, so a gate that
//! goes blind shows as a changed count, not a silent pass. The weight
//! corruption moves only a face interior (every boundary stays
//! consistent); only the geometry gate can see it.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::exact_bounded_parameter_range;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_part::prelude::*;

/// Native samples must lie on the re-imported surface to this
/// distance: the export writes coordinates at full `f64` precision, so
/// only the importer's linear tolerance separates the two.
const GEOMETRY_TOL: f64 = 1.0e-6;
/// Samples per trim-box side for the geometry gate.
const GRID: usize = 9;

fn native_plate() -> Part {
    Profile::on(Plane::XY)
        .rect(60.0, 40.0)
        .extrude(8.0)
        .unwrap_or_else(|e| panic!("plate: {e:?}"))
        .fillet_with(
            |recipe| recipe.radius(3.0).conic(0.35),
            |e| {
                e.straight()
                    .direction(Dir::X)
                    .chord_length_between(55.0, 65.0)
                    .at_height(8.0)
            },
        )
        .unwrap_or_else(|e| panic!("fillet: {e:?}"))
}

fn exported(part: &Part) -> String {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk509");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let path = dir.join("plate.step");
    part.export_step(&path)
        .unwrap_or_else(|e| panic!("export: {e:?}"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{e}"))
}

/// The audit gates that rejected `text` (empty = the audit passes).
fn audit(native: &BRep, text: &str) -> Vec<&'static str> {
    let options = popkern_io::StepImportOptions {
        healing: popkern_io::step::healing::HealingLevel::None,
        ..Default::default()
    };
    let Ok(result) = popkern_io::import_step(text.as_bytes(), &options) else {
        return vec!["import"];
    };
    let mut rejected = Vec::new();
    if !result.diagnostics.is_empty() {
        rejected.push("diagnostics");
    }
    let solids: Vec<&BRep> = result
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .collect();
    let [imported] = solids.as_slice() else {
        rejected.push("closure");
        return rejected;
    };
    if !imported.validate_all().is_empty() {
        rejected.push("validate");
    }
    let counts = |b: &BRep| (b.faces.len(), b.edges.len(), b.vertices.len());
    if counts(imported) != counts(native) {
        rejected.push("closure");
    }
    if geometry_deviation(native, imported) > GEOMETRY_TOL {
        rejected.push("geometry");
    }
    rejected
}

/// The re-imported face with the native face's outer-loop vertex set.
fn matched_face(native: &BRep, face: FaceHandle, imported: &BRep) -> Option<FaceHandle> {
    let corners = native.face_outer_loop_positions(face);
    imported.faces.keys().find(|&candidate| {
        let other = imported.face_outer_loop_positions(candidate);
        other.len() == corners.len()
            && corners
                .iter()
                .all(|p| other.iter().any(|q| p.distance_to(q) <= GEOMETRY_TOL))
    })
}

/// Largest distance from a native surface sample (over each face's
/// exact trim box) to the matched re-imported surface; infinite when a
/// face has no match or a sample cannot be projected.
fn geometry_deviation(native: &BRep, imported: &BRep) -> f64 {
    let mut worst = 0.0_f64;
    for (face, data) in native.faces.iter() {
        let Some(native_surface) = data.surface.as_ref() else {
            continue;
        };
        let Some(surface) =
            matched_face(native, face, imported).and_then(|m| imported.faces[m].surface.as_ref())
        else {
            return f64::INFINITY;
        };
        let corners = native.face_outer_loop_positions(face);
        let ((u0, u1), (v0, v1)) =
            exact_bounded_parameter_range(&corners, &native_surface.geometry);
        for i in 0..=GRID {
            for j in 0..=GRID {
                let (s, t) = (i as f64 / GRID as f64, j as f64 / GRID as f64);
                let p: Point3 = native_surface
                    .geometry
                    .point_at(u0 + s * (u1 - u0), v0 + t * (v1 - v0));
                let Some((u, v)) = surface.geometry.closest_point(&p) else {
                    return f64::INFINITY;
                };
                worst = worst.max(surface.geometry.point_at(u, v).distance_to(&p));
            }
        }
    }
    worst
}

/// The byte range of the first `open ... close` span after `anchor`.
fn span_after(text: &str, anchor: &str, open: &str, close: &str) -> Option<(usize, usize)> {
    let at = text.find(anchor)? + anchor.len();
    let start = at + text[at..].find(open)? + open.len();
    let end = start + text[start..].find(close)?;
    Some((start, end))
}

fn splice(text: &str, (start, end): (usize, usize), body: &str) -> String {
    format!("{}{}{}", &text[..start], body, &text[end..])
}

/// One vertex's `CARTESIAN_POINT` x coordinate shifted by 0.5.
fn corrupt_coordinate(text: &str) -> Option<String> {
    let vp = text.find("=VERTEX_POINT('',#")? + "=VERTEX_POINT('',#".len();
    let id: String = text[vp..]
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect();
    let (start, end) = span_after(text, &format!("#{id}=CARTESIAN_POINT"), "('',(", ",")?;
    let x: f64 = text[start..end].trim().parse().ok()?;
    Some(splice(text, (start, end), &format!("{:?}", x + 0.5)))
}

/// The first rational surface's middle-row interior weight doubled.
/// The row sits at the fillet's mid-length, inside the trimmed face;
/// weights on the boundary rows are untouched, so every edge stays on
/// the corrupted surface and only the face interior moves.
fn corrupt_weight(text: &str) -> Option<String> {
    let (start, end) = span_after(text, "RATIONAL_B_SPLINE_SURFACE(", "((", "))")?;
    let mut rows: Vec<Vec<f64>> = text[start..end]
        .split("),(")
        .map(|row| row.split(',').map(|w| w.trim().parse::<f64>()).collect())
        .collect::<Result<_, _>>()
        .ok()?;
    let middle = rows.len() / 2;
    let row = rows.get_mut(middle)?;
    let inner = row.len() / 2;
    row[inner] *= 2.0;
    let body: Vec<String> = rows
        .iter()
        .map(|r| {
            r.iter()
                .map(|w| format!("{w:?}"))
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect();
    Some(splice(text, (start, end), &body.join("),(")))
}

/// The first surface's u-knot just below the parameter 0.25, where the
/// face trim ends, moved 30% toward its successor. The fillet's
/// u-direction is a ruled straight run, so a knot away from the trim
/// only reparametrizes it; this knot moves where the end edges' pcurves
/// land.
fn corrupt_knot(text: &str) -> Option<String> {
    let entity = text.find("B_SPLINE_SURFACE_WITH_KNOTS((")?;
    let mults_end = entity + text[entity..].find("),(")?;
    let second_end = mults_end + 3 + text[mults_end + 3..].find("),(")?;
    let start = second_end + 3;
    let end = start + text[start..].find(')')?;
    let mut knots: Vec<f64> = text[start..end]
        .split(',')
        .map(|k| k.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .ok()?;
    let k = knots.iter().rposition(|&k| k < 0.25)?;
    knots[k] += 0.3 * (knots.get(k + 1)? - knots[k]);
    let body: Vec<String> = knots.iter().map(|k| format!("{k:?}")).collect();
    Some(splice(text, (start, end), &body.join(",")))
}

/// The first oriented edge's sense flipped.
fn corrupt_orientation(text: &str) -> Option<String> {
    let at = text.find("=ORIENTED_EDGE('',*,*,#")?;
    let end = at + text[at..].find(");")?;
    let entity = &text[at..end];
    let flipped = if entity.ends_with(".T.") {
        entity.replace(".T.", ".F.")
    } else {
        entity.replace(".F.", ".T.")
    };
    Some(splice(text, (at, end), &flipped))
}

/// The last face dropped from the closed shell.
fn corrupt_shell_closure(text: &str) -> Option<String> {
    let (start, end) = span_after(text, "=CLOSED_SHELL", "('',(", "))")?;
    let comma = start + text[start..end].rfind(',')?;
    Some(splice(text, (comma, end), ""))
}

/// A targeted edit of the Part 21 text (`None`: no target found).
type Corruption = fn(&str) -> Option<String>;

#[test]
fn the_step_audit_rejects_one_corruption_per_class() {
    let part = native_plate();
    let native = part.brep();
    let text = exported(&part);
    assert_eq!(
        audit(native, &text),
        Vec::<&str>::new(),
        "positive control: the plate round-trips through every gate"
    );
    let corruptions: [(&str, Corruption); 5] = [
        ("coordinate", corrupt_coordinate),
        ("weight", corrupt_weight),
        ("knot", corrupt_knot),
        ("orientation", corrupt_orientation),
        ("shell_closure", corrupt_shell_closure),
    ];
    let mut verdicts = Vec::new();
    for (class, corrupt) in corruptions {
        let corrupted = corrupt(&text).unwrap_or_else(|| panic!("{class}: no target in the file"));
        assert_ne!(corrupted, text, "{class}: the corruption changed nothing");
        verdicts.push((class, audit(native, &corrupted)));
    }
    // The per-gate rejection counts: the weight corruption is seen only
    // by the geometry gate, and the knot corruption (the surface is the
    // same point set, reparametrized) only by the pcurve fidelity
    // diagnostics.
    assert_eq!(
        verdicts,
        vec![
            ("coordinate", vec!["diagnostics", "validate", "geometry"]),
            ("weight", vec!["geometry"]),
            ("knot", vec!["diagnostics"]),
            ("orientation", vec!["diagnostics", "validate", "geometry"]),
            ("shell_closure", vec!["diagnostics", "closure", "geometry"]),
        ],
        "per-class rejecting gates"
    );
}
