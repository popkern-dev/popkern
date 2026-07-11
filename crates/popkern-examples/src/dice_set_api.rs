//! Polyhedral dice set, written against the `popkern-part` modeling
//! API — the acceptance gate for the API's phase 2 (see
//! `docs/architecture/part-modeling-api-design.md`).
//!
//! Same dice, same gates as `dice_set`: d4/d6/d8/d10/d12/d20 as
//! exact polyhedra with truncated corners and engraved font numbers
//! (the d4 reading per-corner, as real d4s do). Every die must
//! validate, tessellate watertight, and export STEP.

use popkern_core::geometry::Point3;
use popkern_part::poly::PolyFace;
use popkern_part::prelude::*;
use popkern_part::text::Font;
use popkern_part::{EngravePlacement, engrave_polygons, face_frame, polygons_bbox};

/// Canonical unit-scale vertex coordinates for one die.
fn die_vertices(sides: usize) -> Vec<Point3> {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut v = Vec::new();
    match sides {
        4 => {
            for p in [(1., 1., 1.), (1., -1., -1.), (-1., 1., -1.), (-1., -1., 1.)] {
                v.push(Point3::new(p.0, p.1, p.2));
            }
        }
        6 => {
            for x in [-1.0, 1.0] {
                for y in [-1.0, 1.0] {
                    for z in [-1.0, 1.0] {
                        v.push(Point3::new(x, y, z));
                    }
                }
            }
        }
        8 => {
            for s in [-1.0, 1.0] {
                v.push(Point3::new(s, 0.0, 0.0));
                v.push(Point3::new(0.0, s, 0.0));
                v.push(Point3::new(0.0, 0.0, s));
            }
        }
        10 => {
            let h = 0.12;
            for k in 0..10 {
                let th = std::f64::consts::TAU * k as f64 / 10.0;
                let z = if k % 2 == 0 { h } else { -h };
                v.push(Point3::new(th.cos(), th.sin(), z));
            }
            let (a, b, c) = (v[0], v[1], v[2]);
            let n = (b - a).cross(&(c - a));
            let z0 = (n.dot(&(a - Point3::new(0.0, 0.0, 0.0)))) / n.z;
            v.push(Point3::new(0.0, 0.0, z0.abs()));
            v.push(Point3::new(0.0, 0.0, -z0.abs()));
        }
        12 => {
            for x in [-1.0, 1.0] {
                for y in [-1.0, 1.0] {
                    for z in [-1.0, 1.0] {
                        v.push(Point3::new(x, y, z));
                    }
                }
            }
            let inv = 1.0 / phi;
            for s1 in [-1.0, 1.0] {
                for s2 in [-1.0, 1.0] {
                    v.push(Point3::new(0.0, s1 * inv, s2 * phi));
                    v.push(Point3::new(s1 * inv, s2 * phi, 0.0));
                    v.push(Point3::new(s1 * phi, 0.0, s2 * inv));
                }
            }
        }
        20 => {
            for s1 in [-1.0, 1.0] {
                for s2 in [-1.0, 1.0] {
                    v.push(Point3::new(0.0, s1, s2 * phi));
                    v.push(Point3::new(s1, s2 * phi, 0.0));
                    v.push(Point3::new(s2 * phi, 0.0, s1));
                }
            }
        }
        other => panic!("no die with {other} sides"),
    }
    v
}

/// Engrave one number per corner of a d4 face, each reading toward
/// its corner (the rolled vertex shows its number upright at the top
/// of all three visible faces).
fn engrave_d4_corners(
    points: &mut Vec<Point3>,
    faces: &mut Vec<PolyFace>,
    face_idx: usize,
    corners: &[(Point3, usize)],
    font: &Font,
    depth: f64,
) -> Result<(), PartError> {
    let frame = face_frame(points, &faces[face_idx]);
    let outline: Vec<Point3> = faces[face_idx].outer.iter().map(|&i| points[i]).collect();
    for &(corner, number) in corners {
        let toward = corner - frame.centroid;
        let toward = toward - frame.normal * toward.dot(&frame.normal);
        let up = toward
            .normalized()
            .ok_or_else(|| PartError::NotApplicable("degenerate corner direction".into()))?;
        let right = up.cross(&frame.normal);
        // Reach: centroid → truncated outline along the corner ray.
        let mut dist = f64::INFINITY;
        let m = outline.len();
        for i in 0..m {
            let (a, b) = (outline[i], outline[(i + 1) % m]);
            let ab = b - a;
            let denom = up.cross(&ab).dot(&frame.normal);
            if denom.abs() < 1e-12 {
                continue;
            }
            let ac = a - frame.centroid;
            let t = ac.cross(&ab).dot(&frame.normal) / denom;
            let u = ac.cross(&up).dot(&frame.normal) / denom;
            if t > 0.0 && (-1e-9..=1.0 + 1e-9).contains(&u) {
                dist = dist.min(t);
            }
        }
        let polys = font.polygons(&number.to_string())?;
        let (_, _, min_y, max_y) = polygons_bbox(&polys);
        let placement = EngravePlacement {
            center: frame.centroid + up * (0.55 * dist),
            right,
            up,
            scale: 0.30 * dist / (max_y - min_y),
        };
        engrave_polygons(points, faces, face_idx, &polys, &placement, depth)?;
    }
    Ok(())
}

fn build_die(sides: usize, font: &Font) -> Result<Part, PartError> {
    let radius = 10.0;
    let trim = 0.18;
    let depth = 0.5;
    let raw = die_vertices(sides);
    let max_r = raw
        .iter()
        .map(|p| (*p - Point3::new(0.0, 0.0, 0.0)).length())
        .fold(0.0_f64, f64::max);
    let scaled: Vec<Point3> = raw
        .iter()
        .map(|p| {
            Point3::new(
                p.x / max_r * radius,
                p.y / max_r * radius,
                p.z / max_r * radius,
            )
        })
        .collect();

    let solid = Part::polyhedron(&scaled, Tolerance::default())?;
    let trimmed = solid.truncate_vertices(trim)?;

    // Numbered faces are the truncation's 2n-gons — one per original
    // hull face, identified by having twice that face's vertex count
    // ... simply: the first `sides` faces extracted correspond to
    // the original faces (truncation emits them first, and
    // extraction preserves construction order).
    if sides == 4 {
        // Per-vertex numbers at each corner (apex-read convention).
        let hull_faces = popkern_part::poly::convex_hull_faces(&scaled);
        trimmed.edit_polyhedron("engrave d4 corner numbers".to_string(), |points, faces| {
            for (f, hull_face) in hull_faces.iter().enumerate() {
                let corners: Vec<(Point3, usize)> =
                    hull_face.iter().map(|&v| (scaled[v], v + 1)).collect();
                engrave_d4_corners(points, faces, f, &corners, font, depth)?;
            }
            Ok(())
        })
    } else {
        trimmed.edit_polyhedron(format!("engrave numbers 1..{sides}"), |points, faces| {
            for f in 0..sides {
                let number = if sides == 10 { f } else { f + 1 };
                let polys = font.polygons(&number.to_string())?;
                let frame = face_frame(points, &faces[f]);
                let placement = EngravePlacement::centered(&frame, &polys, 0.62);
                engrave_polygons(points, faces, f, &polys, &placement, depth)?;
            }
            Ok(())
        })
    }
}

fn main() {
    let font = Font::system_sans().expect("system font");
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("out dir");

    for sides in [4usize, 6, 8, 10, 12, 20] {
        let die = build_die(sides, &font)
            .unwrap_or_else(|e| panic!("d{sides} build failed: {e}"))
            .assert_watertight()
            .unwrap_or_else(|e| panic!("d{sides} watertight gate: {e}"));
        let path = out_dir.join(format!("die_d{sides}_api.step"));
        die.export_step(&path)
            .unwrap_or_else(|e| panic!("d{sides} export: {e}"));
        println!(
            "d{sides}: {} faces ✓ watertight, wrote {}",
            die.face_count(),
            path.display()
        );
    }
}
