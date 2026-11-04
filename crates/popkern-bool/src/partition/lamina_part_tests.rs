//! Sheet ∩ solid and sheet − solid through the partition.
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use super::super::partition;

/// The core's two-face quad over x in [x0, x1], y in [y0, y1] at
/// height `z`, in the one-face sheet form.
fn sheet(x0: f64, x1: f64, y0: f64, y1: f64, z: f64) -> BRep {
    let quad = BRep::sheet_quad(x0, x1, y0, y1, z, Tolerance::default()).expect("quad");
    crate::sheet::lamina::single_face_sheet(&quad).expect("one-face sheet")
}

/// Area of one planar loop (Newell).
fn loop_area(points: &[Point3]) -> f64 {
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for (i, p) in points.iter().enumerate() {
        let q = &points[(i + 1) % points.len()];
        x += (p.y - q.y) * (p.z + q.z);
        y += (p.z - q.z) * (p.x + q.x);
        z += (p.x - q.x) * (p.y + q.y);
    }
    0.5 * (x * x + y * y + z * z).sqrt()
}

/// A planar sheet's area: its outer loops less its holes.
fn area(b: &BRep) -> f64 {
    b.faces
        .values()
        .map(|f| {
            let ring = |l: popkern_core::topology::LoopHandle| {
                let start = b.loops[l].half_edge;
                let mut out = vec![b.vertices[b.half_edges[start].origin].position];
                let mut h = b.half_edges[start].next;
                while h != start {
                    out.push(b.vertices[b.half_edges[h].origin].position);
                    h = b.half_edges[h].next;
                }
                out
            };
            let holes: f64 = f.inner_loops.iter().map(|&l| loop_area(&ring(l))).sum();
            loop_area(&ring(f.outer_loop)) - holes
        })
        .sum()
}

/// (area inside, area outside) of `s` against the unit cube, each
/// result a valid sheet facing the input's way.
fn parts(s: &BRep) -> (f64, f64) {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let p = partition(&[&cube, s], tol).expect("partition");
    let mut out = [0.0; 2];
    for (k, inside) in [true, false].into_iter().enumerate() {
        let part = p
            .lamina_part(1, |b| b.contains(0) == inside)
            .expect("lamina part");
        if part.faces.is_empty() {
            continue;
        }
        assert_eq!(part.body_kind(), BodyKind::Sheet);
        let findings = part.validate_all();
        assert!(findings.is_empty(), "{findings:?}");
        for f in part.faces.values() {
            let n = f
                .surface
                .as_ref()
                .and_then(|s| s.folded_normal(0.0, 0.0))
                .expect("normal");
            assert!(n.z > 0.99, "the sheet's sense is kept: {n:?}");
        }
        out[k] = area(&part);
    }
    (out[0], out[1])
}

fn close(got: (f64, f64), want: (f64, f64)) {
    assert!(
        (got.0 - want.0).abs() < 1e-9 && (got.1 - want.1).abs() < 1e-9,
        "{got:?} != {want:?}"
    );
}

#[test]
fn a_sheet_through_a_box_parts_inside_and_outside() {
    close(parts(&sheet(-1.0, 2.0, -1.0, 2.0, 0.5)), (1.0, 8.0));
}

#[test]
fn a_sheet_ending_inside_a_box_parts_at_the_box() {
    close(parts(&sheet(-1.0, 0.5, -1.0, 2.0, 0.5)), (0.5, 4.0));
}

#[test]
fn a_sheet_inside_a_box_is_all_inside() {
    close(parts(&sheet(0.25, 0.75, 0.25, 0.75, 0.5)), (0.25, 0.0));
}

#[test]
fn a_sheet_missing_a_box_is_all_outside() {
    close(parts(&sheet(2.0, 3.0, 0.0, 1.0, 0.5)), (0.0, 1.0));
}

#[test]
fn a_solid_is_no_sheet_part() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let s = sheet(-1.0, 2.0, -1.0, 2.0, 0.5);
    let p = partition(&[&cube, &s], tol).expect("partition");
    assert!(matches!(
        p.lamina_part(0, |_| true),
        Err(crate::BoolError::DegenerateInput(_))
    ));
}
