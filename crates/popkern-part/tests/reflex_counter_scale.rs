//! A counter through a reflex glyph prism is cut at every scale.
//!
//! The clip router sends reflex faces to the UV-space clipper (the
//! Cyrus–Beck line clipper is convex-only). Its reflex test compared
//! `(e1 × e2) · Newell` — a size^4 quantity — against an absolute
//! 1e-9, so at the dice's millimetre glyphs every turn read collinear,
//! the "4" and "8" outlines read convex, the counter walls were
//! clipped away, the subtraction found no records, and the counters
//! (the islands in 4 and 8) were never cut. Scale sweep: the same
//! cut at 1e-3, 1 and 1e3 must give the same topology and volume.
use popkern_part::{Part, Plane, Profile};

/// Inter Bold "4", cap height 1: the reflex outline and its counter.
const OUTLINE: [(f64, f64); 11] = [
    (0.5463, 0.1779),
    (0.5463, 0.0000),
    (0.7430, 0.0000),
    (0.7430, 0.1779),
    (0.8678, 0.1779),
    (0.8678, 0.3430),
    (0.7430, 0.3430),
    (0.7430, 1.0000),
    (0.4859, 1.0000),
    (0.0678, 0.3396),
    (0.0678, 0.1779),
];
const COUNTER: [(f64, f64); 3] = [(0.2745, 0.3510), (0.5503, 0.7745), (0.5503, 0.3430)];

fn shoelace(poly: &[(f64, f64)]) -> f64 {
    let n = poly.len();
    0.5 * (0..n)
        .map(|i| {
            let ((x0, y0), (x1, y1)) = (poly[i], poly[(i + 1) % n]);
            x0 * y1 - x1 * y0
        })
        .sum::<f64>()
        .abs()
}

fn signed_mesh_volume(part: &Part) -> f64 {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

#[test]
fn a_glyph_counter_is_cut_at_every_scale() {
    for k in [7e-3, 7.0, 7e3] {
        let scaled = |poly: &[(f64, f64)]| -> Vec<(f64, f64)> {
            poly.iter().map(|&(x, y)| (k * x, k * y)).collect()
        };
        let h = 0.17 * k;
        let ink = Profile::on(Plane::XY)
            .polygon(&scaled(&OUTLINE))
            .extrude(h)
            .expect("ink");
        let counter = Profile::on(Plane::XY)
            .polygon(&scaled(&COUNTER))
            .extrude(3.0 * h)
            .expect("counter")
            .translate(0.0, 0.0, -h)
            .expect("counter placement");
        let cut = ink.subtract(&counter).expect("counter cut");
        assert_eq!(
            cut.face_count(),
            OUTLINE.len() + 2 + COUNTER.len(),
            "scale {k}: outline walls, two caps and the counter's walls"
        );
        assert_eq!(
            cut.open_edge_count().expect("mesh"),
            0,
            "scale {k}: watertight"
        );
        let want = (shoelace(&OUTLINE) - shoelace(&COUNTER)) * k * k * h;
        let got = signed_mesh_volume(&cut);
        assert!(
            (got - want).abs() < 1e-9 * want,
            "scale {k}: volume {got:e}, prism with its counter removed {want:e}"
        );
    }
}
