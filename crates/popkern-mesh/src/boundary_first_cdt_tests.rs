use super::super::boundary_first_delaunay::bowyer_watson;
use super::*;

fn has_edge(tris: &[[usize; 3]], a: usize, b: usize) -> bool {
    tris.iter()
        .any(|t| (0..3).any(|k| key(t[k], t[(k + 1) % 3]) == key(a, b)))
}

fn area(pts: &[(f64, f64)], tris: &[[usize; 3]]) -> f64 {
    tris.iter()
        .map(|t| 0.5 * orient(pts[t[0]], pts[t[1]], pts[t[2]]))
        .sum()
}

/// A long segment across a column of points that sit inside its
/// diametral disk: Bowyer–Watson drops it, recovery restores it
/// without changing the triangulated area, every triangle stays
/// counter-clockwise, and the edges off the segment stay Delaunay.
#[test]
fn a_segment_crossed_by_encroachers_is_recovered_by_flips() {
    let mut pts = vec![(0.0, 0.0), (10.0, 0.1)];
    for k in 0..9 {
        let x = 1.0 + k as f64;
        pts.push((x, 0.5 + 0.03 * k as f64));
        pts.push((x + 0.5, -0.4));
    }
    pts.extend([(5.0, 4.0), (5.0, -4.0)]);
    let tris = bowyer_watson(&pts).expect("delaunay");
    assert!(
        !has_edge(&tris, 0, 1),
        "the fixture must encroach the segment"
    );
    let recovered = recover_segments(&pts, &tris, &[(0, 1)]).expect("recovered");
    assert!(has_edge(&recovered, 0, 1));
    assert_eq!(recovered.len(), tris.len());
    assert!((area(&pts, &recovered) - area(&pts, &tris).abs()).abs() < 1e-9);
    assert!(
        recovered
            .iter()
            .all(|t| orient(pts[t[0]], pts[t[1]], pts[t[2]]) > 0.0)
    );
    let mesh = Mesh::new(&pts, &recovered).expect("mesh");
    for &(p, q) in mesh.edges.keys() {
        if key(p, q) == key(0, 1) {
            continue;
        }
        if let Some((r, s, _, _)) = mesh.quad(p, q) {
            assert!(
                !in_circle(pts[p], pts[q], pts[r], pts[s]),
                "edge {p}-{q} is not locally Delaunay"
            );
        }
    }
}

/// A point on the segment's interior (a T-point) is the conforming
/// split's job: recovery declines rather than flip around it.
#[test]
fn a_point_on_the_segment_declines_recovery() {
    let pts = vec![(0.0, 0.0), (4.0, 0.0), (2.0, 0.0), (2.0, 1.0), (2.0, -1.0)];
    let tris = bowyer_watson(&pts).expect("delaunay");
    assert!(recover_segments(&pts, &tris, &[(0, 1)]).is_none());
}

/// Two crossing segments cannot both be edges: recovery declines.
#[test]
fn crossing_segments_decline_recovery() {
    let pts = vec![(0.0, 0.0), (4.0, 0.0), (2.0, -2.0), (2.0, 2.0), (2.0, 0.3)];
    let tris = bowyer_watson(&pts).expect("delaunay");
    assert!(recover_segments(&pts, &tris, &[(0, 1), (2, 3)]).is_none());
}

/// Segments already present leave the triangulation as it was.
#[test]
fn present_segments_leave_the_triangulation_unchanged() {
    let pts = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.1)];
    let tris = bowyer_watson(&pts).expect("delaunay");
    let (a, b) = (tris[0][0], tris[0][1]);
    let recovered = recover_segments(&pts, &tris, &[(a, b)]).expect("present");
    let canon = |ts: &[[usize; 3]]| {
        let mut v: Vec<[usize; 3]> = ts
            .iter()
            .map(|t| {
                let mut s = *t;
                s.sort_unstable();
                s
            })
            .collect();
        v.sort_unstable();
        v
    };
    assert_eq!(canon(&recovered), canon(&tris));
}

/// Incremental insertion keeps every constraint, the covered area and
/// counter-clockwise winding, leaves unconstrained edges locally
/// Delaunay, and declines a point on a constraint.
#[test]
fn steiner_insertion_keeps_constraints_and_the_empty_circle_property() {
    let mut pts = vec![(0.0, 0.0), (10.0, 0.1), (10.0, 4.0), (0.0, 4.2)];
    for k in 0..9 {
        pts.push((1.0 + k as f64, 0.5 + 0.03 * k as f64));
    }
    let tris = bowyer_watson(&pts).expect("delaunay");
    let segments = [(0, 1), (1, 2), (2, 3), (3, 0)];
    let mut cdt = ConstrainedMesh::new(&pts, &tris, &segments).expect("mesh");
    let before = area(&pts, &cdt.triangles());
    let inserted: Vec<(f64, f64)> = (0..20)
        .map(|k| (0.5 + 0.45 * k as f64, 0.2 + 0.17 * k as f64 % 3.5))
        .collect();
    for &p in &inserted {
        let v = cdt.insert(p).expect("inserted");
        assert_eq!(v, pts.len());
        pts.push(p);
    }
    // A point on the constrained bottom edge is the conforming split's.
    assert!(cdt.insert((5.0, 0.05)).is_none());
    let tris = cdt.triangles();
    for &(a, b) in &segments {
        assert!(has_edge(&tris, a, b), "constraint {a}-{b} kept");
    }
    assert!((area(&pts, &tris) - before).abs() < 1e-9);
    assert!(
        tris.iter()
            .all(|t| orient(pts[t[0]], pts[t[1]], pts[t[2]]) > 0.0)
    );
    let mesh = Mesh::new(&pts, &tris).expect("mesh");
    for &(p, q) in mesh.edges.keys() {
        if segments.iter().any(|&(a, b)| key(a, b) == key(p, q)) {
            continue;
        }
        if let Some((r, s, _, _)) = mesh.quad(p, q) {
            assert!(
                !in_circle(pts[p], pts[q], pts[r], pts[s]),
                "edge {p}-{q} is not locally Delaunay"
            );
        }
    }
}
