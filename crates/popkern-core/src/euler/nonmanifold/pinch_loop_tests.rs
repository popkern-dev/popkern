use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A cylinder's lateral loop passes its bottom seam vertex twice (along
/// its bottom arc and up its seam). Parting it there gives two loops;
/// joining them again gives back the one loop, every fin linked both
/// ways and the face's loop count restored.
#[test]
fn a_loop_parted_and_joined_at_a_vertex_is_restored() {
    let mut c = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    let (face, l, v) = c
        .faces
        .iter()
        .find_map(|(f, face)| {
            let l = face.outer_loop;
            let fins = c.loop_fins(l);
            fins.iter()
                .map(|&h| c.half_edges[h].origin)
                .find(|&v| {
                    fins.iter()
                        .filter(|&&h| c.half_edges[h].origin == v)
                        .count()
                        == 2
                })
                .map(|v| (f, l, v))
        })
        .expect("a loop through a vertex twice");
    let before = c.loop_fins(l);
    let new = c.split_loop_at_vertex(face, l, v).expect("parts");
    assert_eq!(c.faces[face].inner_loops, vec![new]);
    let (a, b) = (c.loop_fins(l), c.loop_fins(new));
    assert_eq!(a.len() + b.len(), before.len());
    for &h in a.iter().chain(&b) {
        let he = &c.half_edges[h];
        assert_eq!(c.half_edges[he.next].prev, h, "next/prev agree");
    }
    c.join_loops_at_vertex(face, l, new, v).expect("joins");
    assert!(c.faces[face].inner_loops.is_empty());
    let after = c.loop_fins(l);
    let mut x = before.clone();
    let mut y = after.clone();
    x.sort();
    y.sort();
    assert_eq!(x, y, "the same fins");
    for &h in &after {
        assert_eq!(c.half_edges[c.half_edges[h].next].prev, h);
        assert_eq!(c.half_edges[h].loop_, l);
    }
    assert!(c.loops.get(new).is_none());
}

/// Joining declines, changing nothing, when a loop misses the vertex.
#[test]
fn joining_declines_when_a_loop_misses_the_vertex() {
    let mut c = BRep::unit_cube(Tolerance::default()).expect("cube");
    let (f, face) = c
        .faces
        .iter()
        .next()
        .map(|(f, x)| (f, x.clone()))
        .expect("a face");
    let other = c
        .faces
        .iter()
        .find(|(g, _)| *g != f)
        .map(|(_, x)| x.outer_loop)
        .expect("another face");
    let v = c.half_edges[c.loops[face.outer_loop].half_edge].origin;
    let snapshot = format!("{:?}", c.faces[f]);
    assert!(
        c.join_loops_at_vertex(f, face.outer_loop, other, v)
            .is_err()
    );
    assert_eq!(format!("{:?}", c.faces[f]), snapshot);
}

/// The pass finds the two loops a parting left at their shared vertex
/// and joins them back, once; a body with no such face is untouched.
#[test]
fn the_pass_joins_loops_sharing_a_vertex() {
    let mut c = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    assert_eq!(c.join_loops_sharing_vertices(), 0);
    let (face, l, v) = c
        .faces
        .iter()
        .find_map(|(f, face)| {
            let l = face.outer_loop;
            let fins = c.loop_fins(l);
            fins.iter()
                .map(|&h| c.half_edges[h].origin)
                .find(|&v| {
                    fins.iter()
                        .filter(|&&h| c.half_edges[h].origin == v)
                        .count()
                        == 2
                })
                .map(|v| (f, l, v))
        })
        .expect("a loop through a vertex twice");
    c.split_loop_at_vertex(face, l, v).expect("parts");
    assert_eq!(c.join_loops_sharing_vertices(), 1);
    assert!(c.faces[face].inner_loops.is_empty());
}
