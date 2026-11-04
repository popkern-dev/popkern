use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn a_solid_has_no_mirror_loops() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    for (h, f) in &cube.faces {
        assert!(!cube.is_mirror_loop(h, f.outer_loop));
    }
}

#[test]
fn a_one_face_sheet_holds_its_boundary_twice() {
    // The two-face quad's backside loop moved onto the front: the
    // single-face sheet form.
    let mut quad = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, Tolerance::default()).expect("quad");
    let mut faces: Vec<_> = quad.faces.keys().collect();
    faces.sort();
    let (front, back) = (faces[0], faces[1]);
    let ring = quad.faces[back].outer_loop;
    quad.loops[ring].face = front;
    quad.faces[front].inner_loops.push(ring);
    quad.faces.remove(back);
    assert!(quad.is_mirror_loop(front, ring));
    assert!(!quad.is_mirror_loop(front, quad.faces[front].outer_loop));
}
