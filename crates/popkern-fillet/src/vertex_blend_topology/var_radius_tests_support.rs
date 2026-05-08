#![allow(unused_imports)]

pub(super) use crate::FilletOptions;
pub(super) use crate::types::{FilletRadius, FilletSpec};
pub(super) use popkern_core::geometry::Point3;
pub(super) use popkern_core::tolerance::Tolerance;
pub(super) use popkern_core::topology::BRep;

pub(super) fn build_pyramid_4() -> (
    popkern_core::topology::BRep,
    popkern_core::topology::VertexHandle,
    [popkern_core::topology::FaceHandle; 4],
) {
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let mut brep = BRep::new(Tolerance::default());

    // Build square A(0,0) → B(1,0) → C(1,1) → D(0,1) in F_orig.
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap(); // A
    let mev_b = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap(); // B
    let mev_c = brep
        .mev(mev_b.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap(); // C
    let mev_d = brep
        .mev(mev_c.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap(); // D
    // mef D→A: closes the square. The "old face" (F_orig from
    // mvfs) keeps half-edges from new_he onward — i.e., the
    // interior square loop {DA, AB, BC, CD}.  The new face
    // (mef_da.face) gets the loop walking the opposite way
    // {AD, DC, CB, BA} — that's F_back.
    let mef_da = brep.mef(mev_d.half_edge_twin, mev_b.half_edge).unwrap();
    let f_inside = mvfs.face; // contains interior square loop
    let _f_back = mef_da.face; // back/outside face

    // Insert V_apex at (0.5, 0.5, 0) via mev from A's outgoing
    // (= mev_b.half_edge, A→B in F_inside).  The spur lives in
    // F_inside.
    let mev_v = brep
        .mev(mev_b.half_edge, Point3::new(0.5, 0.5, 0.0))
        .unwrap();
    let v_apex = mev_v.vertex;
    // mev_v.half_edge_twin is V→A (in F_inside)

    // mef V_apex → B: split F_inside into triangle {A,B,V} and pentagon {V,B,C,D,A}.
    // mev_v.half_edge_twin originates at V (V→A), and mev_c.half_edge originates at B (B→C).
    // Both in F_inside's loop after the spur insertion.
    let mef_apex_b = brep.mef(mev_v.half_edge_twin, mev_c.half_edge).unwrap();
    // mef_apex_b.half_edge is V→B (in F_inside, now pentagon)
    // mef_apex_b.face is the new triangle face (A,B,V).

    // mef V_apex → C: split pentagon into triangle {B,C,V} and quad {V,C,D,A}.
    // mef_apex_b.half_edge originates at V (V→B in pentagon).
    // mev_d.half_edge originates at C (C→D in pentagon).
    let mef_apex_c = brep.mef(mef_apex_b.half_edge, mev_d.half_edge).unwrap();
    // mef_apex_c.half_edge is V→C (in F_inside, now quad).
    // mef_apex_c.face is the new triangle face (B,C,V).

    // mef V_apex → D: split quad into triangle {C,D,V} and triangle {D,A,V}.
    // mef_apex_c.half_edge originates at V (V→C in quad).
    // mef_da.half_edge originates at D (D→A in quad).
    let mef_apex_d = brep.mef(mef_apex_c.half_edge, mef_da.half_edge).unwrap();
    // mef_apex_d.half_edge is V→D (in F_inside, now last triangle (D,A,V)).
    // mef_apex_d.face is the new triangle face (C,D,V).

    // Faces in cyclic order around V_apex going CCW (V→A is at
    // angle 225°, V→B at 315°, V→C at 45°, V→D at 135°):
    //   tri0 = (V,A,B) face  (angles 225°→315°)
    //   tri1 = (V,B,C) face  (angles 315°→45°)
    //   tri2 = (V,C,D) face  (angles 45°→135°)
    //   tri3 = (V,D,A) face  (angles 135°→225°)
    //
    //   tri0 = mef_apex_b.face (the (A,B,V) triangle face)
    //   tri1 = mef_apex_c.face (the (B,C,V) triangle face)
    //   tri2 = mef_apex_d.face (the (C,D,V) triangle face)
    //   tri3 = f_inside        (the leftover (D,A,V) triangle face)
    let tri0 = mef_apex_b.face;
    let tri1 = mef_apex_c.face;
    let tri2 = mef_apex_d.face;
    let tri3 = f_inside;
    (brep, v_apex, [tri0, tri1, tri2, tri3])
}
