//! Membership of a shell class no face of an operand bounds is the
//! classification of a point on one of its faces: an offset probe can
//! cross a gap thinner than the offset into the other operand.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

/// A 100 × 100 × 10 plate bored through by an r = 10 hole, and an
/// r = 9.998 pin through the bore from z = −10 to 20: a 0.002 radial
/// gap, so the bodies share no material. The pin's wall crosses the
/// plate's top and bottom planes inside the bore (records that split
/// only the pin).
fn plate_and_pin(tol: Tolerance) -> (BRep, BRep) {
    let plate = BRep::prism(
        &[
            Point3::new(-50.0, -50.0, 0.0),
            Point3::new(50.0, -50.0, 0.0),
            Point3::new(50.0, 50.0, 0.0),
            Point3::new(-50.0, 50.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        tol,
    )
    .expect("plate");
    let bore = BRep::cylinder(10.0, 12.0, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -1.0))
        .expect("placed bore");
    let plate =
        crate::boolean(BoolOperation::Subtraction, plate, bore, tol, None).expect("bored plate");
    let pin = BRep::cylinder(9.998, 30.0, tol)
        .expect("pin")
        .transform(&Transform3::translation(0.0, 0.0, -10.0))
        .expect("placed pin");
    (plate, pin)
}

#[test]
fn a_pin_clear_of_its_bore_intersects_to_nothing() {
    let tol = Tolerance::default();
    let (plate, pin) = plate_and_pin(tol);
    let ix = crate::intersect::find_intersections(&plate, &pin, &tol).expect("intersections");
    let served =
        crate::complex::serve::rescue(&plate, &pin, &ix, BoolOperation::Intersection, &tol)
            .expect("manifold");
    let (body, _) = served.expect("the complex serves the intersection");
    let boxes = crate::intersect::compute_face_aabbs(&body);
    for p in [
        Point3::new(0.0, 0.0, 5.0),
        Point3::new(9.99, 0.0, 5.0),
        Point3::new(30.0, 30.0, 5.0),
    ] {
        let got = crate::classify::point_in_solid_final(&p, &body, &boxes, &tol);
        assert!(
            !matches!(got, Ok(crate::PointClassification::Inside)),
            "no material at {p:?}: {} faces",
            body.faces.len()
        );
    }
}

/// A unit box with an r = 0.5 rod lying along its top face's centre
/// line, the rod's end caps flush with the box's x = 0 and x = 1 faces:
/// each cap's rim circle touches the box's top edge at one point.
fn rod_on_box(place: &Transform3, tol: Tolerance) -> (BRep, BRep) {
    let block = BRep::unit_cube(tol)
        .expect("cube")
        .transform(place)
        .expect("placed cube");
    let rod = BRep::cylinder(0.5, 1.0, tol)
        .expect("rod")
        .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
        .expect("along x")
        .transform(&Transform3::translation(0.0, 0.5, 1.5))
        .expect("resting")
        .transform(place)
        .expect("placed rod");
    (block, rod)
}

/// The rim circles touch the end faces' top edges tangentially, so no
/// arc of them is a boundary of an end face: the complex intersects
/// the touching bodies to nothing and reports the union's contacts,
/// instead of declining on a face whose two sides joined (an end face
/// that took the whole rim circle as a hole, and a phantom disc of it).
#[test]
fn a_rod_resting_on_a_box_with_flush_ends_takes_the_complex() {
    let tol = Tolerance::default();
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.3, -0.2, 0.1),
        Vec3::new(1.0, 2.0, 3.0),
        0.37,
    )
    .expect("rotation");
    for place in [Transform3::identity(), turned] {
        let (block, rod) = rod_on_box(&place, tol);
        let ix = crate::intersect::find_intersections(&block, &rod, &tol).expect("intersections");
        let merge = |op| {
            crate::cellular::boolean::cellular_graph_merge(
                &block,
                &rod,
                &ix,
                op,
                &tol,
                None,
                crate::complex::serve::SelectionSource::Complex,
            )
        };
        let (body, _) = merge(BoolOperation::Intersection)
            .unwrap_or_else(|d| panic!("the complex intersects: {}", d.into_error()));
        assert!(body.faces.is_empty(), "{} faces", body.faces.len());
        let union = merge(BoolOperation::Union).map(|_| ());
        let why = union
            .err()
            .map(|d| d.into_error().to_string())
            .unwrap_or_default();
        assert!(why.contains("non-manifold"), "union: {why}");
    }
}

/// A 20 × 20 × 4 plate with an r = 5 round notch in its x = 10 edge,
/// and an r = 1.6 pin from z = `pin_z.0` to `pin_z.1` standing in the
/// notch `depth` proud of its wall at x = 5: their intersection is a
/// lens `depth` deep. The notch is part of the plate faces' outer
/// loops, as a cycloid disc's window flank is.
fn pin_proud_of_notch(
    depth: f64,
    pin_z: (f64, f64),
    place: &Transform3,
    tol: Tolerance,
) -> (BRep, BRep) {
    let plate = BRep::prism(
        &[
            Point3::new(-10.0, -10.0, 0.0),
            Point3::new(10.0, -10.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(-10.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        4.0,
        tol,
    )
    .expect("plate");
    let notch = BRep::cylinder(5.0, 6.0, tol)
        .expect("notch")
        .transform(&Transform3::rotation_z(std::f64::consts::PI))
        .expect("seam at the lens")
        .transform(&Transform3::translation(10.0, 0.0, -1.0))
        .expect("through");
    let plate = crate::boolean(BoolOperation::Subtraction, plate, notch, tol, None)
        .expect("notched plate")
        .transform(place)
        .expect("placed plate");
    let pin = BRep::cylinder(1.6, pin_z.1 - pin_z.0, tol)
        .expect("pin")
        .transform(&Transform3::translation(
            10.0 - (5.0 - 1.6 + depth),
            0.0,
            pin_z.0,
        ))
        .expect("proud")
        .transform(place)
        .expect("placed pin");
    (plate, pin)
}

/// A shallow lens with the notch wall's seam vertex at its middle: the
/// pin's cap circles cross the notch arc twice each and the complex
/// intersects the bodies to the lens, under a general placement too.
/// (The cycloid disc's flank, many short edges, gives the cap circles
/// crossing noise at every junction; that case is pinned in
/// `popkern-rhai`'s `a_rotor_geared_at_minus_n_interferes`, PK-1022.)
#[test]
fn a_pin_proud_of_its_notch_wall_intersects_to_the_lens() {
    let tol = Tolerance::default();
    let depth = 2.7e-3;
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.3, -0.2, 0.1),
        Vec3::new(1.0, 2.0, 3.0),
        0.37,
    )
    .expect("rotation");
    // A pin through the plate, and one flush with both its faces
    // (PK-1023: the plate's faces were not split by the flush pin's
    // coplanar rims and radial ordering tied them).
    for (place, pin_z) in [
        (Transform3::identity(), (-1.0, 5.0)),
        (turned, (-1.0, 5.0)),
        (Transform3::identity(), (0.0, 4.0)),
        (turned, (0.0, 4.0)),
    ] {
        let (plate, pin) = pin_proud_of_notch(depth, pin_z, &place, tol);
        let ix = crate::intersect::find_intersections(&plate, &pin, &tol).expect("intersections");
        let (mut body, _) = crate::cellular::boolean::cellular_graph_merge(
            &plate,
            &pin,
            &ix,
            BoolOperation::Intersection,
            &tol,
            None,
            crate::complex::serve::SelectionSource::Complex,
        )
        .unwrap_or_else(|d| panic!("the complex intersects: {}", d.into_error()));
        body.validate().expect("valid lens");
        let _ = body.recompute_pcurves_where_missing();
        let boxes = crate::intersect::compute_face_aabbs(&body);
        for (x, inside) in [
            (5.0 - 0.5 * depth, true),
            (5.0 + 0.5 * depth, false),
            (5.0 - 1.5 * depth, false),
        ] {
            let p = place.apply_point(&Point3::new(x, 0.0, 2.0));
            let got =
                crate::classify::point_in_solid_final(&p, &body, &boxes, &tol).expect("classifies");
            let want = if inside {
                crate::PointClassification::Inside
            } else {
                crate::PointClassification::Outside
            };
            assert_eq!(got, want, "x = {x}: {} faces", body.faces.len());
        }
    }
}

fn block(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    BRep::prism(
        &[
            Point3::new(lo.0, lo.1, lo.2),
            Point3::new(hi.0, lo.1, lo.2),
            Point3::new(hi.0, hi.1, lo.2),
            Point3::new(lo.0, hi.1, lo.2),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        hi.2 - lo.2,
        tol,
    )
    .expect("block")
}

/// A block slotted through its x = 0 end, and a block flush against
/// that end: the second block's x = 0 face meets the first's bottom
/// face along the bottom's x = 0 boundary line, which the slot breaks
/// for z in [1, 3]. The record along that line touches the bottom face's
/// boundary only at the slot corners; the slot mouth is a gap, not part
/// of the face (PK-1025: a boss cap flush with a notched web end served
/// the mouth as a slit on the web's bottom, and the complex declined).
/// The touching blocks intersect to nothing.
#[test]
fn a_flush_block_across_a_slot_mouth_intersects_to_nothing() {
    let tol = Tolerance::default();
    let slotted = crate::boolean(
        BoolOperation::Subtraction,
        block((-2.0, 0.0, 0.0), (0.0, 1.0, 4.0), tol),
        block((-1.0, -0.5, 1.0), (0.5, 1.5, 3.0), tol),
        tol,
        None,
    )
    .expect("slotted block");
    let flush = block((0.0, -0.5, -0.5), (1.0, 1.5, 4.5), tol);
    let ix = crate::intersect::find_intersections(&slotted, &flush, &tol).expect("intersections");
    let (body, _) = crate::cellular::boolean::cellular_graph_merge(
        &slotted,
        &flush,
        &ix,
        BoolOperation::Intersection,
        &tol,
        None,
        crate::complex::serve::SelectionSource::Complex,
    )
    .unwrap_or_else(|d| panic!("the complex intersects: {}", d.into_error()));
    assert!(body.faces.is_empty(), "{} faces", body.faces.len());
}

/// An r = 6.5 socket head with a hex socket 4.4 deep (side 3.46) in
/// one cap, and a block whose walls cut that cap's rim circle. `top`
/// picks the cap: the cylinder's two caps walk their rims in opposite
/// senses against the rim curve.
fn socket_head_and_block(top: bool, tol: Tolerance) -> (BRep, BRep) {
    use popkern_core::geometry::PortableMath;
    let head = BRep::cylinder(6.5, 8.0, tol)
        .expect("head")
        .transform(&Transform3::translation(0.0, 0.0, -8.0))
        .expect("placed head");
    // The socket and block from `z0`, spanning the cap at z = 0 or -8.
    let cap = if top { 0.0 } else { -8.0 };
    let socket_z0 = if top { cap - 4.4 } else { cap - 1.0 };
    let hex: Vec<Point3> = (0..6)
        .map(|k| {
            let a = std::f64::consts::FRAC_PI_3 * k as f64 - std::f64::consts::FRAC_PI_6;
            Point3::new(
                3.464_101_615 * a.pcos(),
                3.464_101_615 * a.psin(),
                socket_z0,
            )
        })
        .collect();
    let socket = BRep::prism(&hex, Vec3::new(0.0, 0.0, 1.0), 5.4, tol).expect("socket");
    let head =
        crate::boolean(BoolOperation::Subtraction, head, socket, tol, None).expect("socketed head");
    let block_z0 = if top { cap - 1.0 } else { cap - 2.0 };
    let block = BRep::prism(
        &[
            Point3::new(-1.833, -6.236, block_z0),
            Point3::new(10.0, -6.236, block_z0),
            Point3::new(10.0, 6.236, block_z0),
            Point3::new(-1.833, 6.236, block_z0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        3.0,
        tol,
    )
    .expect("block");
    (head, block)
}

/// A block's walls split a socket head's cap rim into arcs: the complex
/// serves the cut and the overlap, the socket kept as the cap's hole,
/// on either cap (they walk their rims in opposite senses). The arm's
/// socket-head screws lost the hole when the rim's arcs were sampled
/// against their half-edges; that measure is pinned in
/// `walk_loop_measures::tests`.
#[test]
fn a_socket_head_cut_across_its_rim_keeps_the_socket_as_a_hole() {
    let tol = Tolerance::default();
    for (top, op) in [true, false]
        .into_iter()
        .flat_map(|t| [BoolOperation::Subtraction, BoolOperation::Intersection].map(|op| (t, op)))
    {
        let (head, block) = socket_head_and_block(top, tol);
        let ix = crate::intersect::find_intersections(&head, &block, &tol).expect("intersections");
        let (body, _) = crate::cellular::boolean::cellular_graph_merge(
            &head,
            &block,
            &ix,
            op,
            &tol,
            None,
            crate::complex::serve::SelectionSource::Complex,
        )
        .unwrap_or_else(|d| panic!("the complex serves {op:?} (top {top}): {}", d.into_error()));
        body.validate()
            .unwrap_or_else(|e| panic!("{op:?} (top {top}): {e}"));
    }
}

/// Each face's surface kind, sorted vertex positions and orientation
/// flag.
fn oriented_faces(b: &BRep) -> Vec<(String, Vec<[u64; 3]>, bool)> {
    let mut out: Vec<_> = b
        .faces
        .iter()
        .map(|(fh, f)| {
            let s = f.surface.as_ref().expect("surface");
            let mut pts: Vec<[u64; 3]> = b
                .face_outer_loop_positions(fh)
                .iter()
                .map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
                .collect();
            pts.sort_unstable();
            (s.geometry.kind_name().to_string(), pts, s.orientation)
        })
        .collect();
    out.sort();
    out
}

/// A body of two solids that fold different ways (a cylinder and a
/// block united apart), cut at the block: the complex keeps each
/// solid's own folded-normal convention, so every face serves with
/// the orientation flag the cellular assembler gives it. One global
/// convention (the cylinder's) flipped the block's faces and the
/// pocket's walls.
#[test]
fn each_solid_of_a_two_solid_operand_keeps_its_face_orientation() {
    use crate::complex::serve::SelectionSource;
    let tol = Tolerance::default();
    let disc = BRep::cylinder(3.0, 1.5, tol).expect("disc");
    let far = block((100.0, 0.0, 0.0), (104.0, 4.0, 4.0), tol);
    let two = crate::boolean(BoolOperation::Union, disc, far, tol, None).expect("two solids");
    assert_eq!(two.solid_regions().len(), 2);
    let senses: std::collections::BTreeSet<_> = two
        .faces
        .keys()
        .map(|f| format!("{:?}", two.component_material_sense(f)))
        .collect();
    assert_eq!(
        senses.len(),
        2,
        "the solids fold different ways: {senses:?}"
    );
    let pocket = block((101.0, 1.0, -1.0), (103.0, 3.0, 5.0), tol);
    let ix = crate::intersect::find_intersections(&two, &pocket, &tol).expect("intersections");
    let (complex, _) =
        crate::complex::serve::rescue(&two, &pocket, &ix, BoolOperation::Subtraction, &tol)
            .expect("manifold")
            .expect("the complex serves");
    let (cellular, _) = crate::cellular::boolean::cellular_graph_merge(
        &two,
        &pocket,
        &ix,
        BoolOperation::Subtraction,
        &tol,
        None,
        SelectionSource::Cellular,
    )
    .expect("the cellular assembler serves");
    assert_eq!(oriented_faces(&complex), oriented_faces(&cellular));
}

/// A union keeping a tool clear of the part (a block beside a
/// cylinder; they fold different ways): the tool's faces form an output
/// shell with no face of A, and keep the tool's own convention. A's
/// convention applied to it flipped every tool face.
#[test]
fn a_tool_shell_with_no_face_of_a_keeps_the_tool_s_orientation() {
    use crate::complex::serve::SelectionSource;
    let tol = Tolerance::default();
    let part = BRep::cylinder(3.0, 1.5, tol).expect("part");
    let tool = block((10.0, 0.0, 0.0), (14.0, 4.0, 4.0), tol);
    let ix = crate::intersect::find_intersections(&part, &tool, &tol).expect("intersections");
    let (complex, _) = crate::complex::serve::rescue(&part, &tool, &ix, BoolOperation::Union, &tol)
        .expect("manifold")
        .expect("the complex serves");
    let (cellular, _) = crate::cellular::boolean::cellular_graph_merge(
        &part,
        &tool,
        &ix,
        BoolOperation::Union,
        &tol,
        None,
        SelectionSource::Cellular,
    )
    .expect("the cellular assembler serves");
    assert_eq!(complex.faces.len(), part.faces.len() + tool.faces.len());
    assert_eq!(oriented_faces(&complex), oriented_faces(&cellular));
}
