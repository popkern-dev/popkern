use popkern_core::naming::OperationId;

use crate::{concatenate_solids, name_from_operands, split};

use super::*;

/// A box `[x0, x0 + w] x [y0, y0 + d] x [z0, z0 + h]` whose faces are
/// named as operation `op`'s laterals, in slot order.
fn named_box(op: u64, at: (f64, f64, f64), size: (f64, f64, f64)) -> BRep {
    let tol = Tolerance::default();
    let (x1, y1) = (at.0 + size.0, at.1 + size.1);
    let ring = |z: f64| {
        [(at.0, at.1), (x1, at.1), (x1, y1), (at.0, y1)].map(|(x, y)| Point3::new(x, y, z))
    };
    let mut b = BRep::polygon_frustum(&ring(at.2), &ring(at.2 + size.2), tol).expect("box");
    let faces: Vec<FaceHandle> = b.faces.keys().collect();
    for (k, f) in faces.into_iter().enumerate() {
        let name = EntityName {
            operation: OperationId(op),
            role: EntityRole {
                category: RoleCategory::Lateral,
                index: k as u32,
            },
            ancestors: Vec::new(),
        };
        let id = b.name_pool.intern(name, None);
        b.faces[f].name = Some(id);
    }
    b
}

/// The plate cut by a slab at `x`, its pieces named as one under
/// `scheme`: the name of every face and edge, with where it lies in x
/// relative to the slab (left of, on or right of it) and which way a
/// face points in x, sorted.
fn names_at(x: f64, scheme: u32) -> Vec<(String, i8, i8)> {
    let tol = Tolerance::default();
    let plate = named_box(1, (-1.2, -0.8, 0.0), (2.4, 1.6, 0.4));
    let slab = named_box(2, (x - 0.2, -1.5, -0.5), (0.4, 3.0, 2.0));
    let cut = split(&plate, &[&slab], tol).expect("split");
    let mut pieces = cut.pieces.into_iter().map(|p| p.body);
    let first = pieces.next().expect("a piece");
    let mut body = pieces
        .try_fold(first, concatenate_solids)
        .expect("concatenated");
    let mut ctx = TrackingContext::new(OperationId(9)).with_naming_scheme(scheme);
    name_from_operands(&mut body, &plate, &[&slab], &tol, &mut ctx);
    let side = |px: f64| {
        if px < x - 0.2 - 1e-6 {
            -1
        } else if px > x + 0.2 + 1e-6 {
            1
        } else {
            0
        }
    };
    let spelled = |b: &BRep, id: NameId| format!("{:?}", b.name_pool.get(id));
    let mut out = Vec::new();
    for (f, face) in &body.faces {
        let fc = capture_face(&body, f).expect("face");
        let id = face.name.expect("named");
        out.push((
            spelled(&body, id),
            side(fc.centroid.x),
            fc.normal.x.round() as i8,
        ));
    }
    for (e, edge) in &body.edges {
        let (mid, _) = mid_tangent(&body, e).expect("curve");
        let id = edge.name.expect("named");
        out.push((spelled(&body, id), side(mid.x), 2));
    }
    out.sort();
    out
}

/// Under scheme 5 the plate's pieces carry the same names wherever the
/// slab cuts them: every face and edge is named by what it is a piece of
/// and by its neighbours. Under scheme 4 (the negative control) moving
/// the slab renames faces the cut made, so a reference recorded on one
/// no longer resolves.
#[test]
fn a_moved_tool_keeps_every_name_under_scheme_5() {
    let here = names_at(0.0, 5);
    let there = names_at(0.5, 5);
    assert_eq!(here, there);
    let distinct: BTreeSet<_> = here.iter().map(|(n, ..)| n).collect();
    assert_eq!(
        distinct.len(),
        here.len(),
        "every face and edge named apart"
    );
    assert_ne!(names_at(0.0, 4), names_at(0.5, 4), "scheme 4 renames them");
}
