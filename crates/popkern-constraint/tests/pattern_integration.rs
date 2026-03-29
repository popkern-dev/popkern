//! Integration tests for pattern operations.

use popkern_constraint::{Constraint, ElementHandle, Sketch};
use std::f64::consts::PI;

const TOL: f64 = 1e-6;

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() < TOL
}

// ── 1. Linear pattern basic ──────────────────────────────────────────────

#[test]
fn linear_pattern_basic_points() {
    let mut s = Sketch::new();
    let p = s.add_point(0.0, 0.0);
    s.add_constraint(Constraint::Fixed(p));

    let ph = s
        .linear_pattern(&[ElementHandle::Point(p)], 4, (1.0, 0.0), 2.0)
        .unwrap();

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    for (i, inst) in pat.instances.iter().enumerate() {
        let k = (i + 1) as f64;
        match inst.elements[0] {
            ElementHandle::Point(iph) => {
                let (x, y) = s.point_position(iph).unwrap();
                assert!(
                    approx(x, k * 2.0),
                    "instance {i}: x={x}, expected {}",
                    k * 2.0
                );
                assert!(approx(y, 0.0), "instance {i}: y={y}, expected 0.0");
            }
            _ => panic!("expected point"),
        }
    }
}

// ── 2. Linear with internal constraints ──────────────────────────────────

#[test]
fn linear_pattern_with_horizontal() {
    let mut s = Sketch::new();
    let p1 = s.add_point(0.0, 0.0);
    let p2 = s.add_point(1.0, 0.1);
    let l = s.add_line(p1, p2);
    s.add_constraint(Constraint::Fixed(p1));
    s.add_constraint(Constraint::Horizontal(l));

    let ph = s
        .linear_pattern(
            &[
                ElementHandle::Point(p1),
                ElementHandle::Point(p2),
                ElementHandle::Line(l),
            ],
            3,
            (0.0, 1.0),
            3.0,
        )
        .unwrap();

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    for (i, inst) in pat.instances.iter().enumerate() {
        for elem in &inst.elements {
            if let ElementHandle::Line(lh) = elem {
                let (start, end) = s.line_endpoints(*lh).unwrap();
                let (_, y1) = s.point_position(start).unwrap();
                let (_, y2) = s.point_position(end).unwrap();
                assert!(
                    approx(y1, y2),
                    "instance {i}: line not horizontal (y1={y1}, y2={y2})"
                );
            }
        }
    }
}

// ── 3. Circular pattern ──────────────────────────────────────────────────

#[test]
fn circular_pattern_point() {
    let mut s = Sketch::new();
    let center = s.add_point(0.0, 0.0);
    let seed = s.add_point(1.0, 0.0);
    s.add_constraint(Constraint::Fixed(center));
    s.add_constraint(Constraint::Fixed(seed));

    let ph = s
        .circular_pattern(&[ElementHandle::Point(seed)], 4, center, 2.0 * PI)
        .unwrap();

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    let expected_angles: [f64; 3] = [2.0 * PI / 3.0, 4.0 * PI / 3.0, 2.0 * PI];
    for (i, (inst, angle)) in pat.instances.iter().zip(expected_angles.iter()).enumerate() {
        match inst.elements[0] {
            ElementHandle::Point(iph) => {
                let (x, y) = s.point_position(iph).unwrap();
                let ex = angle.cos();
                let ey = angle.sin();
                assert!(
                    approx(x, ex) && approx(y, ey),
                    "instance {i}: ({x}, {y}), expected ({ex}, {ey})"
                );
            }
            _ => panic!("expected point"),
        }
    }
}

// ── 4. Circular with arc ─────────────────────────────────────────────────

#[test]
fn circular_pattern_arc() {
    let mut s = Sketch::new();
    let center = s.add_point(0.0, 0.0);
    s.add_constraint(Constraint::Fixed(center));

    let pin = s.add_point(2.0, 0.0);
    let arc = s.add_arc(2.0, 0.0, 0.5, 0.0, PI / 2.0);
    s.add_constraint(Constraint::PinToCenter(pin, ElementHandle::Arc(arc)));
    s.add_constraint(Constraint::Fixed(pin));

    let ph = s
        .circular_pattern(&[ElementHandle::Arc(arc)], 3, center, PI)
        .unwrap();

    s.solve().unwrap();

    // Instance 0 is at angle pi/2 from seed.
    // Arc centre rotated from (2,0) by pi/2 around (0,0) = (0,2).
    // AngleOffsetArc enforces: arc2.start - arc1.start = theta,
    //                          arc2.end   - arc1.end   = theta.
    // The seed arc angles are free, so we verify the offset rather
    // than absolute values.
    let (_scx, _scy, _sr, seed_sa, seed_ea) = s.arc_params(arc).unwrap();
    let angle_k = PI / 2.0; // first instance angle

    let pat = s.pattern(ph).unwrap();
    let inst = &pat.instances[0];
    match inst.elements[0] {
        ElementHandle::Arc(ah) => {
            let (cx, cy, _r, sa, ea) = s.arc_params(ah).unwrap();
            assert!(
                approx(cx, 0.0) && approx(cy, 2.0),
                "centre ({cx}, {cy}), expected (0, 2)"
            );
            assert!(
                approx(sa - seed_sa, angle_k),
                "angle offset: sa diff={}, expected {angle_k}",
                sa - seed_sa
            );
            assert!(
                approx(ea - seed_ea, angle_k),
                "angle offset: ea diff={}, expected {angle_k}",
                ea - seed_ea
            );
        }
        _ => panic!("expected arc"),
    }
}

// ── 5. Mirror pattern ────────────────────────────────────────────────────

#[test]
fn mirror_pattern_point() {
    let mut s = Sketch::new();
    let seed = s.add_point(1.0, 0.0);
    let l1 = s.add_point(0.0, 0.0);
    let l2 = s.add_point(0.0, 1.0);
    let line = s.add_line(l1, l2);
    s.add_constraint(Constraint::Fixed(seed));
    s.add_constraint(Constraint::FixedLine(line));

    let ph = s
        .mirror_pattern(&[ElementHandle::Point(seed)], line)
        .unwrap();

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    match pat.instances[0].elements[0] {
        ElementHandle::Point(iph) => {
            let (x, y) = s.point_position(iph).unwrap();
            assert!(approx(x, -1.0), "reflected x={x}, expected -1.0");
            assert!(approx(y, 0.0), "reflected y={y}, expected 0.0");
        }
        _ => panic!("expected point"),
    }
}

// ── 6. Mirror arc ────────────────────────────────────────────────────────

#[test]
fn mirror_pattern_arc() {
    let mut s = Sketch::new();
    let arc = s.add_arc(2.0, 0.0, 1.0, 0.0, PI / 2.0);
    let l1 = s.add_point(0.0, 0.0);
    let l2 = s.add_point(0.0, 1.0);
    let line = s.add_line(l1, l2);
    s.add_constraint(Constraint::FixedLine(line));
    let pin = s.add_point(2.0, 0.0);
    s.add_constraint(Constraint::PinToCenter(pin, ElementHandle::Arc(arc)));
    s.add_constraint(Constraint::Fixed(pin));

    let ph = s.mirror_pattern(&[ElementHandle::Arc(arc)], line).unwrap();

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    match pat.instances[0].elements[0] {
        ElementHandle::Arc(ah) => {
            let (cx, cy, r, sa, ea) = s.arc_params(ah).unwrap();
            assert!(approx(cx, -2.0), "reflected cx={cx}, expected -2.0");
            assert!(approx(cy, 0.0), "reflected cy={cy}");
            assert!(approx(r, 1.0), "r={r}");
            // Mirror across y-axis: line_angle = pi/2
            // sa = 2*(pi/2) - (pi/2) = pi/2
            // ea = 2*(pi/2) - 0 = pi
            assert!(approx(sa, PI / 2.0), "sa={sa}, expected pi/2");
            assert!(approx(ea, PI), "ea={ea}, expected pi");
        }
        _ => panic!("expected arc"),
    }
}

// ── 7. Phantom centre points ─────────────────────────────────────────────

#[test]
fn phantom_center_points_linear_arc() {
    let mut s = Sketch::new();
    let arc = s.add_arc(1.0, 0.0, 0.5, 0.0, PI);
    let pin = s.add_point(1.0, 0.0);
    s.add_constraint(Constraint::PinToCenter(pin, ElementHandle::Arc(arc)));
    s.add_constraint(Constraint::Fixed(pin));

    let ph = s
        .linear_pattern(&[ElementHandle::Arc(arc)], 3, (1.0, 0.0), 3.0)
        .unwrap();

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    for (i, inst) in pat.instances.iter().enumerate() {
        assert!(
            !inst.phantom_points.is_empty(),
            "instance {i} has no phantom points"
        );
        assert!(
            !inst.pin_constraints.is_empty(),
            "instance {i} has no pin constraints"
        );
        let phantom = inst.phantom_points[0];
        let (px, py) = s.point_position(phantom).unwrap();
        match inst.elements[0] {
            ElementHandle::Arc(ah) => {
                let (cx, cy, _, _, _) = s.arc_params(ah).unwrap();
                assert!(
                    approx(px, cx) && approx(py, cy),
                    "instance {i}: phantom ({px}, {py}) != centre ({cx}, {cy})"
                );
            }
            _ => panic!("expected arc"),
        }
    }
}

// ── 8. Suppression ───────────────────────────────────────────────────────

#[test]
fn suppress_and_unsuppress() {
    let mut s = Sketch::new();
    let p = s.add_point(0.0, 0.0);
    s.add_constraint(Constraint::Fixed(p));

    let ph = s
        .linear_pattern(&[ElementHandle::Point(p)], 4, (1.0, 0.0), 2.0)
        .unwrap();

    s.solve().unwrap();

    // Suppress instance 1 (second copy, at x=4)
    s.suppress_instance(ph, 1).unwrap();
    s.solve().unwrap();

    // Unsuppress
    s.unsuppress_instance(ph, 1).unwrap();
    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    match pat.instances[1].elements[0] {
        ElementHandle::Point(iph) => {
            let (x, _y) = s.point_position(iph).unwrap();
            assert!(approx(x, 4.0), "after unsuppress: x={x}, expected 4.0");
        }
        _ => panic!("expected point"),
    }
}

// ── 9. Update count ──────────────────────────────────────────────────────

#[test]
fn update_count_increase_and_decrease() {
    let mut s = Sketch::new();
    let p = s.add_point(0.0, 0.0);
    s.add_constraint(Constraint::Fixed(p));

    let ph = s
        .linear_pattern(&[ElementHandle::Point(p)], 3, (1.0, 0.0), 2.0)
        .unwrap();

    s.solve().unwrap();
    assert_eq!(s.pattern(ph).unwrap().instances.len(), 2);

    // Increase to 5
    s.update_pattern_count(ph, 5).unwrap();
    s.solve().unwrap();
    assert_eq!(s.pattern(ph).unwrap().instances.len(), 4);

    // Verify instance at x=8 (k=4, 4*2=8)
    match s.pattern(ph).unwrap().instances[3].elements[0] {
        ElementHandle::Point(iph) => {
            let (x, _) = s.point_position(iph).unwrap();
            assert!(approx(x, 8.0), "new instance x={x}, expected 8.0");
        }
        _ => panic!("expected point"),
    }

    // Decrease to 2
    s.update_pattern_count(ph, 2).unwrap();
    assert_eq!(s.pattern(ph).unwrap().instances.len(), 1);
}

// ── 10. Refresh ──────────────────────────────────────────────────────────

#[test]
fn refresh_adds_new_seed_constraints() {
    let mut s = Sketch::new();
    let p1 = s.add_point(0.0, 0.0);
    let p2 = s.add_point(1.0, 0.1);
    let l = s.add_line(p1, p2);
    s.add_constraint(Constraint::Fixed(p1));

    let ph = s
        .linear_pattern(
            &[
                ElementHandle::Point(p1),
                ElementHandle::Point(p2),
                ElementHandle::Line(l),
            ],
            3,
            (0.0, 1.0),
            3.0,
        )
        .unwrap();

    // Add Horizontal to seed line AFTER pattern creation
    s.add_constraint(Constraint::Horizontal(l));

    let old_count = s.pattern(ph).unwrap().instances[0]
        .internal_constraints
        .len();

    s.refresh_pattern(ph).unwrap();

    let new_count = s.pattern(ph).unwrap().instances[0]
        .internal_constraints
        .len();
    assert!(
        new_count > old_count,
        "refresh should add constraints: old={old_count}, new={new_count}"
    );

    s.solve().unwrap();

    let pat = s.pattern(ph).unwrap();
    for (i, inst) in pat.instances.iter().enumerate() {
        for elem in &inst.elements {
            if let ElementHandle::Line(lh) = elem {
                let (start, end) = s.line_endpoints(*lh).unwrap();
                let (_, y1) = s.point_position(start).unwrap();
                let (_, y2) = s.point_position(end).unwrap();
                assert!(
                    approx(y1, y2),
                    "instance {i}: line not horizontal after refresh"
                );
            }
        }
    }
}

// ── 11. Cannot pattern instance element ──────────────────────────────────

#[test]
fn cannot_pattern_instance_element() {
    let mut s = Sketch::new();
    let p = s.add_point(0.0, 0.0);

    let ph = s
        .linear_pattern(&[ElementHandle::Point(p)], 3, (1.0, 0.0), 2.0)
        .unwrap();

    let instance_elem = s.pattern(ph).unwrap().instances[0].elements[0];
    let result = s.linear_pattern(&[instance_elem], 2, (0.0, 1.0), 1.0);
    assert!(
        result.is_err(),
        "should error when patterning instance element"
    );
}

// ── 12. Delete pattern ───────────────────────────────────────────────────

#[test]
fn delete_pattern_removes_all() {
    let mut s = Sketch::new();
    let p = s.add_point(0.0, 0.0);

    let ph = s
        .linear_pattern(&[ElementHandle::Point(p)], 3, (1.0, 0.0), 2.0)
        .unwrap();

    let points_before = s.point_count();
    let constraints_before = s.constraint_count();

    s.delete_pattern(ph).unwrap();

    assert!(s.point_count() < points_before, "points not removed");
    assert!(
        s.constraint_count() < constraints_before,
        "constraints not removed"
    );

    // Seed point should still exist
    assert!(s.point_position(p).is_some(), "seed point was removed");
}
