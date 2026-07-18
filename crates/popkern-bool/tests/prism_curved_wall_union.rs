//! PK-150 — prism-through-curved-wall union structuring (stage 1:
//! measurement + driving gate).
//!
//! Measured 2026-07-18 (via the PK-141 seam fixture): the union of a
//! rectangular prism piercing a cylinder wall validates clean but is
//! structurally garbage — the prism's side walls are missing, its
//! top/bottom planes appear twice, the wall is carved by FULL
//! circles instead of the boss footprint arcs, and single-edge loops
//! degenerate to points at the boss corners.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn fixture() -> (BRep, BRep) {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let outline = [
        Point3::new(0.6, -0.25, 0.75),
        Point3::new(1.3, -0.25, 0.75),
        Point3::new(1.3, 0.25, 0.75),
        Point3::new(0.6, 0.25, 0.75),
    ];
    let boss = BRep::prism(&outline, Vec3::new(0.0, 0.0, 1.0), 0.5, tol).unwrap();
    (cyl, boss)
}

/// Driving gate: the union must carry the boss's ±y side walls and
/// must not degenerate any edge to a point.
#[test]
fn union_preserves_boss_side_walls() {
    let (cyl, boss) = fixture();
    let united = popkern_bool::union(cyl, boss, Tolerance::default(), None).expect("union");
    // Boss side walls: planes with normal ±Y at y = ±0.25.
    let mut side_walls = 0usize;
    for (_, f) in united.faces.iter() {
        if let Some(sof) = f.surface.as_ref()
            && let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry
            && normal.y.abs() > 0.99
            && (origin.y.abs() - 0.25).abs() < 1e-9
        {
            side_walls += 1;
        }
    }
    assert_eq!(side_walls, 2, "boss ±y side walls must survive the union");
    // No point-degenerate edges.
    for (e, d) in united.edges.iter() {
        let he = &united.half_edges[d.half_edge];
        let a = united.vertices[he.origin].position;
        let b = united.vertices[united.half_edges[he.twin].origin].position;
        let f1 = united.loops[he.loop_].face;
        let f2 = united.loops[united.half_edges[he.twin].loop_].face;
        if f1 == f2 {
            continue; // seam edges legitimately close on themselves
        }
        let closed = (a - b).length() < 1e-9;
        if closed {
            // A closed non-seam edge must be a genuine ring (its
            // curve must sweep, not collapse): measure curve length
            // via chordal samples.
            let len = united.edges[e]
                .curve
                .as_ref()
                .map(|c| {
                    use popkern_core::geometry::CurveEval;
                    let mut l = 0.0;
                    let mut prev = c.geometry.point_at(c.t_start);
                    for k in 1..=16 {
                        let t = c.t_start + (c.t_end - c.t_start) * k as f64 / 16.0;
                        let p = c.geometry.point_at(t);
                        l += (p - prev).length();
                        prev = p;
                    }
                    l
                })
                .unwrap_or(0.0);
            assert!(
                len > 1e-6,
                "edge {e:?} is a point-degenerate loop (len {len:.3e})"
            );
        }
    }
}

/// Measurement: record the face census and which pipeline produced
/// it (run with PK_BOOL_TRACE=1 for pipeline attribution).
#[test]
fn union_census_measured() {
    let (cyl, boss) = fixture();
    let united = popkern_bool::union(cyl, boss, Tolerance::default(), None).expect("union");
    let mut planes = 0usize;
    let mut cylinders = 0usize;
    for (_, f) in united.faces.iter() {
        match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Cylinder { .. }) => cylinders += 1,
            _ => {}
        }
    }
    eprintln!(
        "[pk150] union census: {} faces ({planes} planes, {cylinders} cylinders), {} edges, {} vertices",
        united.faces.len(),
        united.edges.len(),
        united.vertices.len()
    );
    for (fh, f) in united.faces.iter() {
        let kind = match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { origin, normal, .. }) => format!(
                "Plane o({:.2},{:.2},{:.2}) n({:.2},{:.2},{:.2})",
                origin.x, origin.y, origin.z, normal.x, normal.y, normal.z
            ),
            Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder".into(),
            Some(_) => "other".into(),
            None => "NONE".into(),
        };
        eprintln!("[pk150]   face {fh:?}: {kind}");
    }
    let defects = united.validate_all();
    eprintln!("[pk150] validate_all: {} defects", defects.len());
    for d in defects.iter().take(6) {
        eprintln!("[pk150]   {d:?}");
    }
    // Correct structure: 3 cylinder-side faces (wall + 2 caps) and 5
    // boss faces (4 sides + outer cap; inner cap swallowed) = wall
    // carve pending; this print is the ledger evidence line.
}
