//! PK-168 gate: the PLAIN skewed-corner chamfer — each bevel ending
//! at a body corner where the yawed+tilted boolean-cut rest plane
//! meets the end face obliquely (the laptop-stand arm-root class).

use popkern_part::{EdgeQuery, Plane, Profile};

/// PIN → GATE (PK-168, armed check 9): the PLAIN skewed-corner
/// chamfer — no chain junctions, each bevel ending at a body corner
/// where the yawed+tilted rest plane meets the end face obliquely
/// (the laptop-stand arm-root class). The historical end treatment
/// folded the strip loop through the UNCUT corner vertex (d/√2 off
/// the bevel plane); every mesh gate stayed green because the fold
/// closed the shell, and only the armed curve-surface validation
/// exposed it. The terminal chord trim now recognizes skewed planar
/// stubs (chart-direction gates are surface-scoped) and settles the
/// end to the rail-to-rail chord. Part::chamfer validates the
/// output, so this gate enforces check 9 end to end.
#[test]
fn plain_skewed_corner_chamfer_settles_terminal_chords() {
    let tilt = 18.0_f64.to_radians();
    let w = 12.0_f64.to_radians();
    let bar = Profile::on(Plane::XY)
        .rect(40.0, 60.0)
        .extrude(96.0)
        .expect("bar");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("upper")
        .rotate_x(tilt)
        .expect("rx")
        .rotate_z(-w)
        .expect("rz")
        .translate(0.0, 0.0, 80.0)
        .expect("place");
    let bar = bar.subtract(&upper).expect("rest cut");
    let rest_z = move |x: f64, y: f64| 80.0 + (x * w.sin() + y * w.cos()) * tilt.tan();
    let dressed = bar
        .chamfer(1.0, move |e: EdgeQuery| {
            e.straight().filter(move |i| {
                i.direction.map(|d| d.y.abs() > 0.85).unwrap_or(false)
                    && (i.midpoint.z - rest_z(i.midpoint.x, i.midpoint.y)).abs() < 0.5
            })
        })
        .expect("plain skewed chamfer must validate (check 9 armed)");
    assert_eq!(dressed.open_edge_count().expect("mesh"), 0);
    dressed
        .brep()
        .validate()
        .expect("plain skewed chamfer output validates");
}
