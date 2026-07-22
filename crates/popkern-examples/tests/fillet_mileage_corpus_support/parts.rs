//! Corpus part builders for the fillet mileage suite — one
//! representative CAD part per fillet/boolean interaction class.
//!
//! Split from `fillet_mileage_corpus.rs` (via `#[path]`, in a
//! subdirectory so cargo does not compile it as its own test crate)
//! to keep the suite under the structural 500-line cap.

use popkern_part::prelude::*;

pub(crate) fn build_mounting_bracket() -> Result<Part, PartError> {
    // L-bracket: base plate + upright wall, mounting holes, filleted
    // outer verticals and chamfered hole rims.
    let base = Profile::on(Plane::XY).rect(60.0, 40.0).extrude(6.0)?;
    let wall = Profile::on(Plane::XY)
        .at(-27.0, 0.0)
        .rect(6.0, 40.0)
        .extrude(35.0)?;
    let body = base.union(&wall)?;
    let hole = Profile::on(Plane::XY)
        .at(15.0, 0.0)
        .circle(4.0)
        .extrude(6.0)?;
    let holed = body.subtract(&hole)?;
    holed.fillet(3.0, |e| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(5.0, 7.0)
    })
}

pub(crate) fn build_ribbed_plate() -> Result<Part, PartError> {
    // Plate with a stiffening rib; fillet the rib-plate junction
    // (the classic concave line fillet).
    let plate = Profile::on(Plane::XY).rect(80.0, 50.0).extrude(5.0)?;
    let rib = Profile::on(Plane::xy_at(5.0))
        .rect(80.0, 6.0)
        .extrude(12.0)?;
    let body = plate.union(&rib)?;
    body.fillet(2.0, |e| {
        e.straight()
            .direction(Dir::X)
            .chord_length_between(70.0, 90.0)
            .allow_empty()
    })
}

pub(crate) fn build_bossed_housing() -> Result<Part, PartError> {
    // Housing block with a cylindrical boss on top; fillet the boss
    // base seam (curved-seam fillet — the flagship class).
    let block = Profile::on(Plane::XY).rect(50.0, 50.0).extrude(20.0)?;
    let boss = Profile::on(Plane::xy_at(20.0)).circle(10.0).extrude(15.0)?;
    let body = block.union(&boss)?;
    body.fillet(3.0, |e| {
        e.circular()
            .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
            .allow_empty()
    })
}

pub(crate) fn build_turned_bushing() -> Result<Part, PartError> {
    // Turned part: cylinder stack with a rim fillet through the
    // periodic seam vertex (the PK-141 rim class).
    let shaft = Profile::on(Plane::XY).circle(12.0).extrude(30.0)?;
    let flange = Profile::on(Plane::XY).circle(20.0).extrude(6.0)?;
    let body = flange.union(&shaft)?;
    body.fillet(2.0, |e| {
        e.circular()
            .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
            .allow_empty()
    })
}

pub(crate) fn build_pocketed_tray() -> Result<Part, PartError> {
    // Pocketed tray with chamfered outer verticals (the bracket
    // example's class) plus a pocket-floor fillet attempt.
    let plate = Profile::on(Plane::XY).rect(60.0, 36.0).extrude(10.0)?;
    let tool = Profile::on(Plane::xy_at(4.0))
        .rect(40.0, 20.0)
        .extrude(10.0)?;
    let pocketed = plate.subtract(&tool)?;
    let chamfered = pocketed.chamfer(2.5, |e| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(8.0, 12.0)
            .allow_empty()
    })?;
    chamfered.fillet(1.5, |e| {
        e.straight()
            .direction(Dir::X)
            .chord_length_between(35.0, 45.0)
            .allow_empty()
    })
}

pub(crate) fn build_slotted_block() -> Result<Part, PartError> {
    // Through-slot interrupting long edges, then chamfer the long
    // top edges (the slot-interrupted class — recorded survivors).
    let block = Profile::on(Plane::XY).rect(70.0, 30.0).extrude(12.0)?;
    let slot = Profile::on(Plane::XY).rect(12.0, 30.0).extrude(12.0)?;
    let slotted = block.subtract(&slot)?;
    slotted.chamfer(1.5, |e| {
        e.straight()
            .direction(Dir::X)
            .chord_length_between(20.0, 80.0)
            .allow_empty()
    })
}

pub(crate) fn build_countersunk_plate() -> Result<Part, PartError> {
    // Perpendicular countersunk bore (the PK-147 class): through-bore
    // plus a coaxial drafted-frustum sink crossing the bore wall.
    let plate = Profile::on(Plane::XY).rect(60.0, 60.0).extrude(30.0)?;
    let bore = Profile::on(Plane::xy_at(-1.0)).circle(10.0).extrude(32.0)?;
    let plate = plate.subtract(&bore)?;
    let sink = Profile::on(Plane::xy_at(24.0))
        .circle(8.0)
        .extrude(8.0)?
        .draft(-30.0_f64.to_radians(), Dir::Z, |f: FaceQuery| {
            f.kind(SurfaceKind::Cylinder)
        })?;
    plate.subtract(&sink)
}

pub(crate) fn build_oblique_countersunk_plate() -> Result<Part, PartError> {
    // The PK-166 stand class: cone countersink through a tilted rest
    // plane over a blind seat.
    let plate = Profile::on(Plane::XY).rect(200.0, 200.0).extrude(140.0)?;
    let shave = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(18.0_f64.to_radians())?
        .translate(0.0, 0.0, 100.0)?;
    let plate = plate.subtract(&shave)?;
    let seat = Profile::on(Plane::xy_at(70.0)).circle(35.2).extrude(60.0)?;
    let plate = plate.subtract(&seat)?;
    let sink = Profile::on(Plane::xy_at(82.0))
        .circle(34.6)
        .extrude(36.0)?
        .draft(-8.5_f64.to_radians(), Dir::Z, |f: FaceQuery| {
            f.kind(SurfaceKind::Cylinder)
        })?;
    plate.subtract(&sink)
}

pub(crate) fn build_drafted_cup() -> Result<Part, PartError> {
    // Drafted cylinder wall — the exact-frustum class.
    Profile::on(Plane::XY).circle(30.0).extrude(40.0)?.draft(
        6.0_f64.to_radians(),
        Dir::Z,
        |f: FaceQuery| f.kind(SurfaceKind::Cylinder),
    )
}

pub(crate) fn build_cross_drilled_boss() -> Result<Part, PartError> {
    // Radial hole through a cylinder wall (cylinder-x-cylinder rims).
    let body = Profile::on(Plane::XY).circle(16.0).extrude(40.0)?;
    let drill = Profile::on(Plane::XY)
        .circle(6.0)
        .extrude(80.0)?
        .rotate_y(90.0_f64.to_radians())?
        .translate(-40.0, 0.0, 20.0)?;
    body.subtract(&drill)
}

pub(crate) fn build_coaxial_bushing() -> Result<Part, PartError> {
    // Tube: both walls full-period periodic faces.
    let outer = Profile::on(Plane::XY).circle(20.0).extrude(25.0)?;
    let bore = Profile::on(Plane::xy_at(-1.0)).circle(8.0).extrude(27.0)?;
    outer.subtract(&bore)
}

pub(crate) fn build_tilted_hole_plate() -> Result<Part, PartError> {
    // Oblique through-hole: cylinder-x-plane ellipse rims both faces.
    let plate = Profile::on(Plane::XY).rect(60.0, 60.0).extrude(15.0)?;
    let drill = Profile::on(Plane::XY)
        .circle(6.0)
        .extrude(80.0)?
        .rotate_x(25.0_f64.to_radians())?
        .translate(0.0, -15.0, -20.0)?;
    plate.subtract(&drill)
}

pub(crate) fn build_offset_bracket() -> Result<Part, PartError> {
    // L-bracket with the wall OVERLAPPING the base interior (not
    // flush) so the union stays off the coincident class; fillet the
    // outer verticals.
    let base = Profile::on(Plane::XY).rect(60.0, 40.0).extrude(8.0)?;
    let wall = Profile::on(Plane::xy_at(4.0))
        .at(-24.0, 0.0)
        .rect(8.0, 40.0)
        .extrude(35.0)?;
    let body = base.union(&wall)?;
    body.fillet(4.0, |e| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(6.0, 10.0)
            .allow_empty()
    })
}

pub(crate) fn build_rim_chamfered_disc() -> Result<Part, PartError> {
    // Large-leg rim chamfer: the conical band class.
    Profile::on(Plane::XY)
        .circle(44.0)
        .extrude(20.0)?
        .chamfer(3.0, |e| e.circular().at_height(20.0))
}

pub(crate) fn build_star_socket_plate() -> Result<Part, PartError> {
    // Hexagram socket through a plate (the star-membrane class).
    let star: Vec<(f64, f64)> = (0..12)
        .map(|i| {
            let a = std::f64::consts::PI / 6.0 * i as f64;
            let r = if i % 2 == 0 { 14.0 } else { 7.0 };
            (r * a.cos(), r * a.sin())
        })
        .collect();
    let plate = Profile::on(Plane::XY).rect(60.0, 60.0).extrude(12.0)?;
    let tool = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(14.0)?;
    plate.subtract(&tool)
}

pub(crate) fn build_mixed_corner_block() -> Result<Part, PartError> {
    // Three fillets meeting at one cube corner in a single operation
    // — the V-free DefaultOn corner class via the Part facade.
    let block = Profile::on(Plane::XY).rect(30.0, 30.0).extrude(30.0)?;
    block.fillet(3.0, |e| {
        e.straight().filter(|i| {
            let near_x = i.midpoint.x > 14.9;
            let near_y = i.midpoint.y > 14.9;
            let near_z = i.midpoint.z > 29.9;
            (near_x && near_y) || (near_x && near_z) || (near_y && near_z)
        })
    })
}

pub(crate) fn build_chamfered_peg() -> Result<Part, PartError> {
    // Cylinder peg with chamfer lead-ins on both end rims (PK-135).
    Profile::on(Plane::XY)
        .circle(10.0)
        .extrude(30.0)?
        .chamfer(1.5, |e| e.circular().at_height(30.0))?
        .chamfer(1.5, |e| e.circular().at_height(0.0))
}

pub(crate) fn build_skewed_dress_bar() -> Result<Part, PartError> {
    // Yawed+tilted boolean-cut rest plane with side chamfers — the
    // PK-168 skewed terminal-chord class.
    let bar = Profile::on(Plane::XY).rect(40.0, 60.0).extrude(96.0)?;
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(18.0_f64.to_radians())?
        .rotate_z(-12.0_f64.to_radians())?
        .translate(0.0, 0.0, 80.0)?;
    let bar = bar.subtract(&upper)?;
    let tilt = 18.0_f64.to_radians();
    let w = -12.0_f64.to_radians();
    let rest_z = move |x: f64, y: f64| 80.0 + (-x * w.sin() + y * w.cos()) * tilt.tan();
    bar.chamfer(1.0, move |e| {
        e.straight().filter(move |i| {
            i.direction.map(|d| d.y.abs() > 0.85).unwrap_or(false)
                && (i.midpoint.z - rest_z(i.midpoint.x, i.midpoint.y)).abs() < 0.5
        })
    })
}

pub(crate) fn build_counterbored_plate() -> Result<Part, PartError> {
    // Bore + larger coaxial counterbore: concentric rims on the top
    // face (the nested-hole walk class).
    let plate = Profile::on(Plane::XY).rect(50.0, 50.0).extrude(20.0)?;
    let bore = Profile::on(Plane::xy_at(-1.0)).circle(5.0).extrude(22.0)?;
    let plate = plate.subtract(&bore)?;
    let cbore = Profile::on(Plane::xy_at(12.0)).circle(9.0).extrude(9.0)?;
    plate.subtract(&cbore)
}

pub(crate) fn build_angled_rib_plate() -> Result<Part, PartError> {
    // Diagonal stiffening rib; concave junction fillet with a
    // non-axis-aligned direction filter.
    let plate = Profile::on(Plane::XY).rect(80.0, 50.0).extrude(5.0)?;
    let rib = Profile::on(Plane::xy_at(5.0))
        .rect(70.0, 6.0)
        .extrude(10.0)?
        .rotate_z(20.0_f64.to_radians())?;
    let body = plate.union(&rib)?;
    body.fillet(1.5, |e| {
        e.straight()
            .chord_length_between(40.0, 80.0)
            .filter(|i| (i.midpoint.z - 5.0).abs() < 0.1)
            .allow_empty()
    })
}

pub(crate) fn build_double_boss_plate() -> Result<Part, PartError> {
    // Two cylindrical bosses, both base seams filleted (multi-seam
    // curved fillets in one operation).
    let plate = Profile::on(Plane::XY).rect(80.0, 40.0).extrude(10.0)?;
    let boss_a = Profile::on(Plane::xy_at(10.0))
        .at(-20.0, 0.0)
        .circle(8.0)
        .extrude(12.0)?;
    let boss_b = Profile::on(Plane::xy_at(10.0))
        .at(20.0, 0.0)
        .circle(6.0)
        .extrude(15.0)?;
    let body = plate.union(&boss_a)?.union(&boss_b)?;
    body.fillet(2.0, |e| {
        e.circular()
            .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
            .allow_empty()
    })
}

pub(crate) fn build_notched_flange() -> Result<Part, PartError> {
    // Turned flange with a keyway notch through the rim, then a rim
    // chamfer over the interrupted circle (chain class).
    let flange = Profile::on(Plane::XY).circle(25.0).extrude(8.0)?;
    let key = Profile::on(Plane::xy_at(-1.0))
        .at(25.0, 0.0)
        .rect(10.0, 8.0)
        .extrude(10.0)?;
    let notched = flange.subtract(&key)?;
    notched.chamfer(1.0, |e| e.circular().at_height(8.0).allow_empty())
}
