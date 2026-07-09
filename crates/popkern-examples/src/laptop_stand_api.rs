//! Expandable laptop stand — the first *product* modeled with the
//! `popkern-part` API. Three of the four pieces build today: the
//! hexagram joint peg, the stacking column, and the two-leg base.
//! The angled holder needs the rotated-slab boolean (W2, in
//! progress) and ships separately.
//!
//! Joint system: loose 6-point hexagram pegs in blind star sockets
//! on every stacking face — anti-rotation by shape, no threads, no
//! supports. Sockets are cut at NOMINAL size; the peg body is
//! offset −0.2 mm per side and its pilot tips −0.8 mm, so every
//! joint gets a printed clearance fit with a stepped lead-in.
//!
//! Print orientation (strength-driven, encoded in each builder's
//! log): the peg lies on a side face (shear across layers), the
//! column stands on end (pure compression), the base prints flat.
//!
//! Every piece is gated: kernel validation at each feature step
//! (built into the API), watertight tessellation, STEP round-trip
//! with face census, and a build-plate fit assert. Outputs land in
//! `examples/output/laptop_stand/` as STEP + STL + OBJ, plus a
//! translated-copy assembly preview OBJ (meshes only, no booleans).

use popkern_part::prelude::*;

/// All product dimensions in one place, in millimetres.
struct StandParams {
    /// Build plate square side (Bambu A1 class).
    plate: f64,
    /// Hexagram outer-point radius (24 mm across points).
    star_r: f64,
    /// Per-side clearance for the peg body against a nominal socket.
    fit_clearance: f64,
    /// Extra per-side clearance on the peg pilot tips (lead-in).
    pilot_clearance: f64,
    /// Total peg length.
    peg_len: f64,
    /// Pilot tip length at each end of the peg.
    pilot_len: f64,
    /// Blind socket depth on every stacking face.
    socket_depth: f64,
    /// Column outer radius and height.
    column_r: f64,
    column_h: f64,
    /// Base plate width, thickness, and reach: the counterweight pad
    /// extends `base_front` behind the tower (+y), and two legs run
    /// `base_legs_reach` under the laptop (−y) so the rear tipping
    /// axis sits past the laptop's centre of gravity. Total depth
    /// must fit the plate.
    base_w: f64,
    base_front: f64,
    base_legs_reach: f64,
    base_t: f64,
    /// Optional ballast: two open-top trays in the front pad (the
    /// maximum-lever position from the rear tipping axis), closed by
    /// flush drop-in lids after filling with steel shot or sand.
    ballast_chambers: bool,
    /// Hub boss radius and height above the base plate.
    boss_r: f64,
    boss_h: f64,
    /// Per-side inset of the arm's swoop gusset from the plate run's
    /// side walls. 2.0 reads as a distinct rib (the default look);
    /// 0.0 makes the web fully FLUSH — a coplanar-face union,
    /// supported since the coincident-boolean campaign closed
    /// (flush_gusset_union_is_watertight_and_volume_exact).
    swoop_inset: f64,
    /// Dovetail tab relief style. `false` (default): a margined
    /// rectangular cutter that overshoots the tab by 2–3 mm per side,
    /// leaving a visible notch in the root face (the historical
    /// transversal-boolean margin rule). `true`: the cutter is the
    /// tab's own dovetail profile, every side face coplanar with an
    /// existing arm face — a crisp flush-edged relief, supported
    /// since the coincident-boolean subtract leg closed
    /// (flush_tab_relief_subtract_is_watertight_and_volume_exact).
    tab_relief_flush: bool,
    /// Arm plan length from the root face to the tip. 275 supports a
    /// ~290 mm laptop with the original CG margin (~35 mm); the
    /// SPLIT arm (root + tip with the sliding-dovetail lap) is the
    /// printable form — a one-piece arm at this length only just
    /// fits a 240 plate diagonally.
    arm_len: f64,
    /// Half-width step between the arm's root body and its plate
    /// run, per side. 2.0 (default): the root is 44 mm wide against
    /// the 40 mm run — a visible shoulder at the dovetail joint.
    /// 0.0: uniform 40 mm width from root face to tip, the walls one
    /// continuous plane (pairs with swoop_inset = 0 for the fully
    /// flush arm).
    root_shoulder: f64,
    /// Central seat recess profile on the block. `false` (default):
    /// the 76 × 76 square. `true`: a Ø76 cylindrical recess —
    /// supported since the star-membrane campaign closed
    /// (round_seat_star_socket_is_watertight_and_void).
    seat_round: bool,
}

impl StandParams {
    fn default_stand() -> Self {
        StandParams {
            plate: 240.0,
            star_r: 12.0,
            fit_clearance: 0.2,
            pilot_clearance: 0.8,
            peg_len: 30.0,
            pilot_len: 6.0,
            socket_depth: 16.0,
            column_r: 35.0,
            column_h: 60.0,
            base_w: 200.0,
            base_front: 70.0,
            base_legs_reach: 162.0,
            base_t: 22.0,
            ballast_chambers: false,
            boss_r: 30.0,
            boss_h: 22.0,
            swoop_inset: 2.0,
            arm_len: 275.0,
            tab_relief_flush: false,
            root_shoulder: 2.0,
            seat_round: false,
        }
    }

    /// Web of material left between the column's two opposed sockets.
    fn column_web(&self) -> f64 {
        self.column_h - 2.0 * self.socket_depth
    }

    fn assert_fits_plate(&self, name: &str, w: f64, d: f64) {
        assert!(
            w <= self.plate && d <= self.plate,
            "{name} footprint {w:.0}x{d:.0} exceeds {:.0} plate",
            self.plate
        );
    }
}

/// Hexagram (6-point star): 12 vertices alternating between the
/// outer-point radius R and the inner radius R/√3, optionally shrunk
/// by offsetting every edge inward along its outward normal —
/// per-FACE clearance, which radial scaling cannot give on a star.
fn hexagram(r_outer: f64, edge_offset: f64) -> Vec<(f64, f64)> {
    let raw: Vec<(f64, f64)> = (0..12)
        .map(|k| {
            let a = std::f64::consts::TAU * k as f64 / 12.0;
            let r = if k % 2 == 0 {
                r_outer
            } else {
                r_outer / 3.0_f64.sqrt()
            };
            (r * a.cos(), r * a.sin())
        })
        .collect();
    if edge_offset == 0.0 {
        return raw;
    }
    // Shift each edge line inward by `edge_offset`, then intersect
    // consecutive shifted lines to recover the vertices.
    let n = raw.len();
    let shifted: Vec<((f64, f64), (f64, f64))> = (0..n)
        .map(|i| {
            let (x0, y0) = raw[i];
            let (x1, y1) = raw[(i + 1) % n];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let len = (dx * dx + dy * dy).sqrt();
            // Outward normal of a CCW polygon edge is (dy, -dx).
            let (nx, ny) = (dy / len, -dx / len);
            ((x0 - nx * edge_offset, y0 - ny * edge_offset), (dx, dy))
        })
        .collect();
    (0..n)
        .map(|i| {
            let (p_prev, d_prev) = shifted[(i + n - 1) % n];
            let (p_cur, d_cur) = shifted[i];
            // Intersect p_prev + s·d_prev with p_cur + t·d_cur.
            let det = d_prev.0 * d_cur.1 - d_prev.1 * d_cur.0;
            let s = ((p_cur.0 - p_prev.0) * d_cur.1 - (p_cur.1 - p_prev.1) * d_cur.0) / det;
            (p_prev.0 + s * d_prev.0, p_prev.1 + s * d_prev.1)
        })
        .collect()
}

/// Joint peg, plain hexagram prism at fit clearance.
///
/// TEMPORARY CONTAINMENT (recorded): the designed peg adds stepped
/// pilot tips at BOTH ends. The single star-boss union is closed
/// (W3, `star_boss_union_is_watertight`), but the SECOND pilot union
/// onto the already-stepped body leaks (pinned
/// `peg_double_pilot_gap_is_pinned`) — the second-boolean-on-a-
/// unioned-star-body class. The plain prism remains a complete
/// printable joint.
fn build_peg(p: &StandParams) -> Result<Part, PartError> {
    let body_star = hexagram(p.star_r, p.fit_clearance);
    let _ = p.pilot_clearance; // returns with the double-pilot pin
    let _ = p.pilot_len;

    let peg = Profile::on(Plane::XY)
        .polygon(&body_star)
        .extrude(p.peg_len)?;
    p.assert_fits_plate("peg", 2.0 * p.star_r, p.peg_len);
    println!("  peg: print lying on a side face");
    Ok(peg)
}

/// Stacking column: plain cylinder with a blind NOMINAL star socket
/// in each end. Both sockets accept the same peg; the web between
/// them keeps the column solid.
fn build_column(p: &StandParams) -> Result<Part, PartError> {
    let socket_star = hexagram(p.star_r, 0.0);
    assert!(
        p.column_web() >= 20.0,
        "column web {:.0} too thin",
        p.column_web()
    );

    let blank = Profile::on(Plane::XY)
        .circle(p.column_r)
        .extrude(p.column_h)?;
    // Tools run 1 mm proud of each end face for clean transversal
    // crossings.
    let bottom_tool = Profile::on(Plane::xy_at(-1.0))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?;
    let top_tool = Profile::on(Plane::xy_at(p.column_h - p.socket_depth))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?;

    let column = blank.subtract(&bottom_tool)?.subtract(&top_tool)?;
    p.assert_fits_plate("column", 2.0 * p.column_r, 2.0 * p.column_r);
    println!("  column: print standing on end (loads pure compression)");
    Ok(column)
}

/// Base, v1: a convex trapezoid plate — 90 mm front edge under the
/// holder's tower widening to the full 200 mm at the rear — with the
/// star socket cut straight into the hub.
///
/// TEMPORARY CONTAINMENT (recorded, two layers):
/// 1. The designed base is a 10 mm two-leg deck with a Ø60×22 socket
///    boss. The boss union on a reflex plate is pinned as
///    `cylinder_boss_on_reflex_plate_union_gap_is_pinned`; v1 uses a
///    22 mm plate with the socket cut directly.
/// 2. The two-leg outline itself is deferred: an interior socket
///    ring on a REFLEX face fails at the split stage ("trim curve
///    end point does not match any vertex"), and carving the leg
///    V-cutout from a convex plate is a through + side-overlap
///    boolean (both open frontiers). The trapezoid gives the same
///    stance and stability with proven operations; the sculpted
///    legs return with the W2 campaign.
fn build_base(p: &StandParams) -> Result<Part, PartError> {
    let socket_star = hexagram(p.star_r, 0.0);
    let _ = (p.boss_r, p.boss_h); // return with the boss (see above)

    let half_w = p.base_w / 2.0;
    let front = p.base_front; // counterweight pad, behind the tower (+y)
    let legs = p.base_legs_reach; // under the laptop (−y)
    let pad_edge = -30.0; // pad's laptop-side edge; legs continue from here
    let socket_center_y = 35.0;
    // Two 44 mm legs under the arms (arm centres at x = ±28); the
    // notch between them keeps the outline printable and saves
    // plastic while the rear tip line (the tipping axis) sits past
    // the laptop CG.
    let outline: Vec<(f64, f64)> = vec![
        (-half_w, front),
        (-half_w, pad_edge),
        (-54.0, pad_edge),
        (-54.0, -legs),
        (-10.0, -legs),
        (-10.0, pad_edge),
        (10.0, pad_edge),
        (10.0, -legs),
        (54.0, -legs),
        (54.0, pad_edge),
        (half_w, pad_edge),
        (half_w, front),
    ];
    p.assert_fits_plate("base", 2.0 * half_w, front + legs);

    let plate = Profile::on(Plane::XY).polygon(&outline).extrude(p.base_t)?;

    // Dress the outer vertical edges BEFORE cutting the socket
    // (measured rule: pockets-after-dressing).
    let verticals = |e: EdgeQuery| e.straight().direction(Dir::Z);
    let dressed = plate.chamfer(2.0, verticals)?;

    // Full-depth star socket straight into the hub (blind: 6 mm
    // floor at 22 mm plate).
    assert!(p.base_t - p.socket_depth >= 4.0, "socket floor too thin");
    let socket_tool = Profile::on(Plane::xy_at(p.base_t - p.socket_depth))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?
        .translate(0.0, socket_center_y, 0.0)?;
    let mut base = dressed.subtract(&socket_tool)?;

    if p.ballast_chambers {
        // Two open-top trays flanking the socket in the front pad —
        // the farthest usable position from the rear tipping axis,
        // and clear of the socket's 6 mm floor. Each tray: 14 mm
        // deep over a 4 mm floor, with a 4 mm-deep rebate ledge so
        // the lid drops in flush with the base top.
        // Feature order is load-bearing: the shallow rebate first,
        // then the tray into its floor — the shipping block's proven
        // pocket-into-pocket-floor order (the reverse is a
        // counterbore over an existing hole, a pinned class).
        for side in [1.0_f64, -1.0] {
            let rebate = Profile::on(Plane::xy_at(p.base_t - 4.0))
                .rect(51.0, 66.0)
                .at(side * 62.5, 30.0)
                .extrude(8.0)?;
            base = base.subtract(&rebate)?;
            let tray = Profile::on(Plane::xy_at(4.0))
                .rect(45.0, 60.0)
                .at(side * 62.5, 30.0)
                .extrude(p.base_t - 6.0)?;
            base = base.subtract(&tray)?;
        }
        let shot_cc = 2.0 * (45.0 * 60.0 * (p.base_t - 8.0)) / 1000.0;
        println!(
            "  base: ballast trays hold ≈ {shot_cc:.0} cm³ → ~{:.0} g of steel shot",
            shot_cc * 4.6
        );
    }
    println!("  base: print flat (chamfered outline, socket up)");
    Ok(base)
}

/// Flush drop-in lid for the ballast trays (print two): sits on the
/// tray's rebate ledge at fit clearance, top flush with the base.
fn build_ballast_lid(p: &StandParams) -> Result<Part, PartError> {
    let c = p.fit_clearance;
    let lid = Profile::on(Plane::XY)
        .rect(51.0 - 2.0 * c, 66.0 - 2.0 * c)
        .extrude(4.0)?;
    println!("  ballast lid: print two, glue or press onto the rebate after filling");
    Ok(lid)
}

/// Holder v3 BLOCK — the stacking piece. Modular architecture
/// (designed with the user): the block carries the column seat, both
/// stacking star sockets, and two VERTICAL dovetail sockets in its
/// rear face where the arms hang; the arms are separate identical
/// pieces so their length is not capped by the build plate, which is
/// what restores stacking clearance (the laptop tops out short of
/// the seat). Chain pinned by
/// `holder_v3_block_chain_is_watertight_and_chi_clean`.
/// The front is a SMOOTH ARC — a true Cylinder wall from the
/// mixed-path sketch API. Booleans pinned χ-clean
/// (arc_block_arm_sockets_are_watertight); the rational-curve STEP
/// export uses the standard Part-21 complex form (OCC-verified).
fn build_block(p: &StandParams) -> Result<Part, PartError> {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let rest_z_hub = 100.0;
    let seat_floor = 70.0;
    let socket_star = hexagram(p.star_r, 0.0);

    let _ = hub_r;
    p.assert_fits_plate("block", 116.0, 114.0);

    // SMOOTH arc front — one arc segment in the sketch, extruded to
    // a true Cylinder wall.
    let block = Profile::on(Plane::XY)
        .path(58.0, 0.0)
        .arc_to((0.0, 58.0), (-58.0, 0.0))
        .line_to(-44.0, -56.0)
        .line_to(44.0, -56.0)
        .close()
        .extrude(125.0)?;
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, rest_z_hub)?;
    let block = block.subtract(&upper)?;

    // Seat recess (square or round), then the star socket into its
    // floor. The round seat is supported since the star-membrane
    // campaign closed (round_seat_star_socket_is_watertight_and_void:
    // periodic-chart containment + loop-winding closure).
    let seat = if p.seat_round {
        Profile::on(Plane::xy_at(seat_floor))
            .circle(38.0)
            .extrude(60.0)?
    } else {
        Profile::on(Plane::xy_at(seat_floor))
            .rect(76.0, 76.0)
            .extrude(60.0)?
    };
    let block = block.subtract(&seat)?;
    let up = Profile::on(Plane::xy_at(seat_floor - p.socket_depth))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?;
    let block = block.subtract(&up)?;
    let down = Profile::on(Plane::xy_at(-1.0))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?;
    let mut block = block.subtract(&down)?;

    // Vertical dovetail arm sockets: the arm's bending moment
    // resolves as a horizontal couple in the dovetail; insertion is
    // from above, and the laptop's weight locks it.
    for side in [1.0_f64, -1.0] {
        let tool = Profile::on(Plane::xy_at(20.0))
            .polygon(&[(-6.0, -7.0), (6.0, -7.0), (9.0, 4.0), (-9.0, 4.0)])
            .at(side * 28.0, -50.0)
            .extrude(110.0)?;
        block = block.subtract(&tool)?;
    }
    println!("  block: print REST-FACE DOWN");
    Ok(block)
}

/// Holder v3 ARM — one identical piece, printed twice. One plan
/// outline (root body + dovetail tab + plate run) extruded, the rest
/// plane cut across it, the plate run thinned to 13 mm by a
/// tool-shaped lower cutter, a true-arc swoop gusset for rigidity,
/// a THROUGH dovetail lip slot near the tip (through ⇒ genus 1 ⇒
/// χ' = 0 by design), and the tab relieved below the socket floor.
/// Chain pinned by `holder_v3_arm_chain_is_watertight_and_chi_clean`.
fn build_arm(p: &StandParams) -> Result<Part, PartError> {
    build_arm_shape(p, p.arm_len, true)
}

/// The arm chain with a parametric plan tip: the SPLIT pieces reuse
/// it (the root half ends at the joint shoulder and carries no
/// slot).
fn build_arm_shape(p: &StandParams, plan_tip: f64, with_slot: bool) -> Result<Part, PartError> {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8; // block rest plane at its rear face
    let root_half = 20.0 + p.root_shoulder;
    let plan: Vec<(f64, f64)> = if p.root_shoulder > 0.0 {
        vec![
            (root_half, 0.0),
            (5.8, 0.0),
            (8.8, 9.8),
            (-8.8, 9.8),
            (-5.8, 0.0),
            (-root_half, 0.0),
            (-root_half, -24.0),
            (-20.0, -24.0),
            (-20.0, -plan_tip),
            (20.0, -plan_tip),
            (20.0, -24.0),
            (root_half, -24.0),
        ]
    } else {
        // Uniform width: the shoulder vertices vanish and each side
        // wall is ONE plane from the root face to the tip.
        vec![
            (20.0, 0.0),
            (5.8, 0.0),
            (8.8, 9.8),
            (-8.8, 9.8),
            (-5.8, 0.0),
            (-20.0, 0.0),
            (-20.0, -plan_tip),
            (20.0, -plan_tip),
        ]
    };
    // The rest plane descends from z = 81.8 at the root; past
    // y ≈ −252 the plate sits entirely below z = 0, so the run's
    // extrude floor follows the plate depth at the tip while the
    // ROOT body keeps its flat z = 0 underside (cut back below).
    let plate_under = |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let z_lo = (plate_under(-plan_tip) - 2.0).min(-2.0);
    let mut arm = Profile::on(Plane::xy_at(z_lo))
        .polygon(&plan)
        .extrude(96.0 - z_lo)?;
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, plane_at_root)?;
    arm = arm.subtract(&upper)?;
    let root_floor = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(60.0, 46.0)
        .at(0.0, -8.0)
        .extrude(5.0 - z_lo)?;
    arm = arm.subtract(&root_floor)?;
    let lower_box = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(600.0, p.arm_len - 16.0)
        .at(0.0, -(26.0 + p.arm_len + 10.0) / 2.0)
        .extrude(120.0)?;
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())?;
    let tool = lower_box.subtract(&p2)?;
    arm = arm.subtract(&tool)?;
    // TRUE-ARC support swoop. The web's TOP edge runs PARALLEL to
    // the plate underside, embedded 1 mm, so it bears along the full
    // run; the free edge is ONE concave arc from the plate down into
    // the root collar (arc_gusset_union_is_watertight). Profile
    // plane maps (px, py) -> (y = py, z = -px) after the 90-degree
    // Y rotation.
    let (y_root, y_tip) = (-20.0, -72.0);
    let half_run = 20.0 - p.swoop_inset;
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(2.0 * half_run)?
        .rotate_y(90.0_f64.to_radians())?
        .translate(-half_run, 0.0, 0.0)?;
    arm = arm.union(&gusset)?;
    if with_slot {
        let slot = Profile::on(Plane::XY)
            .polygon(&[
                (-2.0, -4.5),
                (1.0, -4.5),
                (6.0, -7.0),
                (6.0, 7.0),
                (1.0, 4.5),
                (-2.0, 4.5),
            ])
            .extrude(44.0)?
            .rotate_y(90.0_f64.to_radians())?
            .rotate_x(tilt)?
            .translate(
                -22.0,
                -(p.arm_len - 16.0),
                plane_at_root - (p.arm_len - 16.0) * tilt.tan(),
            )?;
        arm = arm.subtract(&slot)?;
    }
    // Relieve the tab below the block's socket floor so the arm
    // HANGS on the floor while the root face bears on the wall.
    // Flush style: the cutter IS the tab's dovetail profile, so all
    // four side faces of the cut are coplanar with existing arm
    // faces (two dovetail flanks, the tab front, and the root face
    // plane) — no margin notch. Margined style keeps the historical
    // 2–3 mm transversal overshoot.
    let relief = if p.tab_relief_flush {
        // Flank lines extended 1 mm past the tab front (to y = 10.8)
        // so the cutter's front face exits through air instead of
        // landing exactly on the tab front — one less full-face
        // coincidence, with no change to the visible flush edges.
        Profile::on(Plane::xy_at(-1.0))
            .polygon(&[
                (5.8, 0.0),
                (5.8 + 3.0 * (10.8 / 9.8), 10.8),
                (-5.8 - 3.0 * (10.8 / 9.8), 10.8),
                (-5.8, 0.0),
            ])
            .extrude(21.0)?
    } else {
        Profile::on(Plane::xy_at(-5.0))
            .rect(24.0, 14.0)
            .at(0.0, 5.0)
            .extrude(25.0)?
    };
    let arm = arm.subtract(&relief)?;
    // 44 × 241 footprint prints diagonally on a 240 plate.
    println!("  arm: print two, REST-FACE DOWN, placed diagonally");
    Ok(arm)
}

/// Sliding-dovetail lap joint constants for the SPLIT arm. The seam
/// sits past the swoop where the section is the plain 40 × 13 plate
/// and the bending moment is a fraction of the root value. The top
/// shoulder (visible seam) is the single axial stop; the rail's nose
/// is shortened so it never bottoms first.
const ARM_JOINT_MOUTH_Y: f64 = 120.0; // tunnel mouth / underside seam
const ARM_JOINT_SHOULDER_Y: f64 = 170.0; // top shoulder = the stop
const ARM_RAIL_HALF_BASE: f64 = 11.0; // dovetail half-width at the underside
const ARM_RAIL_HALF_NECK: f64 = 7.0; // half-width at the neck (top)
const ARM_RAIL_HEIGHT: f64 = 7.0; // of the 13 mm plate (6 mm ceiling)

/// The dovetail prism used by both joint halves, aligned to the
/// TILTED plate: profile drawn in XY (x across the arm, y along the
/// plate normal, underside at y = 0), extruded along +z, then
/// rotated so the axis runs down the plate (rotate_x by 90° + tilt)
/// and translated so the profile plane sits at plan station
/// `y_start` on the plate underside. `inset` shrinks the flanks and
/// ceiling by a per-side clearance (0 = the tunnel's nominal shape);
/// `below` extends the flank lines past the underside so the cutter
/// crosses it transversally (a flank ENDING exactly on the underside
/// would put a cutter edge in the face — the pinned corner-contact
/// class).
fn arm_joint_prism(
    tilt: f64,
    y_start: f64,
    length_along_plate: f64,
    inset: f64,
    below: f64,
    plate_under: &dyn Fn(f64) -> f64,
) -> Result<Part, PartError> {
    // Flank line: x(y) = ±(base − slope·y), slope = (base − neck)/h.
    let slope = (ARM_RAIL_HALF_BASE - ARM_RAIL_HALF_NECK) / ARM_RAIL_HEIGHT;
    // Perpendicular inset of a flank shifts x by inset·√(1+slope²).
    let flank_shift = inset * (1.0 + slope * slope).sqrt();
    let x_at = |y: f64| ARM_RAIL_HALF_BASE - flank_shift - slope * y;
    let y_top = ARM_RAIL_HEIGHT - inset;
    let y_bot = -below;
    let profile = vec![
        (x_at(y_bot), y_bot),
        (x_at(y_top), y_top),
        (-x_at(y_top), y_top),
        (-x_at(y_bot), y_bot),
    ];
    Profile::on(Plane::XY)
        .polygon(&profile)
        .extrude(length_along_plate)?
        .rotate_x(std::f64::consts::FRAC_PI_2 + tilt)?
        .translate(0.0, -y_start, plate_under(-y_start))
}

/// ROOT half of the split arm: the arm chain with its plan ending at
/// the joint shoulder (no slot — that lives on the tip), plus the
/// dovetail TUNNEL opening in the underside over the lap zone.
/// Prints rest-face down like the full arm.
fn build_arm_root(p: &StandParams) -> Result<Part, PartError> {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let plate_under = move |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let root = build_arm_shape(p, ARM_JOINT_SHOULDER_Y, false)?;
    let lap_plate = (ARM_JOINT_SHOULDER_Y - ARM_JOINT_MOUTH_Y) / tilt.cos();
    let tunnel = arm_joint_prism(
        tilt,
        ARM_JOINT_MOUTH_Y,
        lap_plate + 5.0 / tilt.cos(),
        0.0,
        2.0,
        &plate_under,
    )?;
    let root = root.subtract(&tunnel)?;
    println!("  arm root: print REST-FACE DOWN; tunnel opens at the underside");
    Ok(root)
}

/// TIP half of the split arm: its own simple chain — a rectangular
/// plan from the tunnel mouth to the arm tip, the rest-plane cut,
/// the plate-thinning tool, the lip slot — then the lap zone is
/// carved down to the joint RAIL (dovetail at fit clearance, nose
/// 0.4 mm short so the top shoulder is the single stop). Slides
/// into the root piece along the plate and glues.
fn build_arm_tip(p: &StandParams) -> Result<Part, PartError> {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let plate_under = move |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let z_lo = (plate_under(-p.arm_len) - 2.0).min(-2.0);
    // Plain plate segment from the shoulder to the tip: rest cut,
    // thinning tool, lip slot — all canonical. The RAIL is unioned
    // LAST (embedded 1 mm past the shoulder; its floor is coplanar
    // with the finished underside — the closed coincident-union
    // class). Building it first and carving the lap zone instead
    // needs record-network T-junction splitting at the shoulder
    // cross-section — the pinned junction-network frontier.
    let mut tip = Profile::on(Plane::xy_at(z_lo))
        .polygon(&[
            (20.0, -ARM_JOINT_SHOULDER_Y),
            (20.0, -p.arm_len),
            (-20.0, -p.arm_len),
            (-20.0, -ARM_JOINT_SHOULDER_Y),
        ])
        .extrude(96.0 - z_lo)?;
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, plane_at_root)?;
    tip = tip.subtract(&upper)?;
    let lower_box = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(600.0, p.arm_len - 16.0)
        .at(0.0, -(26.0 + p.arm_len + 10.0) / 2.0)
        .extrude(120.0)?;
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())?;
    let tool = lower_box.subtract(&p2)?;
    tip = tip.subtract(&tool)?;
    let slot = Profile::on(Plane::XY)
        .polygon(&[
            (-2.0, -4.5),
            (1.0, -4.5),
            (6.0, -7.0),
            (6.0, 7.0),
            (1.0, 4.5),
            (-2.0, 4.5),
        ])
        .extrude(44.0)?
        .rotate_y(90.0_f64.to_radians())?
        .rotate_x(tilt)?
        .translate(
            -22.0,
            -(p.arm_len - 16.0),
            plane_at_root - (p.arm_len - 16.0) * tilt.tan(),
        )?;
    tip = tip.subtract(&slot)?;
    let rail = arm_joint_prism(
        tilt,
        ARM_JOINT_MOUTH_Y + 0.4,
        (ARM_JOINT_SHOULDER_Y - ARM_JOINT_MOUTH_Y - 0.4 + 1.0) / tilt.cos(),
        p.fit_clearance,
        0.0,
        &plate_under,
    )?;
    let tip = tip.union(&rail)?;
    println!("  arm tip: print REST-FACE DOWN; rail slides into the root and glues");
    Ok(tip)
}

/// The slide-in lip: dovetail tab at fit clearance + blade, one
/// 10-gon profile × one extrude (no booleans). The slot is a
/// through-slot so identical arms serve both sides; the lip inserts
/// from the outboard edge.
fn build_lip(p: &StandParams) -> Result<Part, PartError> {
    let c = p.fit_clearance;
    let lip = Profile::on(Plane::XY)
        .polygon(&[
            (0.0, -(4.5 - c)),
            (1.0, -(4.5 - c)),
            (5.6, -(6.6 - c)),
            (5.6, 6.6 - c),
            (1.0, 4.5 - c),
            (0.0, 4.5 - c),
            (0.0, 4.0),
            (-17.0, 4.0),
            (-17.0, -4.0),
            (0.0, -4.0),
        ])
        .extrude(40.0)?
        .rotate_y(90.0_f64.to_radians())?;
    println!("  lip: print on a side face; two per holder");
    Ok(lip)
}

/// Write a merged multi-piece OBJ (translated meshed copies — a
/// visual stack preview, not a boolean assembly).
fn write_assembly_preview(
    pieces: &[(&Part, [f64; 3])],
    path: &std::path::Path,
) -> std::io::Result<()> {
    use std::io::Write as _;
    let mut out = String::new();
    let mut base_index = 1usize;
    let options = popkern_mesh::TessellationOptions::default();
    for (part, offset) in pieces {
        let mesh = part.mesh(&options).expect("preview mesh");
        for v in &mesh.vertices {
            out.push_str(&format!(
                "v {} {} {}\n",
                v.position.x + offset[0],
                v.position.y + offset[1],
                v.position.z + offset[2]
            ));
        }
        for t in &mesh.triangles {
            out.push_str(&format!(
                "f {} {} {}\n",
                t[0] + base_index,
                t[1] + base_index,
                t[2] + base_index
            ));
        }
        base_index += mesh.vertices.len();
    }
    let mut f = std::fs::File::create(path)?;
    f.write_all(out.as_bytes())
}

fn gate_and_export(part: Part, name: &str, out_dir: &std::path::Path) -> Part {
    // An EMPTY body meshes trivially watertight — a silently swallowed
    // boolean must not pass the gate (measured: an arm tip whose run
    // had wedged out to nothing exported as a valid 0-face part).
    assert!(part.face_count() > 0, "{name}: empty body");
    let part = part.assert_watertight().expect("watertight gate");
    let step = out_dir.join(format!("{name}.step"));
    part.export_step(&step).expect("STEP export");
    part.export_stl(out_dir.join(format!("{name}.stl")))
        .expect("STL export");
    part.export_obj(out_dir.join(format!("{name}.obj")))
        .expect("OBJ export");
    let part = part.assert_step_roundtrip(&step).expect("round-trip gate");
    println!(
        "  {name}: watertight ✓ round-trip ✓ ({} faces) → {}",
        part.face_count(),
        step.display()
    );
    part
}

/// Variant pieces are exploratory: gate them like production pieces
/// but report failures as frontier signatures instead of aborting
/// the standard export run.
fn gate_and_export_variant(part: Result<Part, PartError>, name: &str, out_dir: &std::path::Path) {
    let built = match part {
        Ok(p) => p,
        Err(e) => {
            println!("  {name}: BUILD FAILED — {e:?}");
            return;
        }
    };
    match built.assert_watertight() {
        Ok(p) => {
            let step = out_dir.join(format!("{name}.step"));
            if let Err(e) = p.export_step(&step) {
                println!("  {name}: STEP export failed — {e:?}");
                return;
            }
            match p.assert_step_roundtrip(&step) {
                Ok(p) => println!(
                    "  {name}: watertight ✓ round-trip ✓ ({} faces) → {}",
                    p.face_count(),
                    step.display()
                ),
                Err(e) => println!("  {name}: round-trip FAILED — {e:?}"),
            }
        }
        Err(e) => println!("  {name}: NOT WATERTIGHT — {e:?}"),
    }
}

fn main() {
    let out_dir =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output/laptop_stand");
    std::fs::create_dir_all(&out_dir).expect("out dir");
    let p = StandParams::default_stand();

    println!("laptop stand pieces (plate {:.0} mm):", p.plate);
    let peg = gate_and_export(build_peg(&p).expect("peg"), "laptop_stand_peg", &out_dir);
    let column = gate_and_export(
        build_column(&p).expect("column"),
        "laptop_stand_column",
        &out_dir,
    );
    let base = gate_and_export(build_base(&p).expect("base"), "laptop_stand_base", &out_dir);
    let block = gate_and_export(
        build_block(&p).expect("block"),
        "laptop_stand_block",
        &out_dir,
    );
    let arm_root = gate_and_export(
        build_arm_root(&p).expect("arm root"),
        "laptop_stand_arm_root",
        &out_dir,
    );
    let arm_tip = gate_and_export(
        build_arm_tip(&p).expect("arm tip"),
        "laptop_stand_arm_tip",
        &out_dir,
    );
    let lip = gate_and_export(build_lip(&p).expect("lip"), "laptop_stand_lip", &out_dir);

    // Requested variants — exploratory, reported honestly.
    println!("variants:");
    let flush_swoop = StandParams {
        swoop_inset: 0.0,
        ..StandParams::default_stand()
    };
    gate_and_export_variant(
        build_arm(&flush_swoop),
        "laptop_stand_arm_flush_swoop",
        &out_dir,
    );
    let flush_arm = StandParams {
        swoop_inset: 0.0,
        tab_relief_flush: true,
        root_shoulder: 0.0,
        ..StandParams::default_stand()
    };
    let round_seat = StandParams {
        seat_round: true,
        ..StandParams::default_stand()
    };
    gate_and_export_variant(
        build_block(&round_seat),
        "laptop_stand_block_round_seat",
        &out_dir,
    );
    let ballast = StandParams {
        ballast_chambers: true,
        ..StandParams::default_stand()
    };
    gate_and_export_variant(build_base(&ballast), "laptop_stand_base_ballast", &out_dir);
    gate_and_export_variant(
        build_ballast_lid(&ballast),
        "laptop_stand_ballast_lid",
        &out_dir,
    );
    gate_and_export_variant(build_arm(&flush_arm), "laptop_stand_arm_flush", &out_dir);

    // Stack preview: base, peg into the boss, column onto the peg,
    // second peg, second column. Pegs sit half-proud of each joint.
    let joint = p.peg_len - 2.0 * p.socket_depth; // peg length not in sockets (gap)
    let c1 = p.base_t + joint; // column-1 bottom z
    let c2 = c1 + p.column_h + joint;
    // The holder's own z=0 is its flat underside; it stacks on the
    // second column like another tier (its hub sits over the stack
    // axis, arms sweeping toward −y = under the laptop).
    let h1 = c2 + p.column_h + joint;
    let tilt = 18.0_f64.to_radians();
    write_assembly_preview(
        &[
            (&base, [0.0, 0.0, 0.0]),
            (&peg, [0.0, 35.0, p.base_t - p.socket_depth]),
            (&column, [0.0, 35.0, c1]),
            (&peg, [0.0, 35.0, c1 + p.column_h - p.socket_depth]),
            (&column, [0.0, 35.0, c2]),
            (&peg, [0.0, 35.0, c2 + p.column_h - p.socket_depth]),
            (&block, [0.0, 35.0, h1]),
            (&arm_root, [28.0, 35.0 - 56.0, h1]),
            (&arm_tip, [28.0, 35.0 - 56.0, h1]),
            (&arm_root, [-28.0, 35.0 - 56.0, h1]),
            (&arm_tip, [-28.0, 35.0 - 56.0, h1]),
            (
                &lip,
                [
                    28.0 - 20.0,
                    35.0 - 56.0 - (p.arm_len - 16.0),
                    h1 + 81.8 - (p.arm_len - 16.0) * tilt.tan(),
                ],
            ),
            (
                &lip,
                [
                    -28.0 - 20.0,
                    35.0 - 56.0 - (p.arm_len - 16.0),
                    h1 + 81.8 - (p.arm_len - 16.0) * tilt.tan(),
                ],
            ),
        ],
        &out_dir.join("laptop_stand_assembly_preview.obj"),
    )
    .expect("assembly preview");
    println!(
        "  stack preview: base + 2 columns + block + arms + lips → seat mouth at z ≈ {:.0} mm",
        h1 + 100.0
    );
    // Stability FYI: laptop CG plan position vs the base's rear
    // reach. The base v1 trapezoid does not yet implement the
    // splayed rear legs from the design doc; a negative margin here
    // means the base needs them before shipping long arms.
    let lip_y = 35.0 - 56.0 - (p.arm_len - 16.0);
    let laptop_depth = 290.0_f64;
    let cg_y = lip_y + (laptop_depth / 2.0) * tilt.cos();
    let base_rear = -p.base_legs_reach;
    println!(
        "  stability: {laptop_depth:.0} mm laptop CG at plan y ≈ {cg_y:.0} vs base rear {base_rear:.0} → margin {:.0} mm",
        cg_y - base_rear
    );
}
