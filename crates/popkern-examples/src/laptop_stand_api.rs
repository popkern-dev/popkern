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
    /// Base plate footprint and thickness.
    base_w: f64,
    base_d: f64,
    base_t: f64,
    /// Hub boss radius and height above the base plate.
    boss_r: f64,
    boss_h: f64,
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
            base_d: 220.0,
            base_t: 22.0,
            boss_r: 30.0,
            boss_h: 22.0,
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

/// Joint peg, v1: plain hexagram prism at fit clearance.
///
/// TEMPORARY CONTAINMENT (recorded): the designed peg adds a
/// stepped, smaller-star pilot tip at each end (2 mm lead-in for
/// blind insertion). That union — a star boss piercing a star end
/// face — is an open kernel frontier pinned as
/// `star_boss_union_gap_is_pinned` (W3, ledger PK-102 follow-up);
/// when it closes, restore the pilots from this function's history.
/// The plain prism is a complete, printable joint: the step only
/// eases insertion.
fn build_peg(p: &StandParams) -> Result<Part, PartError> {
    let body_star = hexagram(p.star_r, p.fit_clearance);
    let _ = p.pilot_clearance; // used by the pilot tips (see above)
    let _ = p.pilot_len;

    let peg = Profile::on(Plane::XY)
        .polygon(&body_star)
        .extrude(p.peg_len)?;
    p.assert_fits_plate("peg", 2.0 * p.star_r, p.peg_len);
    println!("  peg: print lying on a side face (loads shear across layers)");
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

    let hub_half = 45.0;
    let rear_half = p.base_w / 2.0;
    let rear_y = p.base_d - 8.0;
    let socket_center_y = 35.0;
    let outline: Vec<(f64, f64)> = vec![
        (-hub_half, 0.0),
        (hub_half, 0.0),
        (rear_half, rear_y),
        (-rear_half, rear_y),
    ];
    p.assert_fits_plate("base", 2.0 * rear_half, rear_y);

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
    let base = dressed.subtract(&socket_tool)?;
    println!("  base: print flat (chamfered outline, socket up)");
    Ok(base)
}

/// The angled holder, v2 (designed with the user): a full-depth
/// STACK BLOCK at the front (column seat + stacking star sockets)
/// unioned with a ROTATED ARM PLATE whose bottom face IS the 12 mm
/// arm underside; the 18° rest plane is one wedge cut across the
/// union — HIGH over the block, descending backward, so the laptop
/// slides INTO the lips. The lips are separate printed pieces that
/// slide into DOVETAIL slots across each arm tip (blind inboard end
/// seats them; the laptop pushes the blade against the slot cheeks,
/// wedging the tab down). Prints REST-FACE-DOWN: the tilted top is
/// one flat plane on the bed, the dovetail cheeks lean ~23°, the
/// blind end walls print vertical — no supports.
///
/// Every stage is pinned watertight + χ'==2 by
/// `holder_v2_full_chain_is_watertight_and_chi_clean`.
/// Recorded refinement (PK-102): a 45°-chopped catch face (print-
/// orientation insurance) is one boolean config away — the chopped
/// tool's subtract still leaks; the flat catch needs no supports in
/// the designated orientation.
fn build_holder(p: &StandParams) -> Result<Part, PartError> {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let rest_z_hub = 75.0;
    let seat_floor = 40.0;
    let socket_star = hexagram(p.star_r, 0.0);

    // Stack block: hub arc + straight rear, full depth.
    let mut block_plan: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        block_plan.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    block_plan.extend_from_slice(&[(44.0, -56.0), (-44.0, -56.0)]);
    let block = Profile::on(Plane::XY).polygon(&block_plan).extrude(95.0)?;

    // Arm plate: fanned plan with wide lip pads, authored flat and
    // rotated onto the slope; its bottom face is the arm underside.
    let arm_plan: Vec<(f64, f64)> = vec![
        (44.0, -40.0),
        (92.0, -130.0),
        (92.0, -170.0),
        (52.0, -170.0),
        (0.0, -90.0),
        (-52.0, -170.0),
        (-92.0, -170.0),
        (-92.0, -130.0),
        (-44.0, -40.0),
    ];
    let plate = Profile::on(Plane::XY)
        .polygon(&arm_plan)
        .extrude(20.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, rest_z_hub - 13.0)?;
    p.assert_fits_plate("holder", 184.0, hub_r + 170.0 * tilt.cos().max(0.9));
    let holder = block.union(&plate)?;

    // One wedge cut = the coplanar rest plane across block and arms.
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)?
        .rotate_x(tilt)?
        .translate(0.0, 0.0, rest_z_hub)?;
    let holder = holder.subtract(&upper)?;

    // Square column seat over the stack axis, flat floor.
    let seat = Profile::on(Plane::xy_at(seat_floor))
        .rect(76.0, 76.0)
        .extrude(60.0)?;
    let holder = holder.subtract(&seat)?;

    // Stacking star sockets (8 mm web between them).
    let up = Profile::on(Plane::xy_at(seat_floor - p.socket_depth))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?;
    let holder = holder.subtract(&up)?;
    let down = Profile::on(Plane::xy_at(-1.0))
        .polygon(&socket_star)
        .extrude(p.socket_depth + 1.0)?;
    let mut holder = holder.subtract(&down)?;

    // Dovetail lip slots across each arm tip: mouth 9 wide (2 mm
    // proud, vertical to 1 mm depth), cheeks leaning to 14 wide at
    // 6 deep, open at the outboard edge, blind at the inboard end.
    for side in [1.0_f64, -1.0] {
        let tool = dovetail_slot_tool(tilt, rest_z_hub, side)?;
        holder = holder.subtract(&tool)?;
    }
    println!("  holder: print REST-FACE DOWN (tilted top on the bed; no supports)");
    Ok(holder)
}

/// The dovetail slot cutter for one arm tip (`side` = ±1).
fn dovetail_slot_tool(tilt: f64, rest_z_hub: f64, side: f64) -> Result<Part, PartError> {
    let slot_y = -145.0;
    let prism = Profile::on(Plane::XY)
        .polygon(&[
            (-2.0, -4.5),
            (1.0, -4.5),
            (6.0, -7.0),
            (6.0, 7.0),
            (1.0, 4.5),
            (-2.0, 4.5),
        ])
        .extrude(41.0)?
        .rotate_y(90.0_f64.to_radians())?;
    let mut tool = prism;
    if side < 0.0 {
        tool = tool.rotate_z(std::f64::consts::PI)?;
    }
    let dz = rest_z_hub + slot_y * tilt.tan();
    tool.rotate_x(tilt)?.translate(side * 52.0, slot_y, dz)
}

/// The slide-in lip: a dovetail tab at fit clearance carrying a
/// blade that rises 16 mm proud of the rest face to stop the
/// laptop's bottom edge. Two identical lips per holder; prints on a
/// side face.
fn build_lip(p: &StandParams) -> Result<Part, PartError> {
    let c = p.fit_clearance;
    // Tab and blade share the prism axis, so the whole lip is ONE
    // profile — dovetail tab below, blade above — extruded once:
    // correct by construction, no boolean. Profile x = depth into
    // the slot (negative = blade, proud), y = width.
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
        .extrude(32.0)?
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
    let holder = gate_and_export(
        build_holder(&p).expect("holder"),
        "laptop_stand_holder",
        &out_dir,
    );
    let _lip = gate_and_export(build_lip(&p).expect("lip"), "laptop_stand_lip", &out_dir);

    // Stack preview: base, peg into the boss, column onto the peg,
    // second peg, second column. Pegs sit half-proud of each joint.
    let joint = p.peg_len - 2.0 * p.socket_depth; // peg length not in sockets (gap)
    let c1 = p.base_t + joint; // column-1 bottom z
    let c2 = c1 + p.column_h + joint;
    // The holder's own z=0 is its flat underside; it stacks on the
    // second column like another tier (its hub sits over the stack
    // axis, arms sweeping toward −y = under the laptop).
    let h1 = c2 + p.column_h + joint;
    write_assembly_preview(
        &[
            (&base, [0.0, 0.0, 0.0]),
            (&peg, [0.0, 35.0, p.base_t - p.socket_depth]),
            (&column, [0.0, 35.0, c1]),
            (&peg, [0.0, 35.0, c1 + p.column_h - p.socket_depth]),
            (&column, [0.0, 35.0, c2]),
            (&peg, [0.0, 35.0, c2 + p.column_h - p.socket_depth]),
            (&holder, [0.0, 35.0, h1]),
        ],
        &out_dir.join("laptop_stand_assembly_preview.obj"),
    )
    .expect("assembly preview");
    println!(
        "  stack preview: base + 2 columns + holder → rest plane front edge at z ≈ {:.0} mm",
        h1 + 40.0
    );
}
