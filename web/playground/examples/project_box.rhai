// A printable project box: the base is shelled to a uniform wall and
// grows four screw bosses from a loop; the lid nests with a lip, is
// vented by another loop, and its screw holes get exact-cone
// countersinks. Two bodies from one script; a configuration turns
// the vents off.
let m = model();
m.param("wall", 0.06, "m");

// ---- Base ----
let base = m.profile(xy()).rect(1.6, 1.0).sketch("base_outline").extrude(blind(0.6));
base = base.fillet(edges().straight().direction("z"), 0.15, #{ label: "base_corners" });
base = base.shell(faces().facing("z").at_height(0.6), p("wall"), #{ label: "cavity" });

// Screw bosses at the four inner corners, pilot holes drilled in the
// same loop. Positions are computed, not typed four times.
for corner in [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]] {
    let cx = corner[0] * 0.61;
    let cy = corner[1] * 0.31;
    let boss = m.tool(|m| m.profile(xy_at(0.06)).circle(0.1).at(cx, cy)
        .sketch("boss_" + corner[0] + "_" + corner[1])
        .extrude(blind(0.5)));
    base = base.union(boss);
    let pilot = m.tool(|m| m.profile(xy_at(0.25)).circle(0.04).at(cx, cy)
        .sketch("pilot_" + corner[0] + "_" + corner[1])
        .extrude(blind(0.4)));
    base = base.subtract(pilot);
}

// ---- Lid ----
let lid = m.profile(xy_at(0.72)).rect(1.6, 1.0).sketch("lid_outline").extrude(blind(0.08));
lid = lid.fillet(edges().straight().direction("z"), 0.15, #{ label: "lid_corners" });
// The nesting lip is four bars (0.02 clearance a side against the
// opening) so the lid's centre stays single-walled and the vents cut
// plain plate.
for bar in [[0.0, -0.39, 1.44, 0.06], [0.0, 0.39, 1.44, 0.06],
            [-0.69, 0.0, 0.06, 0.66], [0.69, 0.0, 0.06, 0.66]] {
    let rib = m.tool(|m| m.profile(xy_at(0.66)).rect(bar[2], bar[3]).at(bar[0], bar[1])
        .sketch("lip_bar_" + bar[0] + "_" + bar[1])
        .extrude(blind(0.06)));
    lid = lid.union(rib);
}

// Screw holes over the bosses, countersunk in one chamfer.
for corner in [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]] {
    let hole = m.tool(|m| m.profile(xy_at(0.64)).circle(0.04).at(corner[0] * 0.61, corner[1] * 0.31)
        .sketch("hole_" + corner[0] + "_" + corner[1])
        .extrude(blind(0.25)));
    lid = lid.subtract(hole);
}
// The countersinks: the hole rims are the CLOSED circles at the top
// plane (chord ~0); the corner fillets' top arcs live at the same
// height and stay sharp.
lid = lid.chamfer(edges().circular().at_height(0.8).chord_length_between(0.0, 0.01), 0.03,
    #{ label: "countersinks" });

// Vent slots across the middle of the lid; each cut's label is read
// back so the "solid_lid" configuration can suppress exactly them.
let vent_cuts = [];
for i in [-2.0, -1.0, 0.0, 1.0, 2.0] {
    let vent = m.tool(|m| m.profile(xy_at(0.64)).rect(0.06, 0.5).at(i * 0.16, 0.0)
        .sketch("vent_" + i)
        .extrude(blind(0.25)));
    lid = lid.subtract(vent);
    let timeline = m.features();
    vent_cuts.push(timeline[timeline.len() - 1].label);
}

m.add_config("solid_lid", "no vent slots");
for cut in vent_cuts {
    m.set_suppression("solid_lid", m.feature(cut), true);
}
m.assert_watertight();
m
