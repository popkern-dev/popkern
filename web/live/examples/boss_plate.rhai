// The boss-plate showcase in the popkern Rhai dialect (phase 3
// acceptance: mesh-identical to the Rust recorder and to the direct
// popkern-part build).
let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch")
             .extrude(blind(0.4));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
             .extrude(blind(0.8)));
plate = plate.union(boss);
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{ continuity: "g1", label: "rounds" });
m.assert_watertight();
print(m.feature("rounds").describe());
m
