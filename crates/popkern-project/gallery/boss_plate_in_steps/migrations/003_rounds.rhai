// Step 3: G1 rounds on the boss rim and the plate's edges.
let m = previous();
m.active().fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{ continuity: "g1", label: "rounds" });
m.assert_watertight();
m
