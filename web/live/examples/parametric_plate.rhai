// A plate driven by parameters and an expression, with a named
// configuration that overrides the wall: popkern live's sliders
// edit `wall` and `span`; the "heavy" configuration doubles the wall.
let m = model();
m.param("wall", 0.3, "m");
m.param("span", 2.0, "m");
m.param_expr("lip", "wall * 0.5", "m");
let plate = m.profile(xy()).rect(p("span"), 1.2).sketch("plate_sketch").extrude(blind(p("wall")));
plate = plate.fillet(edges().straight().direction("z"), 0.1, #{ label: "corners" });
let rib = m.tool(|m| m.profile(xy_at(p("wall"))).rect(expr("span * 0.8"), 0.12).sketch("rib_sketch")
    .extrude(blind(p("lip"))));
plate = plate.union(rib);
m.add_config("heavy", "double wall");
m.set_override("heavy", "wall", 0.6);
m.assert_watertight();
m
