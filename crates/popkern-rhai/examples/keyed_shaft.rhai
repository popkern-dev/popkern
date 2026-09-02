// A shaft with a keyway: the slot is a moved and turned tool body,
// and the shoulder is chamfered.
let m = model();
let shaft = m.profile(xy()).circle(0.3).sketch("shaft_sketch").extrude(blind(1.6));
let key = m.tool(|m| m.profile(xy_at(1.0)).rect(0.12, 0.5).sketch("key_sketch").extrude(blind(0.7)));
key = key.translate(0.0, 0.28, 0.0);
key = key.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 30.0);
shaft = shaft.subtract(key);
shaft = shaft.chamfer(edges().circular().at_height(0.0), 0.04, #{ label: "lead_in" });
m.assert_watertight();
m
