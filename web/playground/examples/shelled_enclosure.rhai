// An enclosure: a rounded block hollowed from the top, leaving a
// uniform wall. The shell keeps the rounds on the outside and the
// inside alike.
let m = model();
let box = m.profile(xy()).rect(1.6, 1.0).sketch("outline").extrude(blind(0.6));
box = box.fillet(edges().straight().direction("z"), 0.15, #{ label: "corners" });
box = box.shell(faces().facing("z").at_height(0.6), 0.06, #{ label: "hollow" });
m.assert_watertight();
m
