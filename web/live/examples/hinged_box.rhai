// A printable box whose lid is hinged to it: the knuckles grow out of
// the box's back wall and the lid's back edge along the seam where they
// meet, and the pin is part of the box's knuckles, printed in place
// through the lid's bored knuckles. Millimetre-magnitude numbers.
let m = model();

let base = m.profile(xy()).rect(80.0, 50.0).sketch("base_outline").extrude(blind(30.0));
base = base.fillet(edges().straight().direction("z"), 6.0, #{ label: "base_corners" });
base = base.shell(faces().facing("z").at_height(30.0), 2.5, #{ label: "cavity" });

let lid = m.profile(xy_at(30.0)).rect(80.0, 50.0).sketch("lid_outline").extrude(blind(3.0));
lid = lid.fillet(edges().straight().direction("z"), 6.0, #{ label: "lid_corners" });

// The seam: the back top edge (y = 25, z = 30), a 50 hinge centred on it.
let h = m.hinge_bodies("print_in_place", base, lid, axis_at([-25.0, 25.0, 30.0], [1.0, 0.0, 0.0]),
    #{ length: 50.0, knuckles: 5, pin: 2.4, wall: 1.2, clearance: 0.35, label: "lid_hinge" });
m.assert_watertight();
m
