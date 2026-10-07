// A piano (continuous) hinge: a butt hinge with many knuckles, here 25
// along 200. The other leaf's knuckle zones go in one cut whatever the
// count, so the hinge builds in about the time of five knuckles' worth
// of booleans per leaf.
let m = model();
m.hinge("butt", 200.0, #{ knuckles: 25, leaf_width: 16.0, thickness: 1.2, pin: 2.0, clearance: 0.15, label: "piano" });
m.assert_watertight();
m
