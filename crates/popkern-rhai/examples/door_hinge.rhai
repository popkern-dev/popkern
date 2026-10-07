// A 76 butt hinge for a cabinet door: two leaves with five knuckles,
// a separate pin, and three countersunk M3 screw holes per leaf.
// Leave out any size and it follows the length; open it with `angle:`.
let m = model();
let h = m.hinge("butt", 76.0, #{ knuckles: 5, screw: "M3", holes: 3, angle: 150.0, label: "door_hinge" });
print(h.leaf_a.is_alive() && h.leaf_b.is_alive() && h.pin.is_alive());
m.assert_watertight();
m
