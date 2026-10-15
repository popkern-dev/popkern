// Fastener kit: an M6 hex bolt with its nut, a socket-head cap screw,
// and a countersunk screw — real modeled ISO threads from the
// fastener feature family (PK-460). Designations carry
// millimetre-magnitude dimensions; heads, tip lead-ins, and thread
// runout come from the ISO nominal tables.
let m = model();

// The hex bolt, laid on its side across the front.
let bolt = m.bolt("M6", 14.0, #{ thread_length: 8.0, label: "hex_bolt" });
bolt = bolt.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 90.0);
bolt = bolt.translate(-9.0, -10.0, 5.0);

// Its nut, printed-fit clearance on the internal thread, threaded on
// 7 mm (seven whole M6 leads) up the bolt's shank, as the gallery kit.
let nut = m.nut("M6", #{ clearance: 0.2, label: "nut" });
nut = nut.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 90.0);
nut = nut.translate(-2.0, -10.0, 5.0);

// A socket-head cap screw standing on its head.
let cap = m.bolt("M4", 10.0, #{ head: "socket", thread_length: 6.0, label: "cap_screw" });
cap = cap.translate(-8.0, 6.0, 4.0);

// A countersunk screw, flush plane up.
let flat = m.bolt("M5", 12.0, #{ head: "countersunk", thread_length: 7.0, label: "flat_head" });
flat = flat.translate(8.0, 6.0, 0.0);

m.assert_watertight();
m
