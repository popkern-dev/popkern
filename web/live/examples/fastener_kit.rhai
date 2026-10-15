// Fastener kit, as an ASSEMBLY: an M6 hex bolt with its nut, a
// socket-head cap screw, and a countersunk screw — real modeled ISO
// threads from the fastener feature family (PK-460). The three screws
// are ONE parametric part: a family declared with its default size,
// and each other size asked for with .variant(…), so every variant is
// built once from the same recipe. Designations carry
// millimetre-magnitude dimensions; heads, tip lead-ins, and thread
// runout come from the ISO nominal tables.
let a = assembly("fastener_kit");

let screw = a.part("screw", #{ size: "M6", length: 14.0, head: "hex", thread_length: 8.0 }, |m, p| {
    m.bolt(p.size, p.length, #{ head: p.head, thread_length: p.thread_length, label: "screw" });
});
let nut = a.part("nut", |m| {
    // Printed-fit clearance on the internal thread.
    m.nut("M6", #{ clearance: 0.2, label: "nut" });
});

// The hex bolt, laid on its side across the front.
a.place(screw, "hex_bolt")
    .rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 90.0)
    .translate(-9.0, -10.0, 5.0);

// Its nut, threaded on: the bolt's own frame, 7 mm (seven whole
// M6 leads) up the shank, so the internal thread mates the external
// one at the nut's clearance (PK-476's phase convention). It stood
// beside the bolt before, and the bolt's tip ran into its hexagon.
a.place(nut)
    .rotate(axis_at([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]), 90.0)
    .translate(-2.0, -10.0, 5.0);

// A socket-head cap screw standing on its head.
a.place(screw.variant(#{ size: "M4", length: 10.0, head: "socket", thread_length: 6.0 }), "cap_screw")
    .translate(-8.0, 6.0, 4.0);

// A countersunk screw, flush plane up.
a.place(screw.variant(#{ size: "M5", length: 12.0, head: "countersunk", thread_length: 7.0 }), "flat_head")
    .translate(8.0, 6.0, 0.0);

a.assert_watertight();
a
