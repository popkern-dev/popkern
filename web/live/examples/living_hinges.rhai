// Two living hinges: a straight flexure for flexible plastics (PP,
// PETG) — a thin web with rounded roots — and a slotted serpentine for
// rigid ones, where alternating kerf slots let a stiff plate bend.
let m = model();
let web = m.living_hinge("straight", 10.0, 40.0, #{ web: 0.4, label: "flexure" });
let kerf = m.living_hinge("slotted", 24.0, 40.0, #{ slot_pitch: 2.0, slot_width: 0.8, bridge: 3.0, label: "kerf" });
kerf = kerf.translate(0.0, 20.0, 0.0);
m.assert_watertight();
m
