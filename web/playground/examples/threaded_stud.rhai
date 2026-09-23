// Threaded stud: the smallest real-thread part — one verb cuts a
// modeled ISO M10×1.5 thread along a stud, with runout from the
// nominal tables. The coarse pitch keeps the tool coil short and the
// tip is left square (tip_chamfer: false), so this is the quickest
// way to watch a thread serve.
let m = model();

let stud = m.threaded_rod("M10x1.5", 8.0, #{ tip_chamfer: false, label: "stud" });

m.assert_watertight();
m
