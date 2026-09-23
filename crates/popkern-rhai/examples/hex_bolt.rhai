// Hex bolt: a single short M6 hex bolt — head, shank and a real
// modeled ISO thread from one verb. Half the shank is left plain
// (thread_length) and the tip is left square (tip_chamfer: false)
// to keep the build quick; see the fastener kit for the full-dress
// versions.
let m = model();

let bolt = m.bolt("M6", 10.0, #{
    thread_length: 5.0,
    tip_chamfer: false,
    label: "bolt" });

m.assert_watertight();
m
