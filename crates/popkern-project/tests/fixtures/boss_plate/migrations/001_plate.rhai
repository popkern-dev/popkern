// The plate.
let m = model();
m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch").extrude(blind(0.4));
m
