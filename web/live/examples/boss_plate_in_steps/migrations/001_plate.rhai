// Step 1: the plate. The first migration builds the document.
let m = model();
m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch").extrude(blind(0.4));
m
