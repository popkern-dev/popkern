// Step 2: a boss, built as a tool body and united with the plate.
// Every later migration starts from previous().
let m = previous();
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch").extrude(blind(0.8)));
m.active().union(boss);
m
