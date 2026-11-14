// A cylinder of radius 5 along the x axis, and one of radius 3 parallel to y on the line x = 0,
// z = 2, both 30 long and extruded both ways from their middles: their tops touch at (0, 0, 5).
// Their common part.
let m = model();
let main = m.profile(yz()).circle(5.0).sketch("main").extrude(mid_plane(15.0));
let branch = m.tool(|m| m.profile(xz()).circle(3.0).at(0.0, 2.0)
    .sketch("branch").extrude(mid_plane(15.0)));
main = main.intersect(branch);
m
