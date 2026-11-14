// Two cylinders of radius 5 and length 30 crossing at right angles: along the x axis, and
// parallel to y on the line x = 0, z = 9.99, each extruded both ways from its middle. Their
// common part, 0.01 thick.
let m = model();
let along_x = m.profile(yz()).circle(5.0).sketch("along_x").extrude(mid_plane(15.0));
let along_y = m.tool(|m| m.profile(xz()).circle(5.0).at(0.0, 9.99)
    .sketch("along_y").extrude(mid_plane(15.0)));
along_x = along_x.intersect(along_y);
m
