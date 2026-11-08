fn ball(m, r) {
    m.sphere(r)
}

// A cylinder of radius 5 on the z axis from z = 0 to 10, with a ball of radius 5 at (0, 0, 10),
// the middle of its top, cut from it.
let m = model();
let dish = m.profile(xy()).circle(5.0).sketch("cylinder").extrude(blind(10.0));
let scoop = m.tool(|m| ball(m, 5.0).translate(0.0, 0.0, 10.0));
dish = dish.subtract(scoop);
m
