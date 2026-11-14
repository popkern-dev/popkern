// The cube [-5, 5]³ united with a cylinder of radius 5 on the x axis from x = -10 to 10, its seam
// turned 45 degrees off the four lines where it touches the cube's faces.
let m = model();
let cube = m.profile(xy_at(-5.0)).rect(10.0, 10.0).sketch("cube").extrude(blind(10.0));
let rod = m.tool(|m| m.profile(yz_at(-10.0)).circle(5.0).sketch("rod").extrude(blind(20.0))
    .rotate(axis_at([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 45.0));
cube = cube.union(rod);
m
