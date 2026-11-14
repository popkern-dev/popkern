// A torus about the z axis, its middle circle of radius 5 and its tube of radius 3, cut by the
// bitangent plane 3x - 4z = 0, keeping what lies on and below it, z <= 3x/4. The tube's section is
// turned from the plane x = 0, on the side y > 0, so the seam is away from the tangent points.
let m = model();
let torus = m.profile(yz()).circle(3.0).at(5.0, 0.0).sketch("tube")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let below = m.tool(|m| m.profile(plane_at([0.0, 0.0, 0.0], [0.6, 0.0, -0.8])).rect(40.0, 40.0)
    .sketch("half_space").extrude(blind(20.0)));
torus = torus.intersect(below, #{ label: "cut" });
m
