// The hinged box as an ASSEMBLY: the box and its lid are separate
// parts, hinged along the seam where they meet — each part grows its
// own knuckles, the pin is printed in place inside the box's — and
// `lid_angle` opens the lid: a slider in popkern live that swings it
// without rebuilding either part. The swing check reports where the
// opening lid first meets the box.
let a = assembly("hinged_case");
let boxp = a.part("box", |m| {
    let b = m.profile(xy()).rect(80.0, 50.0).sketch("base_outline").extrude(blind(30.0));
    b = b.fillet(edges().straight().direction("z"), 6.0, #{ label: "base_corners" });
    b.shell(faces().facing("z").at_height(30.0), 2.5, #{ label: "cavity" });
});
let lidp = a.part("lid", |m| {
    let l = m.profile(xy_at(30.0)).rect(80.0, 50.0).sketch("lid_outline").extrude(blind(3.0));
    l.fillet(edges().straight().direction("z"), 6.0, #{ label: "lid_corners" });
});
let bx = a.place(boxp, "box");
let ld = a.place(lidp, "lid");
a.ground(bx);
let h = a.hinge_bodies("lid", "print_in_place", bx, ld, axis_at([-25.0, 25.0, 30.0], [1.0, 0.0, 0.0]),
    #{ length: 50.0, knuckles: 5, pin: 2.4, wall: 1.2, clearance: 0.35, angle: 70.0 });
let swing = a.swing("lid_angle", #{ from: 1.0, to: 300.0, step: 10.0 });
print("the lid first meets the box at " + swing.first_contact + "°");
a
