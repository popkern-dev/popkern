// A print-in-place hinge as an ASSEMBLY: the two leaves are parts mated
// about the pin axis, and `lid_angle` opens it — a slider in popkern
// live that moves the leaf without rebuilding anything. The swing check
// sweeps the angle and reports where the leaves first collide.
let a = assembly("hinge_assembly");
let h = a.hinge("lid", "print_in_place", 60.0, #{ knuckles: 5, angle: 120.0 });
a.limits("lid_angle", 0.0, 270.0);
let swing = a.swing("lid_angle", #{ from: 180.0, to: 270.0, step: 15.0 });
print("first contact opening past flat: " + swing.first_contact);
a
