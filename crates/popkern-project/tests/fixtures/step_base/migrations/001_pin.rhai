// Adds a pin part beside the imported cube.
let a = previous();
let pin = a.part("pin", |m| {
    m.profile(xy()).circle(0.1).sketch("pin_sketch").extrude(blind(0.5));
});
a.place(pin).translate(2.0, 0.0, 0.0);
a
