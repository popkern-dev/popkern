// The flat-head screw driven into its countersunk seat. The document
// is a part: "body 0" the plate (seat axis at the origin, its top at
// z = 4), "body 1" the M5 x 10 screw standing on its head beside it,
// on the vertical axis x = 18. The screw flips head-up and travels over
// the seat; the seat is a clearance hole, so it drops straight through
// until its tip leaves the plate's underside, then turns as its
// right-hand thread bites the part below: 0.8 a turn for the last 6 mm.
let flip = axis_at([18.0, 0.0, 5.0], [1.0, 0.0, 0.0]);
let spin = axis_at([18.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
let pitch = 0.8;
let drive = 6.0;
let turns = -drive / pitch * 360.0;
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -40.0, elevation: 30.0, zoom: 1.0 } });
film.title("Countersunk joint", "An M5 flat-head screw and its seat", 2.0, |s| {
    s.orbit(15.0);
});
film.step("Head up", "", 2.0, |s| {
    s.turn("body 1", flip, 0.0, 180.0);
    s.move("body 1", #{ to: [0.0, 0.0, 12.0] });
});
film.step("Over the seat", "", 1.8, |s| {
    s.move("body 1", #{ to: [-18.0, 0.0, 0.0] });
});
film.step("Drive it home", "Straight through the seat; clockwise once the thread bites", 4.5, |s| {
    s.move("body 1", #{ to: [0.0, 0.0, -12.0], duration: 1.0 });
    s.move("body 1", #{ to: [0.0, 0.0, -drive], at: 1.0, duration: 3.5, ease: "linear" });
    s.turn("body 1", spin, 0.0, turns, #{ at: 1.0, duration: 3.5, ease: "linear" });
    s.camera(#{ azimuth: -25.0, elevation: 18.0, zoom: 1.5, target: "body 0" });
});
film.step("Flush in its seat", "A section view: the model is not cut", 3.0, |s| {
    s.section(#{ from: "y", from_at: 1.0, to_at: 0.5, duration: 1.0 });
    s.camera(#{ azimuth: 90.0, elevation: 8.0, zoom: 1.9, target: "body 0" });
});
film.step("And back out", "", 4.0, |s| {
    s.section(#{ from: "y", from_at: 0.5, to: "off", duration: 1.0 });
    s.move("body 1", #{ to: [0.0, 0.0, drive], at: 0.5, duration: 2.5, ease: "linear" });
    s.turn("body 1", spin, 0.0, -turns, #{ at: 0.5, duration: 2.5, ease: "linear" });
    s.move("body 1", #{ to: [0.0, 0.0, 12.0], at: 3.0, duration: 1.0 });
    s.camera(#{ azimuth: -40.0, elevation: 30.0, zoom: 1.0 });
});
film
