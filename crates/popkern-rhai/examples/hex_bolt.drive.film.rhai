// The hex bolt driven up its own axis: head below, shank up the z axis
// through the origin. A right-hand thread advances along +z turning
// positively about it, one 1 mm lead a turn.
let axis = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
let travel = 10.0;
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -40.0, elevation: 20.0, zoom: 1.0 } });
film.title("Hex bolt", "M6 x 1, a modelled ISO thread", 2.0, |s| {
    s.orbit(20.0);
});
film.step("Driving in", "", 4.0, |s| {
    s.move("body 0", #{ from: [0.0, 0.0, -travel], ease: "linear" });
    s.turn("body 0", axis, -travel * 360.0, 0.0, #{ ease: "linear" });
});
film.step("The thread", "A section view: the model is not cut", 3.0, |s| {
    s.section(#{ from: "y", from_at: 1.0, to_at: 0.5 });
    s.camera(#{ azimuth: -90.0, elevation: 5.0, zoom: 1.4 });
});
film.step("", "", 1.5, |s| {
    s.section(#{ from: "y", from_at: 0.5, to: "off" });
    s.camera(#{ azimuth: -40.0, elevation: 20.0, zoom: 1.0 });
});
film
