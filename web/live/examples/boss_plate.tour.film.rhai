// A tour of the boss plate. A film reads the built model and never
// changes it: edit this script and press Build to replay it, or export
// it as a video from the film's bar.
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -35.0, elevation: 28.0, zoom: 1.0 } });
film.title("Boss plate", "A plate, a boss and G1 rounds", 2.5, |s| {
    s.orbit(30.0);
});
film.step("All the way round", "", 4.0, |s| {
    s.orbit(150.0);
});
film.step("Through the boss", "A section view: the model is not cut", 3.5, |s| {
    s.section(#{ from: "y", from_at: 1.0, to_at: 0.5 });
    s.camera(#{ azimuth: -75.0, elevation: 18.0, zoom: 1.2 });
});
film.step("The rounds", "Picked by geometry, certified G1", 3.5, |s| {
    s.section(#{ from: "y", from_at: 0.5, to: "off", duration: 1.5 });
    s.camera(#{ azimuth: 30.0, elevation: 42.0, zoom: 1.5 });
});
film.step("", "", 2.0, |s| {
    s.camera(#{ azimuth: -35.0, elevation: 28.0, zoom: 1.0 });
});
film
