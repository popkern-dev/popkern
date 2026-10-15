// A section sweeping through the boss plate: the cut slides across,
// tilts, and lifts away, while the camera turns. The film reads the
// project's built document and never changes it.
let film = storyboard(#{ look: "toon", camera: #{ azimuth: -40.0, elevation: 30.0 } });
film.title("Boss plate", "A plate, a boss, rounds selected by geometry", 2.5, |s| {
    s.orbit(25.0);
});
film.step("Section across the boss", "Display only: the model is not cut", 3.5, |s| {
    s.section(#{ from: "x", from_at: 1.0, to_at: 0.5 });
    s.orbit(25.0);
});
film.step("Tilted", "", 2.5, |s| {
    s.section(#{ from: "x", from_at: 0.5, to: #{ normal: [0.7, 0.0, 0.7], offset: 0.2 } });
    s.camera(#{ azimuth: 20.0, elevation: 20.0, zoom: 1.3 });
});
film.step("Whole again", "", 2.5, |s| {
    s.section(#{ from: #{ normal: [0.7, 0.0, 0.7], offset: 0.2 }, to: "off" });
    s.camera(#{ azimuth: -40.0, elevation: 30.0, zoom: 1.0 });
});
film
