// The two living hinges are single bodies — they bend rather than
// turn on a pin — so the film shows what makes them bend: the thin web
// of the straight flexure ("body 0") and the kerf slots of the
// serpentine ("body 1"), cut through by a section view.
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -35.0, elevation: 35.0, zoom: 1.0 } });
film.title("Living hinges", "A flexure web and a slotted kerf", 2.5, |s| {
    s.orbit(20.0);
});
film.step("The flexure's web", "0.4 thick, rounded at the roots", 3.5, |s| {
    s.section(#{ from: "z", from_at: 1.0, to_at: 0.5 });
    s.camera(#{ target: "body 0", azimuth: 0.0, elevation: 75.0, zoom: 1.6 });
});
film.step("The kerf's slots", "Alternating slots let a rigid plate bend", 3.5, |s| {
    s.camera(#{ target: "body 1", azimuth: -10.0, elevation: 65.0, zoom: 1.3 });
});
film.step("", "", 2.5, |s| {
    s.section(#{ from: "z", from_at: 0.5, to: "off", duration: 1.0 });
    s.camera(#{ target: "all", azimuth: -35.0, elevation: 35.0, zoom: 1.0 });
});
film
