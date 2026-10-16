// The hinged case opening and closing: `lid_angle` is the assembly's
// joint between the box and the lid, so the film sweeps it.
let doc = document();
let open = doc.value("lid_angle");
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -35.0, elevation: 30.0, zoom: 1.0 } });
film.title("Hinged case", "A box and its lid, hinged where they meet", 2.5, |s| {
    s.orbit(25.0);
});
film.step("Shut", "", 2.0, |s| {
    s.joint("lid_angle", open, 1.0);
});
film.step("Open wide", "The pin is printed in place inside the box's knuckles", 3.5, |s| {
    s.joint("lid_angle", 1.0, 110.0);
    s.camera(#{ azimuth: 150.0, elevation: 25.0, zoom: 1.3 });
});
film.step("Through the knuckles", "A section view: the model is not cut", 3.0, |s| {
    s.section(#{ from: "x", from_at: 1.0, to_at: 0.62, duration: 1.0 });
    s.camera(#{ azimuth: 15.0, elevation: 15.0, zoom: 1.6, target: "lid" });
});
film.step("", "", 2.5, |s| {
    s.section(#{ from: "x", from_at: 0.62, to: "off", duration: 1.0 });
    s.joint("lid_angle", 110.0, open);
    s.camera(#{ azimuth: -35.0, elevation: 30.0, zoom: 1.0 });
});
film
