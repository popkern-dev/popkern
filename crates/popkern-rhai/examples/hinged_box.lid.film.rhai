// The printable hinged box: the lid ("body 1") turns about the hinge
// pin, which runs along x at y = 28.1, z = 30 — the seam's back top
// edge, pushed out by the knuckles' radius.
let pin = axis_at([0.0, 28.1, 30.0], [1.0, 0.0, 0.0]);
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -40.0, elevation: 30.0, zoom: 1.0 } });
film.title("Hinged box", "Printed in one piece: the pin prints in place", 2.5, |s| {
    s.orbit(20.0);
});
film.step("Lift the lid", "", 3.0, |s| {
    s.turn("body 1", pin, 0.0, -105.0);
});
film.step("Inside", "A 2.5 shell", 2.5, |s| {
    s.camera(#{ azimuth: -15.0, elevation: 55.0, zoom: 1.15 });
});
film.step("And shut", "", 2.5, |s| {
    s.turn("body 1", pin, 0.0, 105.0);
    s.camera(#{ azimuth: -40.0, elevation: 30.0, zoom: 1.0 });
});
film
