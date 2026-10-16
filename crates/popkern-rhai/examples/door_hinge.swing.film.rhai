// The butt hinge swinging on its pin. The document is a part with
// three bodies — "body 0" the fixed leaf, "body 1" the moving leaf,
// "body 2" the pin — and the pin runs up the z axis through the origin,
// so the moving leaf turns about it.
let pin = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -60.0, elevation: 25.0, zoom: 1.0 } });
film.title("Butt hinge", "Two leaves, five knuckles and a pin", 2.5, |s| {
    s.orbit(20.0);
});
film.step("Closing", "The leaf turns about the pin", 3.0, |s| {
    s.turn("body 1", pin, 0.0, 140.0);
});
film.step("Opening", "", 3.0, |s| {
    s.turn("body 1", pin, 0.0, -170.0);
});
film.step("The pin", "A separate part through every knuckle", 3.0, |s| {
    s.fade("body 0", 1.0, 0.25, #{ duration: 1.0 });
    s.fade("body 1", 1.0, 0.25, #{ duration: 1.0 });
    s.highlight("body 2", #{ at: 0.6, duration: 2.0 });
    s.turn("body 1", pin, 0.0, 30.0, #{ at: 1.0, duration: 2.0 });
});
film.step("", "", 1.5, |s| {
    s.fade("body 0", 0.25, 1.0, #{ duration: 0.8 });
    s.fade("body 1", 0.25, 1.0, #{ duration: 0.8 });
});
film
