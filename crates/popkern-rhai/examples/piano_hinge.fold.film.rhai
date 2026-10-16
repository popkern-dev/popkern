// The piano hinge folding: "body 1" is the moving leaf, turning about
// the pin ("body 2") up the z axis through the origin.
let pin = axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -50.0, elevation: 22.0, zoom: 0.8 } });
film.title("Piano hinge", "Twenty-five knuckles along 200", 2.5, |s| {
    s.orbit(15.0);
});
film.step("Folding", "The leaf turns about the pin", 3.5, |s| {
    s.turn("body 1", pin, 0.0, 165.0);
});
film.step("And open flat again", "", 3.5, |s| {
    s.turn("body 1", pin, 0.0, -165.0);
    s.camera(#{ azimuth: -20.0, elevation: 50.0, zoom: 1.1 });
});
film.step("Knuckle by knuckle", "", 3.0, |s| {
    s.camera(#{ azimuth: -80.0, elevation: 10.0, zoom: 2.4, target: "body 2" });
    s.turn("body 1", pin, 0.0, 90.0, #{ at: 0.5, duration: 2.5 });
});
film
