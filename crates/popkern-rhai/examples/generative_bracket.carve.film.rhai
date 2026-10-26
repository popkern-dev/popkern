// The bracket's design converging (PK-959). The solver starts from the
// block at a quarter density everywhere, and each update moves material
// to where it carries the load: the film shows the design at every
// update its solve recorded (the first dozen, then about every 15%),
// blended between them, over the block it is carved from. The
// intermediate designs are display only; the part is the faceted body
// that fades in at the end. Edit and press Build to replay, or export
// it as a video from the film's bar.
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -50.0, elevation: 22.0, zoom: 1.05 } });
film.title("Generative bracket", "A 50 mm block bolted to the wall, 200 N down on its eye", 3.0, |s| {
    s.show("body 0", false);
    s.orbit(15.0);
});
film.step("Carving", "Each update moves material to where the load runs", 9.0, |s| {
    s.evolve("carve");
    s.orbit(75.0);
});
film.step("The bracket", "The design as a body: planar faces, the wall, bores and eye kept", 4.0, |s| {
    s.fade("design carve", 1.0, 0.0, #{ duration: 1.2 });
    s.show("body 0", true);
    s.fade("body 0", 0.0, 1.0, #{ duration: 1.2 });
    s.orbit(45.0);
});
film.hold(1.0);
film
