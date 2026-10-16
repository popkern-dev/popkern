// The 4:1 gearbox running: the ring held, the sun driving, the carrier
// out. With 9-tooth sun and planets in a 27-tooth ring the carrier
// turns a quarter of the sun; each planet rolls round the ring, spinning
// on its pin three quarters of the sun's angle backwards relative to
// the carrier (half the sun's angle backwards in all). Every turn is
// about +z: the sun's and carrier's through the origin, each planet's
// spin through its pin at radius 0.5.
let z = [0.0, 0.0, 1.0];
let centre = axis_at([0.0, 0.0, 0.0], z);
let pins = [[0.5, 0.0, 0.0], [-0.25, 0.4330127, 0.0], [-0.25, -0.4330127, 0.0]];
// The gears through `sun` degrees of sun rotation (called as a method,
// `s.run(...)`, so `this` is the step). Turns compose in the order they
// are declared, and a planet's spin is about its pin where it was
// drawn, so every spin must come before every orbit: one spin and one
// orbit per planet, spanning the whole run (a second step's spin pivoted
// about the old pin after the first orbit had moved it, and flung the
// planets out of the ring).
fn run(pins, centre, z, sun, options) {
    this.turn("sun_1", centre, 0.0, sun, options);
    this.turn("carrier_1", centre, 0.0, sun / 4.0, options);
    for k in 0..3 {
        this.turn("planet_" + k, axis_at(pins[k], z), 0.0, -0.75 * sun, options);
        this.turn("planet_" + k, centre, 0.0, sun / 4.0, options);
    }
}
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -50.0, elevation: 50.0, zoom: 1.0 } });
film.title("Planetary gearbox", "4:1 — the ring held, the sun drives the carrier", 2.5, |s| {
    s.orbit(15.0);
});
film.step("The sun drives", "Two turns of the sun: half a turn of the carrier", 6.0, |s| {
    // 15 s of running from here: 1710° of sun, easing in and out.
    s.run(pins, centre, z, 1710.0, #{ duration: 15.0 });
});
film.step("Under the carrier", "The planets roll round the ring", 6.0, |s| {
    s.fade("carrier_1", 1.0, 0.2, #{ duration: 1.0 });
    s.camera(#{ azimuth: -20.0, elevation: 70.0, zoom: 1.15 });
});
film.step("", "", 3.0, |s| {
    s.fade("carrier_1", 0.2, 1.0, #{ at: 1.5, duration: 1.0 });
    s.camera(#{ azimuth: -50.0, elevation: 50.0, zoom: 1.0 });
});
film
