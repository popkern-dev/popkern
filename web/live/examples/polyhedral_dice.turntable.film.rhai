// The dice set on a turntable: each die hops in turn as the camera
// goes round. The document is a part with eight bodies (the seven
// numbered dice, then the pipped d6), named "body 0" .. "body 7".
let doc = document();
let film = storyboard(#{ look: "toon", camera: #{ azimuth: -70.0, elevation: 32.0, zoom: 1.05 } });
let names = ["d4", "d6", "d8", "d10", "d%", "d12", "d20", "pipped d6"];

film.title("Polyhedral dice", "Seven classic dice and a pipped d6", 2.5, |s| {
    s.orbit(25.0);
});
film.step("Every die in turn", "Engraved numerals, chamfered edges", 6.0, |s| {
    s.orbit(120.0);
    for k in 0..doc.bodies {
        let at = 0.3 + 0.65 * k;
        s.move("body " + k, #{ to: [0.0, 0.0, 9.0], at: at, duration: 0.35, ease: "out" });
        s.move("body " + k, #{ to: [0.0, 0.0, -9.0], at: at + 0.35, duration: 0.35, ease: "in" });
        s.highlight("body " + k, #{ at: at, duration: 0.7 });
    }
});
film.step("The d20", "Opposite faces sum to 21", 3.0, |s| {
    s.camera(#{ target: "body 6", zoom: 1.0, azimuth: 20.0, elevation: 30.0 });
});
film.step("The whole set", "", 3.0, |s| {
    s.camera(#{ target: "all", zoom: 1.05, azimuth: -40.0, elevation: 32.0 });
    s.orbit(20.0);
});
film
