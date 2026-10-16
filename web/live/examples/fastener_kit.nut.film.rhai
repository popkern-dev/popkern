// The kit's nut spun off its bolt and back on. The hex bolt lies along
// x on the line y = -10, z = 5 (tip at x = 5); the nut (5.2 thick) is
// threaded on it at x = -2. A right-hand nut walks toward the tip
// turning positively about +x, one 1 mm lead a turn, for the 7 mm until
// it leaves the thread; then it lifts clear without turning.
let shank = axis_at([0.0, -10.0, 5.0], [1.0, 0.0, 0.0]);
let travel = 7.0;
let turns = travel / 1.0 * 360.0;
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -60.0, elevation: 28.0, zoom: 1.0 } });
film.title("Fastener kit", "Modelled ISO threads: a bolt and nut, a cap screw, a flat head", 2.5, |s| {
    s.orbit(20.0);
});
film.step("Spin the nut off", "", 4.0, |s| {
    s.move("nut_1", #{ to: [travel, 0.0, 0.0], ease: "linear" });
    s.turn("nut_1", shank, 0.0, turns, #{ ease: "linear" });
    s.camera(#{ azimuth: -100.0, elevation: 20.0, zoom: 1.5, target: "hex_bolt" });
});
film.step("Off the end", "", 1.5, |s| {
    s.move("nut_1", #{ to: [5.0, 0.0, 4.0] });
});
film.step("And back on", "", 4.5, |s| {
    s.move("nut_1", #{ to: [-5.0, 0.0, -4.0], duration: 1.0 });
    s.move("nut_1", #{ to: [-travel, 0.0, 0.0], at: 1.0, duration: 3.5, ease: "linear" });
    s.turn("nut_1", shank, 0.0, -turns, #{ at: 1.0, duration: 3.5, ease: "linear" });
});
film.step("The other two", "A socket-head cap screw and a countersunk screw", 3.0, |s| {
    s.highlight("cap_screw", #{ duration: 1.5 });
    s.highlight("flat_head", #{ at: 1.5, duration: 1.5 });
    s.camera(#{ azimuth: -40.0, elevation: 30.0, zoom: 1.0, target: "all" });
});
film
