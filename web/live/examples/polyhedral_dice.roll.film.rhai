// The dice rolled onto the table one after another. The document is a
// part with eight bodies ("body 0" .. "body 7") in two rows along x;
// each rolls in along +x to where it was drawn, turning about a y axis
// through its centre by the distance over its rolling radius
// (millimetres). The die furthest along each row goes first, so every
// later one stops short of the ones already down (they rolled through
// each other in drawing order).
let centres = [
    [1.4, 1.4, 8.0], [26.0, 0.0, 8.0], [52.0, 0.0, 6.0], [0.0, 28.0, 7.6],
    [26.0, 28.0, 7.6], [52.0, 28.0, 8.6], [78.0, 28.0, 8.8], [78.0, 0.0, 8.0],
];
let distance = 120.0;
let film = storyboard(#{ look: "toon", camera: #{ azimuth: -75.0, elevation: 35.0, zoom: 1.05 } });
film.title("Polyhedral dice", "Seven classic dice and a pipped d6", 2.0, |s| {
    for k in 0..8 { s.show("body " + k, false); }
});
film.step("Roll them out", "", 6.0, |s| {
    let order = [7, 6, 2, 5, 1, 4, 0, 3];
    for i in 0..8 {
        let k = order[i];
        let c = centres[k];
        let at = 0.45 * i;
        let degrees = distance / c[2] * 57.29577951308232;
        s.show("body " + k, true, #{ at: at });
        s.move("body " + k, #{ from: [-distance, 0.0, 0.0], at: at, duration: 2.0, ease: "out" });
        s.turn("body " + k, axis_at(c, [0.0, 1.0, 0.0]), -degrees, 0.0, #{ at: at, duration: 2.0, ease: "out" });
    }
});
film.step("A double six?", "", 3.0, |s| {
    s.camera(#{ target: "body 1", azimuth: -30.0, elevation: 45.0, zoom: 1.0 });
    s.highlight("body 1", #{ at: 0.5, duration: 2.0 });
    s.highlight("body 7", #{ at: 0.5, duration: 2.0 });
});
film.step("", "", 2.5, |s| {
    s.camera(#{ target: "all", azimuth: -40.0, elevation: 32.0, zoom: 1.05 });
});
film
