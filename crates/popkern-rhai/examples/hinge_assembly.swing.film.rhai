// The print-in-place hinge swung through its range by its joint: the
// assembly's `lid_angle` parameter only moves the leaf, so the film
// sweeps it (placement-only, nothing rebuilds).
let doc = document();
let open = doc.value("lid_angle");
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -55.0, elevation: 24.0, zoom: 1.0 } });
film.title("Print-in-place hinge", "Two leaves mated about the pin axis", 2.5, |s| {
    s.orbit(20.0);
});
film.step("Closed", "", 2.5, |s| {
    s.joint("lid_angle", open, 5.0);
});
film.step("Past flat", "The swing check finds where the leaves meet", 3.5, |s| {
    s.joint("lid_angle", 5.0, 230.0);
    s.orbit(-30.0);
});
film.step("Back to the drawn angle", "", 2.5, |s| {
    s.joint("lid_angle", 230.0, open);
});
film
