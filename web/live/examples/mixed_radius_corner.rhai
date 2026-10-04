// One corner, three radii. `fillet_each` gives each group of edges its
// own radius in a single fillet, and `at_vertex` picks the three edges
// meeting at the block's top corner. With radii this different the
// setback G2 cap bends sharply inside, so the kernel serves the
// variational corner fit: three sectors, curvature-continuous across
// their spokes and certified G2 to the rounds. Try equal radii
// to see the setback Gregory cap instead, or drop `continuity`.
// Add `corner_cap: "fitted"` (or "setback", "classic") to choose the
// corner's construction; with `corner_cap_strict: true` the fillet
// fails rather than fall back. The report's `blend.corners` says what
// closed the corner.
let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));

let corner = [0.5, 0.5, 1.0];
block.fillet_each([
    [edges().at_vertex(corner).direction("x"), 0.08],
    [edges().at_vertex(corner).direction("y"), 0.2],
    [edges().at_vertex(corner).direction("z"), 0.3],
], #{ continuity: "g2", label: "mixed_corner" });

m.assert_watertight();
m
