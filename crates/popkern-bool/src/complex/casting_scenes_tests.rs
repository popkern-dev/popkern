//! The near-tangent cap casting: the boss's cap cuts the ring wall
//! next to its extreme, leaving a lens of ring wall inside the boss.

use super::Casting;

#[test]
fn the_near_tangent_cap_casting_serves_through_the_rim_extreme() {
    for cap_x in [
        -368.0,
        -367.5,
        -367.3,
        -367.1,
        -367.05,
        -367.0 + 1e-6,
        -367.001,
        -366.999,
        -366.9,
        // The OCC reference fuses' caps, and the extreme itself.
        -367.002,
        -366.998,
        -367.0 - 1e-6,
        -367.0,
        // The cap's rim arc over the lens, removed from the boss wall
        // as a duplicate of the ring-wall record 1e-6 off it, kept by
        // the cap: one edge of the complex.
        -367.00001,
        -367.0001,
    ] {
        Casting::near_tangent_cap(cap_x).complex_union_serves();
    }
}
