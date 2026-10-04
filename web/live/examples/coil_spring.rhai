// Three springs from one verb: a compression spring, a tapered valve
// spring, and a left-hand spring phased a quarter turn. Every path is
// an exact rational helix — quarter arcs with true weights — so the
// coil radius never drifts, and the circular wire section keeps its
// exact curvature (four smooth faces per coil, not a faceted tube).
let m = model();
m.param("wire_r", 0.04, "m");

// The compression spring: 6.5 turns at constant radius and pitch,
// with GROUND flat ends — each grind plane shaves the last part-turn
// down to the seating flat a real compression spring carries
// (PK-441/PK-444: the cut planes cross the tube wall's seam;
// PK-445: the cut closed wall meshes watertight on the
// curvature-adaptive lattice).
let compression = m.profile(xy()).circle(p("wire_r")).sketch("comp_wire")
    .helical_sweep(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
                   0.3, 0.15, 6.5, #{ label: "compression_coil" });
let grind_bot = m.tool(|m| m.profile(xy_at(-0.08)).rect(1.0, 1.0)
    .sketch("grind_bot").extrude(blind(0.1)));
let compression = compression.subtract(grind_bot);
let grind_top = m.tool(|m| m.profile(xy_at(0.955)).rect(1.0, 1.0)
    .sketch("grind_top").extrude(blind(0.2)));
let compression = compression.subtract(grind_top);

// The valve spring: the coil tightens from 0.3 to 0.16 as it climbs.
let valve = m.tool(|m| m.profile(xy()).circle(p("wire_r")).sketch("valve_wire")
    .helical_sweep(axis_at([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
                   0.3, 0.18, 5.0,
                   #{ taper_to: 0.16, label: "valve_coil" }));

// The left-hand spring: mirrored winding, started a quarter turn
// around so the ends face away from its neighbours.
let left = m.tool(|m| m.profile(xy()).circle(p("wire_r")).sketch("left_wire")
    .helical_sweep(axis_at([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
                   0.3, 0.15, 6.5,
                   #{ left_hand: true, start_angle: 90.0, label: "left_coil" }));

m.assert_watertight();
m
