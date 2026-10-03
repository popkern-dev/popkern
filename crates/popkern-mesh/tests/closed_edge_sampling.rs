//! PK-636 regression pins for imported-solid watertightness.
//!
//! - A closed edge (a full circle bounding a face on its own) smaller
//!   than the chord height must still sample as a ring: it once stayed
//!   one zero-chord interval, and the faces it bounded fell off the
//!   boundary-first path and meshed open (clevis21.stp's drill-point
//!   circles).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate, welded_open_edges};

#[test]
fn a_cylinder_thinner_than_the_chord_height_meshes_watertight() {
    let options = TessellationOptions::default();
    for radius in [
        0.2 * options.chord_height,
        0.45 * options.chord_height,
        2.0 * options.chord_height,
    ] {
        let pin = BRep::cylinder(radius, 0.05, Tolerance::default()).expect("cylinder");
        let mesh = tessellate(&pin, &options).expect("mesh");
        assert_eq!(welded_open_edges(&mesh), 0, "radius {radius}");
        assert!(
            mesh.zipped_t_junctions == 0,
            "radius {radius}: closed by repair, not construction"
        );
    }
}
