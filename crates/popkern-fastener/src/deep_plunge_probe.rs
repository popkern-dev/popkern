//! PK-481 pinned fixture (promoted from the PK-480 tilt-bisection
//! probe): a deep-plunged entry cone on a PERPENDICULAR plate.
//!
//! History: PK-480's angled tap stranded on cone×flank triple points;
//! this probe proved the tilt was a bystander (the interior-wrap
//! junction class reproduces flat) and drove the junction-exactness
//! chain (triple-point refinement, open micro-extension, kept-end
//! weld authority). Those landed; the thread stage now assembles
//! clean through CELLULAR — which leaves the bore wall as ONE
//! multi-wrap helical band face, and the chamfer stage's closed
//! cone∩cylinder circle record clips to nothing against it (the
//! containment vote folds by a single global period shift). That
//! remaining class is PK-481; un-ignore this test when it lands.

use popkern_core::tolerance::Tolerance;
use popkern_sketch::ProfileShape;

use crate::engine::{cut_thread, extrude_shape, revolved_tool, subtract};
use crate::entry_ring::{EntryProbe, entry_ring_profile};
use crate::spec::ThreadSpec;
use crate::thread_form::{self, internal_bore_radius, internal_ridge};

#[test]
#[ignore = "blocked on PK-481: closed records vs multi-wrap helical band faces"]
fn deep_plunge_chamfer_on_perpendicular_plate() {
    let tolerance = Tolerance::default();
    let spec = ThreadSpec::parse("M6").expect("M6");
    let clearance = 0.15;
    let ridge = internal_ridge(&spec, clearance).expect("ridge");
    let r_bore = internal_bore_radius(&spec, clearance);
    let plate = extrude_shape(
        ProfileShape::Rect {
            width: 20.0,
            height: 20.0,
        },
        0.0,
        8.0,
        "probe plate",
        tolerance,
    )
    .expect("plate");
    let up = 0.5 * spec.major_diameter() + spec.lead();
    let bore = extrude_shape(
        ProfileShape::Circle { radius: r_bore },
        -up,
        up + 4.5,
        "probe bore",
        tolerance,
    )
    .expect("bore");
    let bored = subtract(&plate, &bore, "probe bore cut", tolerance).expect("bore cut");
    let threaded = cut_thread(
        &bored,
        &spec,
        &ridge,
        thread_form::internal_phase_start(-spec.lead(), &spec),
        4.5 - 0.75 * spec.lead(),
        None,
        None,
        "probe thread",
        tolerance,
    )
    .expect("thread cut");
    // The angled fixture's probed tilt, forced onto a flat plate:
    // same plunge depth, same interior-wrap crossings, no tilt.
    let forced = EntryProbe {
        tilt_tan: 0.4031,
        hull_height: 0.0,
    };
    let r_mouth = spec.major_diameter() / 2.0 + spec.pitch() / 8.0 + clearance;
    let r_in = r_bore - 2.0 * thread_form::embed_margin(&spec);
    let eps = 0.25 * spec.pitch();
    let profile =
        entry_ring_profile(r_in, r_mouth, eps, &forced, "probe chamfer").expect("profile");
    let ring = revolved_tool(&profile, "probe chamfer", tolerance).expect("ring");
    let out = subtract(&threaded, &ring, "probe chamfer cut", tolerance);
    match out {
        Ok(_) => eprintln!("[deep-plunge-probe] perpendicular deep plunge SERVES"),
        Err(e) => panic!("[deep-plunge-probe] perpendicular deep plunge DECLINES: {e}"),
    }
}
