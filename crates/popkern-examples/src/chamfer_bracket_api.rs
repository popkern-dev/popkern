//! Machined mounting tray via the `popkern-part` API — the chamfer
//! showcase, plus coverage for API surface no other example
//! exercises: boolean SUBTRACT (a penetrating rectangular pocket),
//! the `max_fillet_radius` slider bound, rigid transforms, and
//! OBJ export alongside STEP.
//!
//! The part: a rectangular plate with a deep rectangular pocket
//! (bearing-seat style) and the four outer vertical corners
//! chamfered. Gates: kernel validation at every step (built into
//! the API), watertight tessellation, STEP round-trip with full
//! face census, STEP + OBJ exports.
//!
//! Scope notes, honestly recorded: interior THROUGH-holes and
//! side-overlap notches are open boolean/tessellation frontiers
//! (ledger PK-100), and blend-chain corners (fillet-fillet at 135°,
//! bevel-bevel at box corners) are open blend frontiers — this part
//! uses the configurations the kernel closes watertight today.
//! Fillet coverage lives in the flagship and pentagon examples.

use popkern_part::prelude::*;

const PLATE_W: f64 = 3.0;
const PLATE_D: f64 = 1.8;
const PLATE_T: f64 = 0.5;
const POCKET_W: f64 = 2.0;
const POCKET_D: f64 = 1.0;
const POCKET_FLOOR: f64 = 0.2;

fn build() -> Result<Part, PartError> {
    let plate = Profile::on(Plane::XY)
        .rect(PLATE_W, PLATE_D)
        .extrude(PLATE_T)?;

    // Rectangular pocket: a penetrating tool from above (proud of
    // the top face so the boolean sees clean transversal crossings).
    let tool = Profile::on(Plane::xy_at(POCKET_FLOOR))
        .rect(POCKET_W, POCKET_D)
        .extrude(PLATE_T)?;
    let pocketed = plate.subtract(&tool)?;

    // Chamfer the four OUTER vertical corner edges (full plate
    // height distinguishes them from the pocket's verticals). The
    // blend-size slider bound comes from the API.
    let outer_vertical = |e: EdgeQuery| {
        e.straight()
            .direction(Dir::Z)
            .chord_length_between(0.4, 1.0)
    };
    let rmax = pocketed.max_fillet_radius(outer_vertical).unwrap_or(0.25);
    println!("outer-vertical blend size bound (API slider bound): {rmax:.3}");
    let chamfered = pocketed.chamfer(0.18_f64.min(rmax), outer_vertical)?;

    // Rigid placement (transform coverage): lift the tray to a
    // workspace offset, as an assembly would.
    chamfered.translate(0.0, 0.0, 1.0)
}

fn main() {
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("out dir");

    let part = build().expect("tray build");
    println!(
        "tray: {} faces, {} edges",
        part.face_count(),
        part.edge_count()
    );
    for feature in part.features().entries() {
        println!("  {} — {}", feature.name, feature.summary);
    }

    let part = part.assert_watertight().expect("watertight gate");
    let step = out_dir.join("chamfer_bracket_api.step");
    let obj = out_dir.join("chamfer_bracket_api.obj");
    part.export_step(&step).expect("STEP export");
    part.export_obj(&obj).expect("OBJ export");
    let part = part.assert_step_roundtrip(&step).expect("round-trip gate");
    println!(
        "watertight ✓  round-trip ✓  wrote {} and {} ({} faces)",
        step.display(),
        obj.display(),
        part.face_count()
    );
}
