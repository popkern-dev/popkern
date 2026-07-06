//! Flagship boss plate, written against the `popkern-part` modeling
//! API — the acceptance gate for the API's phase 1 (see
//! `docs/architecture/part-modeling-api-design.md`).
//!
//! Same part, same gates as `flagship_boss_plate`: a rectangular
//! plate and a penetrating cylindrical boss, unioned, with the
//! boss/plate seam and the full 12-edge plate chain filleted; the
//! result must validate, tessellate watertight, and round-trip its
//! STEP export with every face intact.

use popkern_part::prelude::*;

fn build() -> Result<Part, PartError> {
    let plate = Profile::on(Plane::XY).rect(2.4, 1.6).extrude(0.4)?;
    let boss = Profile::on(Plane::xy_at(0.2)).circle(0.35).extrude(0.8)?;

    plate
        .union(&boss)?
        // One blend operation, three selections: the circular
        // boss/plate seam, the four top-perimeter edges, and the
        // four vertical corner edges — blended together so the
        // chain's corner patches resolve in one pass.
        .fillet(0.08, |e| {
            e.between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(0.4)
                .or(|e| e.straight().at_height(0.4))
                .or(|e| e.straight().direction(Dir::Z))
        })
}

fn main() {
    let part = build().expect("boss plate build");
    println!(
        "flagship (API): {} faces, {} edges",
        part.face_count(),
        part.edge_count()
    );
    for feature in part.features().entries() {
        println!("  feature: {} — {}", feature.name, feature.summary);
    }

    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples/output/flagship_boss_plate_api.step");
    let part = part.assert_watertight().expect("watertight gate");
    part.export_step(&out).expect("STEP export");
    let part = part.assert_step_roundtrip(&out).expect("round-trip gate");
    println!(
        "watertight ✓  round-trip ✓  wrote {} ({} faces)",
        out.display(),
        part.face_count()
    );
}
