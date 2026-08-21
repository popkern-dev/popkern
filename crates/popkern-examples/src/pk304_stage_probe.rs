//! PK-304 stage probe: per-stage entity counts for the torus booleans.
//! Transient diagnostic bin; removed at campaign close.

use popkern_bool::progress::{ProgressEvent, ProgressReporter};
use popkern_bool::{BoolOperation, boolean_with_progress};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

struct StderrReporter;
impl ProgressReporter for StderrReporter {
    fn on_phase(&mut self, event: ProgressEvent) {
        eprintln!(
            "  stage {:?}: {:?} ({} entities)",
            event.stage, event.duration, event.entity_count
        );
    }
}

fn main() {
    let tol = Tolerance::default();
    let torus = BRep::torus(1.0, 0.35, tol).expect("torus builds");
    let cyl = BRep::cylinder(0.5, 3.0, tol)
        .expect("cylinder builds")
        .transform(&Transform3::translation(0.8, 0.0, -1.5))
        .expect("translate");
    eprintln!("=== torus x cylinder ===");
    match boolean_with_progress(BoolOperation::Union, &torus, &cyl, tol, &mut StderrReporter) {
        Ok(b) => eprintln!(
            "OK: {} faces validate={:?} validate_all={}",
            b.faces.len(),
            b.validate().err().map(|e| e.to_string()),
            b.validate_all().len()
        ),
        Err(e) => eprintln!("ERR: {e}"),
    }
    let torus_b = BRep::torus(1.0, 0.35, tol)
        .expect("torus builds")
        .transform(&Transform3::rotation_x(0.5).compose(&Transform3::translation(0.7, 0.4, 0.15)))
        .expect("transform");
    eprintln!("=== torus x torus ===");
    match boolean_with_progress(
        BoolOperation::Union,
        &torus,
        &torus_b,
        tol,
        &mut StderrReporter,
    ) {
        Ok(b) => eprintln!(
            "OK: {} faces validate={:?}",
            b.faces.len(),
            b.validate().err().map(|e| e.to_string())
        ),
        Err(e) => eprintln!("ERR: {e}"),
    }
}
