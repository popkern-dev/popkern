use super::*;
use crate::mesh_test_helpers::{face_areas, icosphere, model_of};

/// A one-face sphere model whose kernel area is `exact` (`None` keeps
/// the face freeform).
fn sphere(exact: Option<f64>) -> crate::model::AnalysisModel {
    let mut model = model_of(icosphere([0.0; 3], 1.0, 3), 1);
    if exact.is_some() {
        model.faces[0].surface = Surface::Sphere {
            centre: [0.0; 3],
            radius: 1.0,
        };
    }
    model.faces[0].exact_area = exact;
    model
}

#[test]
fn a_face_tessellated_beyond_its_kernel_area_is_mismatched() {
    let sphere_area = 4.0 * std::f64::consts::PI;
    let tessellated = face_areas(&sphere(None).solids[0], 1)[0];
    // The whole sphere against its own area (chords make it a little
    // smaller): matched.
    let whole = sphere(Some(sphere_area));
    assert!(mismatched_faces(&whole.solids[0], &whole.faces).is_empty());
    // Against an eighth of it, as a corner patch tessellated as the
    // whole sphere (PK-787): mismatched, and every triangle flagged.
    let corner = sphere(Some(sphere_area / 8.0));
    let bad = mismatched_faces(&corner.solids[0], &corner.faces);
    assert_eq!(bad, vec![(0, tessellated, 0.0)]);
    let (flags, warning) = mismatch_warning(&corner.solids[0], &corner.faces);
    assert!(flags.iter().all(|f| *f));
    assert!(matches!(
        warning,
        Some(crate::report::Warning::TessellationMismatch { ref faces, .. }) if faces == &[0]
    ));
    // A kernel area far above the tessellation (a parameter box around
    // a trimmed domain) is only a loose bound: matched.
    let loose = sphere(Some(100.0 * sphere_area));
    assert!(mismatched_faces(&loose.solids[0], &loose.faces).is_empty());
    // Freeform faces are not checked.
    let freeform = sphere(None);
    assert!(
        mismatch_warning(&freeform.solids[0], &freeform.faces)
            .1
            .is_none()
    );
}

#[test]
fn a_face_folded_against_its_surface_is_mismatched() {
    let sphere_area = 4.0 * std::f64::consts::PI;
    let mut model = sphere(Some(sphere_area));
    let tessellated = face_areas(&model.solids[0], 1)[0];
    // Folds under the 1% threshold are chord noise; above it the face is
    // mismatched whatever its area.
    model.faces[0].folded_area = Some(0.005 * tessellated);
    assert!(mismatched_faces(&model.solids[0], &model.faces).is_empty());
    model.faces[0].folded_area = Some(0.02 * tessellated);
    assert_eq!(
        mismatched_faces(&model.solids[0], &model.faces),
        vec![(0, tessellated, 0.02 * tessellated)]
    );
}

/// A coarse sphere's tessellation falls short of 4πr³/3 by more than
/// the threshold and warns with both volumes; a finer one within it
/// does not (the control), nor does a solid without a kernel volume.
#[test]
fn a_tessellation_short_of_the_exact_volume_warns() {
    let exact = 4.0 / 3.0 * std::f64::consts::PI;
    let with_exact = |levels: usize| {
        let mut mesh = icosphere([0.0; 3], 1.0, levels);
        mesh.exact_volume = Some(exact);
        mesh
    };
    let coarse = with_exact(1);
    let short = 1.0 - coarse.enclosed_volume() / exact;
    assert!(short > VOLUME_DEFICIT, "{short}");
    assert!(matches!(
        volume_warning(&coarse),
        Some(crate::report::Warning::TessellationVolume { solid: 0, tessellated, kernel })
            if tessellated == coarse.enclosed_volume() && kernel == exact
    ));
    let fine = with_exact(4);
    assert!(1.0 - fine.enclosed_volume() / exact < VOLUME_DEFICIT);
    assert_eq!(volume_warning(&fine), None);
    assert_eq!(volume_warning(&icosphere([0.0; 3], 1.0, 1)), None);
}
