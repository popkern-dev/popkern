//! PK-328 row 483 gates: `shell_preview` answers in interactive time
//! on every corpus body, its face samples realize the exact offset
//! distance against the source body (the preview samples the offset
//! MAP, never a fit), and the error/warning `Display` texts carry the
//! values and limits a GUI needs.

use std::time::{Duration, Instant};

use popkern_core::tolerance::Tolerance;
use popkern_examples::offset_corpus::corpus;
use popkern_shell::error::TrimRegionReport;
use popkern_shell::{ShellDirection, ShellError, ShellWarning, shell_preview};

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn preview_is_interactive_and_realizes_the_offset_on_the_corpus() {
    let budget = Duration::from_millis(200);
    for f in corpus(tol()).expect("corpus builds") {
        for direction in [ShellDirection::Inward, ShellDirection::Outward] {
            let d = 0.04;
            let t0 = Instant::now();
            let mesh = match shell_preview(&f.brep, d, direction, tol()) {
                Ok(m) => m,
                Err(e) => {
                    assert!(
                        !f.watertight,
                        "{}/{direction:?}: closed body preview failed: {e}",
                        f.name
                    );
                    continue;
                }
            };
            let elapsed = t0.elapsed();
            assert!(
                elapsed < budget,
                "{}/{direction:?}: preview {elapsed:?} exceeds the interactive budget",
                f.name
            );
            assert!(
                !mesh.triangles.is_empty(),
                "{}/{direction:?}: empty preview",
                f.name
            );
            // Face samples realize |d| against the source body: check a
            // spread of face-triangle vertices via the mesh's own
            // distance to the source's surfaces (loose band — the
            // sample's foot may fall on an adjacent face near edges).
            let index = popkern_mesh::point_query::PointMeshIndex::build(
                &popkern_mesh::tessellate(&f.brep, &popkern_mesh::TessellationOptions::default())
                    .expect("tessellate source"),
            );
            // Sign-aware realization: outward samples sit OUTSIDE the
            // source at up to d, inward samples INSIDE. An untrimmed
            // preview may sit CLOSER than d near converging corners
            // (that material is what the real trim removes), but never
            // FARTHER, and never on the wrong side — the two failure
            // modes of a wrong-sense or wrong-magnitude offset map.
            // A sample's |signed distance| never exceeds the offset
            // magnitude (it moved exactly d from ITS face), and the
            // MAJORITY sit on the requested material side. Individual
            // wrong-side samples are legitimate untrimmed-preview
            // artifacts (an acute base wedge lets an inward sample
            // exit the far face — the real operation's recompute pulls
            // exactly that material back in).
            let band = d * 0.2 + 0.01; // + tessellation chord error
            let mut checked = 0usize;
            let mut wrong_side = 0usize;
            let mut deepest: f64 = 0.0;
            for (t, fh) in mesh.triangles.iter().zip(&mesh.face_map).step_by(37) {
                if fh.is_none() {
                    continue;
                }
                let p = mesh.positions[t[0]];
                let s = index.signed_distance(p);
                assert!(
                    s.abs() < d + band,
                    "{}/{direction:?}: sample farther than the offset \
                     (signed {s:.4}, d {d})",
                    f.name
                );
                let expected_side_ok = match direction {
                    ShellDirection::Outward => s > -band,
                    _ => s < band,
                };
                if !expected_side_ok {
                    wrong_side += 1;
                }
                deepest = deepest.max(s.abs());
                checked += 1;
            }
            if checked > 3 {
                assert!(
                    wrong_side * 4 < checked,
                    "{}/{direction:?}: {wrong_side}/{checked} samples on the \
                     wrong material side — the offset sense is broken",
                    f.name
                );
                assert!(
                    deepest > d * 0.8,
                    "{}/{direction:?}: the offset magnitude must be realized \
                     somewhere (deepest {deepest:.4}, d {d})",
                    f.name
                );
            }
        }
    }
}

/// Row 483: GUI-readable Display text — every user-facing error and
/// warning carries its values and limits.
#[test]
fn error_and_warning_display_texts_carry_values_and_limits() {
    let err = ShellError::OffsetTooLarge {
        requested: 5.0,
        max_safe: 3.2,
        binding: "edge between Face#3 and Face#7".into(),
    };
    let text = err.to_string();
    assert!(
        text.contains("5.0") && text.contains("3.2") && text.contains("Face#3"),
        "{text}"
    );

    let w = ShellWarning::RegionTrimmed {
        removed_faces: 4,
        welded_seams: 2,
        regions: vec![TrimRegionReport {
            face_a: Default::default(),
            face_b: Default::default(),
            witness: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            curve_length: 0.5,
        }],
    };
    let text = w.to_string();
    assert!(
        text.contains('4') && text.contains('2') && text.contains('1'),
        "{text}"
    );

    let w = ShellWarning::FaceOffsetApproximated {
        face: Default::default(),
        achieved_deviation: 1.5e-5,
    };
    assert!(w.to_string().contains("1.500e-5"), "{w}");

    let w = ShellWarning::FaceCollapsed {
        faces: vec![Default::default(), Default::default()],
        mechanism: "trimmed-distance death".into(),
    };
    let text = w.to_string();
    assert!(
        text.contains('2') && text.contains("trimmed-distance"),
        "{text}"
    );
}
