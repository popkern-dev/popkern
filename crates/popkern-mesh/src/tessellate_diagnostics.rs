//! Env-gated tessellation diagnostics (PK_MESH_DEEPCHORD_TRACE):
//! per face mesh, flag triangle edges whose midpoints deviate from
//! the face surface by far more than the chord tolerance — the
//! honest "deep chord" signal that caught the PK-146 seam-band
//! cocircular fans. Split out of `tessellate.rs` (500-line cap).

use popkern_core::geometry::Point3;
use popkern_core::geometry::SurfaceEval;
use popkern_core::topology::BRep;

use crate::TessellationOptions;

pub(super) fn deepchord_trace(
    brep: &BRep,
    face_meshes: &[crate::mesh::FaceMesh],
    options: &TessellationOptions,
    stage: &str,
) {
    if std::env::var_os("PK_MESH_DEEPCHORD_TRACE").is_some() {
        for fm in face_meshes {
            let Some(sof) = brep.faces.get(fm.face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let mut deep = 0usize;
            let mut worst = 0.0_f64;
            let mut example = None;
            for t in &fm.triangles {
                for (i, j) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    let a = fm.vertices[i].position;
                    let b = fm.vertices[j].position;
                    let mid = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
                    let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                        &sof.geometry,
                        &mid,
                    ) else {
                        continue;
                    };
                    let dev = sof.geometry.point_at(u, v).distance_to(&mid);
                    if dev > 10.0 * options.chord_height {
                        deep += 1;
                        if dev > worst {
                            worst = dev;
                            example = Some((a, b));
                        }
                    }
                }
            }
            if deep > 0 {
                let ex = example
                    .map(|(a, b)| {
                        format!(
                            "({:.3},{:.3},{:.3})-({:.3},{:.3},{:.3})",
                            a.x, a.y, a.z, b.x, b.y, b.z
                        )
                    })
                    .unwrap_or_default();
                eprintln!(
                    "  [deepchord {stage}] face {:?}: {deep} edges dev>{:.1e} worst {worst:.4} ex {ex}",
                    fm.face,
                    10.0 * options.chord_height
                );
            }
        }
    }
}

pub(crate) fn face_mesh_diagonal(fm: &crate::mesh::FaceMesh) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for v in &fm.vertices {
        for (k, c) in [v.position.x, v.position.y, v.position.z]
            .into_iter()
            .enumerate()
        {
            min[k] = min[k].min(c);
            max[k] = max[k].max(c);
        }
    }
    ((max[0] - min[0]).powi(2) + (max[1] - min[1]).powi(2) + (max[2] - min[2]).powi(2)).sqrt()
}
