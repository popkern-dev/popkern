//! Seam pixel reassignment (PK-503): the second resolve pass that
//! hands a pixel between the display mesh's seam and the true seam to
//! the neighbouring face whose surface its ray actually meets.

use popkern_core::topology::{BRep, FaceHandle};

use super::{SurfaceHit, ZebraSample, locate};

/// Pixel radius searched for a neighbouring face's surface.
const SEAM_REASSIGN_RADIUS: i64 = 2;

/// Second resolve pass: a pixel whose ray missed its claiming face's
/// domain takes the surface of a neighbouring face (within
/// [`SEAM_REASSIGN_RADIUS`] pixels) that its ray DOES meet inside that
/// face's domain, seeded from the neighbour pixel's located parameter.
/// Without this, those pixels carry the claiming face's boundary
/// normal — a flat strip beside every seam that seam readings misread
/// as a G1 break (measured: 2e-2 rad on pentagon5_blunt cap seams
/// whose analytic G1 residual is 7e-6).
pub(super) fn reassign_seam_pixels(
    brep: &BRep,
    width: u32,
    resolved: &[Option<(ZebraSample, f64)>],
    eye: [f64; 3],
    sign_of: &dyn Fn(FaceHandle) -> f64,
) -> Vec<Option<ZebraSample>> {
    let w = width as i64;
    let h = resolved.len() as i64 / w.max(1);
    let mut out: Vec<Option<ZebraSample>> = resolved
        .iter()
        .map(|r| r.as_ref().map(|(s, _)| *s))
        .collect();
    for (idx, entry) in resolved.iter().enumerate() {
        let Some((sample, reach)) = entry else {
            continue;
        };
        if matches!(sample.hit, SurfaceHit::Ray | SurfaceHit::NoSurface) {
            continue;
        }
        let (x, y) = (idx as i64 % w, idx as i64 / w);
        let mut tried: Vec<FaceHandle> = vec![sample.face];
        let mut neighbours: Vec<(i64, ZebraSample)> = Vec::new();
        for dy in -SEAM_REASSIGN_RADIUS..=SEAM_REASSIGN_RADIUS {
            for dx in -SEAM_REASSIGN_RADIUS..=SEAM_REASSIGN_RADIUS {
                let (nx, ny) = (x + dx, y + dy);
                if nx < 0 || ny < 0 || nx >= w || ny >= h {
                    continue;
                }
                if let Some((n, _)) = &resolved[(ny * w + nx) as usize]
                    && n.hit == SurfaceHit::Ray
                    && n.face != sample.face
                {
                    neighbours.push((dx * dx + dy * dy, *n));
                }
            }
        }
        neighbours.sort_by_key(|(d2, _)| *d2);
        for (_, n) in neighbours {
            if tried.contains(&n.face) {
                continue;
            }
            tried.push(n.face);
            let Some(surface) = brep
                .faces
                .get(n.face)
                .and_then(|f| f.surface.as_ref())
                .map(|sof| &sof.geometry)
            else {
                continue;
            };
            let sense = sign_of(n.face);
            let query = locate::PixelQuery {
                eye,
                facet: sample.position,
                reach: *reach,
                sense,
            };
            let (uv, hit) = locate::locate_on_surface(surface, &[n.uv], &query);
            if hit != SurfaceHit::Ray {
                continue;
            }
            let Some(nrm) = popkern_mesh::surface_limit_normal(surface, uv.0, uv.1) else {
                continue;
            };
            out[idx] = Some(ZebraSample {
                surface_face: n.face,
                uv,
                hit: SurfaceHit::NeighbourRay,
                surface_normal: Some([sense * nrm.x, sense * nrm.y, sense * nrm.z]),
                ..*sample
            });
            break;
        }
    }
    out
}
