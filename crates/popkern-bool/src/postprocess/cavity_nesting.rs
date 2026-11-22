//! Cavities wired as solids of their own, nested once the result is
//! complete (corpus SS-X3).
//!
//! The cellular assembly splits its result's faces into shells and nests
//! each enclosed component as a cavity of the solid holding it
//! (`cellular::assemble_cleanup::partition_shells`), deciding by ray
//! parity before the faces' trims and orientation are final. A sphere
//! piece bounded by one circle takes its side from its loop's winding
//! against its orientation, so mid-assembly a point in a void cut from two
//! crossing balls by a third read Outside them: the void was wired as a
//! second solid, its shell enclosing negative volume. Earlier still that
//! void was dropped whole and N-ary (A ∪ B) − C served A ∪ B,
//! validator-clean. Here, after the output repairs, a solid region whose
//! outer shell encloses negative volume (the sign of the sampled volume,
//! as the validator reads it) is moved into the smallest other solid whose
//! outer shell holds a point of its void by the boolean's ray parity, and
//! its emptied region is removed; one nothing holds stays for the
//! validator to report. The sampled winding and volume cannot decide it:
//! they undercount sphere pieces bounded by short arcs (a void of 0.66
//! read 0.50, its host of 19.9 read 10.8, the host's winding 0.5 inside).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, RegionHandle, RegionKind, ShellHandle};

use crate::classify::point_in_solid_final;
use crate::error::PointClassification;

/// One solid region read by its outer shell.
struct SolidCell {
    region: RegionHandle,
    shell: ShellHandle,
    faces: Vec<FaceHandle>,
    volume: f64,
}

/// Nest every solid region whose outer shell encloses negative volume in
/// the smallest other solid holding its void. Returns how many moved.
pub(crate) fn nest_inverted_solid_shells(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let cells = solid_cells(brep);
    if !cells.iter().any(|c| c.volume < 0.0) {
        return 0;
    }
    let boxes = crate::intersect::compute_face_aabbs(brep);
    let boxes_of = |faces: &[FaceHandle]| -> Vec<_> {
        boxes
            .iter()
            .filter(|(f, _)| faces.contains(f))
            .cloned()
            .collect()
    };
    let inside = |q: &Point3, faces: &[FaceHandle]| {
        matches!(
            point_in_solid_final(q, brep, &boxes_of(faces), tolerance),
            Ok(PointClassification::Inside)
        )
    };
    let mut moves: Vec<(ShellHandle, RegionHandle, RegionHandle)> = Vec::new();
    for cell in cells.iter().filter(|c| c.volume < 0.0) {
        let Some(q) = void_point(brep, &cell.faces, &boxes_of(&cell.faces), &|q| {
            inside(q, &cell.faces)
        }) else {
            continue;
        };
        let host = cells
            .iter()
            .filter(|h| h.volume > 0.0 && h.region != cell.region)
            .filter(|h| !h.faces.iter().any(|f| cell.faces.contains(f)))
            .filter(|h| inside(&q, &h.faces))
            .min_by(|a, b| a.volume.total_cmp(&b.volume));
        if let Some(h) = host {
            moves.push((cell.shell, cell.region, h.region));
        }
    }
    for &(shell, from, to) in &moves {
        brep.shells[shell].region = to;
        brep.regions.remove(from);
    }
    moves.len()
}

/// Every solid region with its outer shell's faces and the signed volume
/// they enclose.
fn solid_cells(brep: &BRep) -> Vec<SolidCell> {
    brep.regions
        .iter()
        .filter(|&(r, data)| r != brep.outside && data.kind == RegionKind::Solid)
        .filter_map(|(region, data)| {
            let shell = data.shell;
            let sides: Vec<(FaceHandle, f64)> = brep
                .faces
                .iter()
                .filter_map(|(h, f)| match (f.back == shell, f.front == shell) {
                    (true, false) => Some((h, 1.0)),
                    (false, true) => Some((h, -1.0)),
                    _ => None,
                })
                .collect();
            let volume = brep.enclosed_signed_volume(&sides)?;
            Some(SolidCell {
                region,
                shell,
                faces: sides.iter().map(|s| s.0).collect(),
                volume,
            })
        })
        .collect()
}

/// A point of the void `faces` bound, read inside them by `inside`
/// (parity, blind to their orientation): marched from each face's sample
/// point along its normal, either way, in tenths of the faces' box.
fn void_point(
    brep: &BRep,
    faces: &[FaceHandle],
    boxes: &[(
        FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    inside: &dyn Fn(&Point3) -> bool,
) -> Option<Point3> {
    let (lo, hi) = boxes
        .iter()
        .map(|(_, b)| (b.min, b.max))
        .reduce(|(a0, a1), (b0, b1)| {
            (
                Point3::new(a0.x.min(b0.x), a0.y.min(b0.y), a0.z.min(b0.z)),
                Point3::new(a1.x.max(b1.x), a1.y.max(b1.y), a1.z.max(b1.z)),
            )
        })?;
    let size = lo.distance_to(&hi);
    faces.iter().find_map(|&face| {
        let p = crate::classify::pick_sample_point(brep, face);
        let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
        let (u, v) = surface.closest_point(&p)?;
        let n = surface.normal_at(u, v).normalized()?;
        [1.0, -1.0].into_iter().find_map(|sign| {
            (1..5)
                .map(|k| p + n * (sign * 0.1 * k as f64 * size))
                .find(|q| inside(q))
        })
    })
}
