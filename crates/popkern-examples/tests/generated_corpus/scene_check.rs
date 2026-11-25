//! The scene's independent check (design §6.1; split from `scene.rs` for
//! the file-size limit).

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;

use super::geometry::{closest_approach, closest_approach_except, inside};
use super::scene::{Config, Reject, Scene, place};

/// Probe points of a mesh: every vertex and every triangle centroid,
/// thinned evenly to at most `cap` (a box mesh has eight vertices: a
/// fixed stride would test one of them).
fn probes(mesh: &popkern_mesh::TriMesh, cap: usize) -> Vec<Point3> {
    let mut all: Vec<Point3> = mesh.vertices.iter().map(|v| v.position).collect();
    all.extend(mesh.triangles.iter().map(|t| {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        Point3::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        )
    }));
    let step = all.len().div_ceil(cap).max(1);
    all.into_iter().step_by(step).collect()
}

/// The outward normal of `mesh` at its triangle nearest `p` (by centroid;
/// the meshes wind outward, `generated_corpus_selftest`).
fn outward_normal_near(
    mesh: &popkern_mesh::TriMesh,
    p: Point3,
) -> Option<popkern_core::geometry::Vec3> {
    let corner = |t: &[usize; 3], k: usize| mesh.vertices[t[k]].position;
    let t = mesh.triangles.iter().min_by(|x, y| {
        let d = |t: &[usize; 3]| {
            let (a, b, c) = (corner(t, 0), corner(t, 1), corner(t, 2));
            Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            )
            .distance_to(&p)
        };
        d(x).total_cmp(&d(y))
    })?;
    (corner(t, 1) - corner(t, 0))
        .cross(&(corner(t, 2) - corner(t, 0)))
        .normalized()
}

/// The outward normal of `body`'s faces at `c`: the folded chart normal
/// times the body's material sense, of every face passing within `reach`
/// of `c` there inside its trim. `None` unless they agree (a contact on
/// an edge or a pole decides nothing).
fn exact_outward_normal(
    body: &popkern_core::topology::BRep,
    c: Point3,
    reach: f64,
) -> Option<popkern_core::geometry::Vec3> {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let mut found: Option<popkern_core::geometry::Vec3> = None;
    for (fh, f) in body.faces.iter() {
        let Some(s) = f.surface.as_ref() else {
            continue;
        };
        let Some((u, v)) = s.geometry.closest_point(&c) else {
            continue;
        };
        if s.geometry.point_at(u, v).distance_to(&c) > reach
            || matches!(
                body.face_uv_trim_containment(fh, &s.geometry, u, v),
                None | Some(popkern_core::validation::FaceUvContainment::Outside)
            )
        {
            continue;
        }
        let sense = body.component_material_sense(fh)?.outward_multiplier();
        let n = s.folded_normal(u, v)? * sense;
        match found {
            Some(m) if m.dot(&n) < 0.99 => return None,
            _ => found = Some(n),
        }
    }
    found
}

/// The scene's independent check (design §6.1), on the placed meshes.
pub(crate) fn check(scene: &Scene) -> Result<(), Reject> {
    let tol = Tolerance::default().linear();
    let placed = place(scene);
    let chord = (placed.a_op.chord + placed.b_op.chord * scene.b_scale) * scene.scale;
    let margin = 100.0 * tol + 2.0 * chord;
    match scene.cell.config {
        Config::E | Config::I | Config::V => {
            let c = scene
                .contact
                .map(|p| Point3::new(p[0], p[1], p[2]))
                .expect("contact");
            // No second approach outside a ball about the contact; the
            // ball scales with the operands so a gentle touch's own
            // neighbourhood is not read as a second approach.
            let ball = 0.15 * placed.b_op.size * scene.b_scale * scene.scale;
            if scene.cell.config != Config::V
                && let Some(d) =
                    closest_approach(&placed.a_mesh, &placed.b_mesh, margin, Some((c, ball)))
            {
                return Err(Reject::SecondApproach(d));
            }
            // An approach exists at the contact.
            let near = closest_approach(
                &placed.a_mesh,
                &placed.b_mesh,
                scene.offset.abs() * tol + 4.0 * chord + 1e-12,
                None,
            );
            if near.is_none() {
                return Err(Reject::NoApproach(scene.offset));
            }
            // B lies on its configuration's side of A away from the
            // contact: inside for an internal touch, outside for an
            // external one (a steered contact's normal once came out
            // reversed and placed an internal touch outside A).
            // A tangent crossing's flag: along A's outward normal at the
            // contact, B's material lies inward (normals aligned, B inside
            // A there) or outward (opposed). Unchecked, a scene drawn by a
            // generator whose sampled normal came out reversed kept
            // "opposed" over an internal touch of two cylinders, and its
            // judge expected the union to pinch where the differences do.
            // The operands' exact outward normals at the contact decide
            // first: tangent there, B's material lies inward along A's
            // normal when they agree (the probes below left a cone x cone
            // crossing unchecked whose normals agree, flagged opposed,
            // and its judge expected the union to pinch, CC-V 87).
            let reach = 10.0 * (scene.offset.abs() + 1.0) * tol;
            let exact = exact_outward_normal(&placed.a, c, reach)
                .zip(exact_outward_normal(&placed.b, c, reach))
                .map(|(na, nb)| na.dot(&nb) > 0.0);
            if scene.cell.config == Config::V && exact.is_some_and(|al| al != scene.aligned) {
                return Err(Reject::WrongSide);
            }
            if scene.cell.config == Config::V && exact.is_none() {
                // Shorter probes where both read alike (the contact by B's
                // rim, a probe leaving it through a cap); none deciding
                // leaves the flag unchecked.
                let n = outward_normal_near(&placed.a_mesh, c).ok_or(Reject::WrongSide)?;
                let floor = 10.0 * scene.offset.abs() * tol;
                let side = [20.0, 10.0, 5.0].into_iter().find_map(|k| {
                    let s = (k * chord + floor).min(0.2 * ball);
                    let below = inside(&placed.b_mesh, c - n * s);
                    (below != inside(&placed.b_mesh, c + n * s)).then_some(below)
                });
                if side.is_some_and(|below| below != scene.aligned) {
                    return Err(Reject::WrongSide);
                }
            }
            if scene.cell.config != Config::V {
                let away: Vec<Point3> = probes(&placed.b_mesh, 400)
                    .into_iter()
                    .filter(|p| p.distance_to(&c) > ball)
                    .collect();
                let want = scene.cell.config == Config::I;
                if away.iter().any(|&p| inside(&placed.a_mesh, p) != want) {
                    return Err(Reject::WrongSide);
                }
                // An external touch has neither inside the other: B's
                // points outside A also hold when B holds A touching it
                // inside (a steered sphere of radius 1.15 about one of
                // 0.9, centres 0.25 apart, scene 696 of SS-E).
                if scene.cell.config == Config::E
                    && probes(&placed.a_mesh, 400)
                        .into_iter()
                        .filter(|p| p.distance_to(&c) > ball)
                        .any(|p| inside(&placed.b_mesh, p))
                {
                    return Err(Reject::WrongSide);
                }
            }
        }
        Config::D => {
            if closest_approach(&placed.a_mesh, &placed.b_mesh, margin, None).is_some() {
                return Err(Reject::NotDisjoint);
            }
        }
        Config::N => {
            if closest_approach(&placed.a_mesh, &placed.b_mesh, margin, None).is_some()
                || !probes(&placed.b_mesh, 400)
                    .into_iter()
                    .all(|p| inside(&placed.a_mesh, p))
            {
                return Err(Reject::NotContained);
            }
        }
        Config::O => {}
        Config::K => {
            let curve: Vec<Point3> = scene
                .curve
                .iter()
                .map(|p| Point3::new(p[0], p[1], p[2]))
                .collect();
            let tube = 0.12 * placed.b_op.size * scene.b_scale * scene.scale;
            let balls: Vec<(Point3, f64)> = curve.iter().map(|p| (*p, tube)).collect();
            if let Some(d) = closest_approach_except(&placed.a_mesh, &placed.b_mesh, margin, &balls)
            {
                return Err(Reject::SecondApproach(d));
            }
            // An approach at every sample of the curve.
            let reach = scene.offset.abs() * tol + 4.0 * chord + 1e-12;
            let grid = super::geometry::Grid::new(
                &placed.b_mesh,
                (2.0 * reach).max(placed.b_op.chord * 10.0),
            );
            if curve
                .iter()
                .any(|p| grid.distance_within(*p, reach + 2.0 * chord).is_none())
            {
                return Err(Reject::NoApproach(scene.offset));
            }
        }
        Config::X | Config::X3 => {
            let sample: Vec<bool> = probes(&placed.b_mesh, 400)
                .into_iter()
                .map(|p| inside(&placed.a_mesh, p))
                .collect();
            if !(sample.iter().any(|x| *x) && sample.iter().any(|x| !*x)) {
                return Err(Reject::NotOverlapping);
            }
        }
    }
    Ok(())
}
