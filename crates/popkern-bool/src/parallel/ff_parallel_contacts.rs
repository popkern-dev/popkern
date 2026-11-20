//! Contact discovery in the face-pair pass (W3a/W3b, tangency design
//! §4.8): each unique surface pair the pass solves and finds meeting in
//! no curve is searched for contacts over the same windows. A touch
//! within the linear reach lying on both trimmed faces is a
//! [`FaceContact`], the vertex table's `SourceKind::Contact` candidate.
//! With the shadow log set, one `ContactDiscovery` line per boolean
//! counts the pairs searched, the contacts, those kept, and the search's
//! time beside the solve's.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::contact_classify::{ContactScale, SurfaceContactKind};
use popkern_core::geometry::intersection::contact_discovery::{
    ContactSearch, find_surface_contacts,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::{FacePairTask, PairInputs};

/// Two operand faces touching at an isolated point, neither crossing the
/// other there: found by the face-pair pass where the pair meets in no
/// curve.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FaceContact {
    /// Operand A's face.
    pub face_a: FaceHandle,
    /// Operand B's face.
    pub face_b: FaceHandle,
    /// The touch, on A's face.
    pub point: Point3,
    /// Its parameters on A's face.
    pub uv_a: (f64, f64),
    /// Its parameters on B's face.
    pub uv_b: (f64, f64),
    /// The signed gap there (within the linear reach).
    pub gap: f64,
}

/// What the pass found its faces touch at where they meet in no curve:
/// isolated touches and tangent curves.
pub(crate) type PassContacts = (
    Vec<FaceContact>,
    Vec<popkern_core::geometry::curve::CurveOnEdge>,
);

/// The faces' size: the larger window's diagonal in space.
fn extent(inp: &PairInputs<'_>, tolerance: &Tolerance) -> f64 {
    let diagonal = |s: &SurfaceGeometry, ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64))| {
        use popkern_core::geometry::surface::SurfaceEval as _;
        s.point_at(u0, v0).distance_to(&s.point_at(u1, v1))
    };
    diagonal(inp.surf_a, inp.r1)
        .max(diagonal(inp.surf_b, inp.r2))
        .max(tolerance.linear())
}

/// Whether the pair's NURBS face window is proved farther than the
/// linear reach from the other, analytic, surface: no touch to search
/// for (`self_contact::patches_clear`).
fn proved_apart(inp: &PairInputs<'_>, tolerance: &Tolerance) -> bool {
    use crate::postprocess::self_contact::patches_clear::window_clear_of_surface;
    let pad = tolerance.linear_reach();
    match (inp.surf_a, inp.surf_b) {
        (SurfaceGeometry::Nurbs(n), other) if !matches!(other, SurfaceGeometry::Nurbs(_)) => {
            window_clear_of_surface(n, &inp.r1, other, pad)
        }
        (other, SurfaceGeometry::Nurbs(n)) if !matches!(other, SurfaceGeometry::Nurbs(_)) => {
            window_clear_of_surface(n, &inp.r2, other, pad)
        }
        (
            SurfaceGeometry::Cylinder {
                origin: o1,
                axis: a1,
                radius: r1,
            },
            SurfaceGeometry::Cylinder {
                origin: o2,
                axis: a2,
                radius: r2,
            },
        ) => cylinders_cannot_touch((*o1, *a1, *r1), (*o2, *a2, *r2), pad),
        _ => false,
    }
}

/// Whether two cylinders have no isolated touch: one needs parallel
/// normals, which for skew axes point along their common perpendicular,
/// and the generators there meet only when the axes lie `r1 + r2` or
/// `|r1 - r2|` apart (within `pad`); parallel axes touch along a line or
/// not at all, a tangent curve, never a point (W3c's).
fn cylinders_cannot_touch(
    (o1, a1, r1): (Point3, popkern_core::geometry::Vec3, f64),
    (o2, a2, r2): (Point3, popkern_core::geometry::Vec3, f64),
    pad: f64,
) -> bool {
    let (Some(a1), Some(a2)) = (a1.normalized(), a2.normalized()) else {
        return false;
    };
    let Some(n) = a1.cross(&a2).normalized() else {
        return true;
    };
    let d = (o2 - o1).dot(&n).abs();
    let (r1, r2) = (r1.abs(), r2.abs());
    (d - (r1 + r2)).abs() > pad && (d - (r1 - r2).abs()).abs() > pad
}

/// Search the unique pair `inp` for contacts.
fn search(inp: &PairInputs<'_>, tolerance: &Tolerance) -> ContactSearch {
    find_surface_contacts(
        inp.surf_a,
        inp.r1,
        inp.surf_b,
        inp.r2,
        ContactScale {
            linear: tolerance.linear(),
            extent: extent(inp, tolerance),
        },
    )
}

/// The touches of each task's faces: every unique pair whose solve found
/// no curve is searched (a pair meeting in a curve has its meetings as
/// records, design §4.8), and each of its tasks keeps the isolated
/// touches within the linear reach lying on both of its trimmed faces.
/// Tangent curves are W3c's, never point candidates.
pub(super) fn discover(
    (tasks, unique_of_task): (&[FacePairTask], &[Option<usize>]),
    unique: &[&PairInputs<'_>],
    solved: &[Result<
        popkern_core::geometry::intersection::SurfaceSurfaceIntersectionResult,
        String,
    >],
    (a, b): (&BRep, &BRep),
    (tolerance, solve_ms): (&Tolerance, f64),
) -> PassContacts {
    let clock = popkern_core::time::Stopwatch::start();
    let curveless = |ui: usize| {
        solved[ui]
            .as_ref()
            .is_ok_and(|r| r.intersections.is_empty())
    };
    let searches: Vec<Option<ContactSearch>> = {
        let one = |ui: usize| {
            (curveless(ui) && !proved_apart(unique[ui], tolerance))
                .then(|| search(unique[ui], tolerance))
        };
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            (0..unique.len()).into_par_iter().map(one).collect()
        }
        #[cfg(not(feature = "parallel"))]
        {
            (0..unique.len()).map(one).collect()
        }
    };
    let reach = tolerance.linear_reach();
    let mut out = Vec::new();
    let mut curves = Vec::new();
    for (task, ui) in tasks.iter().zip(unique_of_task) {
        let Some(Some(found)) = ui.map(|ui| &searches[ui]) else {
            continue;
        };
        let (Some(sa), Some(sb)) = (surface(a, task.face_a), surface(b, task.face_b)) else {
            continue;
        };
        if let Some(inp) = ui.map(|ui| unique[ui]) {
            curves.extend(super::ff_parallel_contact_curves::contact_curves(
                &found.contacts,
                (sa, inp.r1, a, task.face_a),
                (sb, inp.r2, b, task.face_b),
                ContactScale {
                    linear: tolerance.linear(),
                    extent: extent(inp, tolerance),
                },
                reach,
            ));
        }
        for c in &found.contacts {
            let on = |brep: &BRep, face: FaceHandle, s: &SurfaceGeometry, uv: (f64, f64)| {
                use popkern_core::geometry::surface::SurfaceEval as _;
                let p = s.point_at(uv.0, uv.1);
                brep.face_interior_holds(face, &p, reach)
                    || brep.point_near_face_boundary(face, &p, reach)
            };
            if c.gap.abs() > reach
                || !matches!(
                    c.classification.kind,
                    SurfaceContactKind::TangentPoint { .. }
                )
                || !on(a, task.face_a, sa, c.uv1)
                || !on(b, task.face_b, sb, c.uv2)
            {
                continue;
            }
            out.push(FaceContact {
                face_a: task.face_a,
                face_b: task.face_b,
                point: c.point,
                uv_a: c.uv1,
                uv_b: c.uv2,
                gap: c.gap,
            });
        }
    }
    if std::env::var_os("POPKERN_COMPLEX_SHADOW_LOG").is_some() && !unique.is_empty() {
        let searched = searches.iter().flatten().count();
        let contacts: usize = searches.iter().flatten().map(|s| s.contacts.len()).sum();
        crate::complex_shadow::append_log_line(&format!(
            "ContactDiscovery pairs {} searched {searched} contacts {contacts} kept {} search_ms {:.2} solve_ms {solve_ms:.2}",
            unique.len(),
            out.len(),
            clock.elapsed_ms().unwrap_or(0.0)
        ));
    }
    (out, curves)
}

fn surface(brep: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}
