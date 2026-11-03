//! Each region of a partition as its own solid, with the faces
//! neighbouring regions share (design 6: conformal meshing needs them).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};

use super::{BodySet, Partition, SetExpr};
use crate::complex::model::SideIndex;
use crate::error::BoolError;
use crate::options::NonManifoldPolicy;

/// One region as a solid.
#[derive(Debug, Clone)]
pub struct RegionSolid {
    /// The solid.
    pub body: BRep,
    /// The bodies holding the region.
    pub bodies: BodySet,
    /// Its index in [`Partition::regions`].
    pub region: usize,
}

/// Two region solids sharing a face: one face in each, on one surface,
/// facing opposite ways.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interface {
    /// Indices into [`RegionSolids::solids`].
    pub solids: (usize, usize),
    /// The shared face in each.
    pub faces: (FaceHandle, FaceHandle),
}

/// Region solids and the faces they share.
#[derive(Debug, Clone)]
pub struct RegionSolids {
    /// One solid per kept region.
    pub solids: Vec<RegionSolid>,
    /// Every face two kept regions share.
    pub interfaces: Vec<Interface>,
}

impl Partition {
    /// Every region `expr` keeps as its own solid, with an [`Interface`]
    /// per face two of them share. A region pinched to itself at an edge
    /// or vertex is served as manifold pieces meeting at recorded
    /// contacts.
    ///
    /// # Errors
    ///
    /// [`BoolError::DegenerateInput`] when `expr` names a body past
    /// [`Partition::body_count`]; extraction's errors per region.
    pub fn region_solids(&self, expr: &SetExpr) -> Result<RegionSolids, BoolError> {
        if let Some(i) = expr.max_body()
            && i >= self.body_count
        {
            return Err(BoolError::DegenerateInput(format!(
                "the expression names body {i}, but the partition holds {}",
                self.body_count
            )));
        }
        let lift = &self.built.complex.lift;
        let c = &lift.complex;
        let mut solids = Vec::new();
        // Kept region index (into `solids`) per complex shell.
        let mut solid_of_shell = vec![None; c.shells.len()];
        for (region, (info, &shell)) in self
            .built
            .regions
            .iter()
            .zip(&self.built.region_shells)
            .enumerate()
        {
            if !expr.holds(&info.bodies, self.body_count) {
                continue;
            }
            let mut only = vec![false; c.shells.len()];
            only[shell] = true;
            let selection =
                crate::complex::select::select_shells(lift, only, self.built.fragments.len());
            let body = self.extract(&selection, NonManifoldPolicy::SplitToManifold, false)?;
            solid_of_shell[shell] = Some(solids.len());
            solids.push(RegionSolid {
                body,
                bodies: info.bodies,
                region,
            });
        }
        let operands: Vec<&BRep> = self.bodies.iter().collect();
        let linear = self.tolerance.linear();
        let mut interfaces = Vec::new();
        for (fi, face) in c.faces.iter().enumerate() {
            let side = |s: SideIndex| {
                face.sides[s as usize]
                    .shell
                    .and_then(|sh| solid_of_shell[sh.index()])
            };
            let (Some(i), Some(j)) = (side(SideIndex::Ahead), side(SideIndex::Behind)) else {
                continue;
            };
            if i == j {
                continue;
            }
            let Some(at) = crate::complex::same_domain_check::inner_point(
                c,
                &operands,
                linear,
                crate::complex::ids::FaceId::from_index(fi),
            ) else {
                continue;
            };
            let (Some(fa), Some(fb)) = (
                face_at(&solids[i].body, &at.point, linear),
                face_at(&solids[j].body, &at.point, linear),
            ) else {
                continue;
            };
            interfaces.push(Interface {
                solids: (i.min(j), i.max(j)),
                faces: if i < j { (fa, fb) } else { (fb, fa) },
            });
        }
        Ok(RegionSolids { solids, interfaces })
    }
}

/// The face of `body` that `p` lies on: within ten tolerances of its
/// surface and inside its trim.
fn face_at(body: &BRep, p: &Point3, linear: f64) -> Option<FaceHandle> {
    body.faces.iter().find_map(|(fh, f)| {
        let surf = f.surface.as_ref()?;
        let (u, v) = surf.geometry.closest_point(p)?;
        if surf.geometry.point_at(u, v).distance_to(p) > linear * 10.0 {
            return None;
        }
        crate::classify_uv_loops::point_in_face_region_exact(body, fh, u, v)
            .unwrap_or(false)
            .then_some(fh)
    })
}
