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

/// Two region solids sharing a patch of surface: the face of each that
/// holds it, on one surface, facing opposite ways, and the patch
/// itself. A face may reach past the patch: the solids' extraction
/// merges co-surfacial faces, so a block offset on another's face holds
/// the patch inside a larger face. `patch` is exactly what they share
/// (its area is the contact's, its boundary where a conformal mesh of
/// the two must agree).
#[derive(Debug, Clone)]
pub struct Interface {
    /// Indices into [`RegionSolids::solids`].
    pub solids: (usize, usize),
    /// The face holding the patch in each.
    pub faces: (FaceHandle, FaceHandle),
    /// The shared patch as a one-face sheet, in the first solid's
    /// sense of its holding face's material (facing out of it).
    pub patch: BRep,
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
        for (region, (info, shells)) in self
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
            for &shell in shells {
                only[shell] = true;
            }
            let selection =
                crate::complex::select::select_shells(lift, only, self.built.fragments.len());
            let body = self.extract(&selection, NonManifoldPolicy::SplitToManifold, false)?;
            for &shell in shells {
                solid_of_shell[shell] = Some(solids.len());
            }
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
            let patch = self.face_patch(&face.source_fragments)?;
            interfaces.push(Interface {
                solids: (i.min(j), i.max(j)),
                faces: if i < j { (fa, fb) } else { (fb, fa) },
                patch,
            });
        }
        Ok(RegionSolids { solids, interfaces })
    }
}

impl Partition {
    /// One complex face as a one-face sheet, from its first source
    /// fragment.
    fn face_patch(&self, fragments: &[usize]) -> Result<BRep, BoolError> {
        let mut selections = vec![
            crate::merge::fragment_classify::FragmentSelection {
                keep: false,
                flip: false,
            };
            self.built.fragments.len()
        ];
        if let Some(sel) = fragments.first().and_then(|&f| selections.get_mut(f)) {
            sel.keep = true;
        }
        let operands: Vec<&BRep> = self.bodies.iter().collect();
        let mut carried = crate::serve_region::Carried::default();
        let body = crate::cellular::assemble::cellular_assemble(
            &self.built.graph,
            &self.built.fragments,
            &selections,
            &operands,
            self.tolerance,
            &mut carried,
            Default::default(),
            None,
        )
        .map_err(|d| d.into_error())?;
        super::lamina_part::open_sheet(body)
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
