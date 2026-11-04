//! Regions bounded by more than one shell.
//!
//! The complex's shells are the connected sets of face sides; a region
//! is bounded by one or more of them. Where every face of the complex
//! is reachable from every other across edges, each shell bounds its
//! own region. Where the complex falls apart into separate face
//! components — a body floating inside another, an operand's internal
//! void, a sheet wholly inside a solid — the outside shell of an inner
//! component and the inside of the cell enclosing it bound ONE region:
//! a split must serve the enclosing piece hollow, not whole with an
//! inside-out copy of the void beside it.
//!
//! A shell bounds its region from outside (the region's outer shell)
//! exactly when the region's interior point lies inside the shell's own
//! closed surface: ray parity, which the shell's orientation does not
//! touch. Any other shell — a void's, or a lamina's holding both sides
//! of its faces — is an inner shell of the innermost outer shell of another
//! component whose surface encloses its interior point. The two must
//! agree on membership; a disagreement, or an inner shell nothing
//! encloses in a bounded region, declines.

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use super::Partition;
use crate::complex::model::CellComplex;
use crate::error::{BoolError, PointClassification};
use crate::options::NonManifoldPolicy;

/// The face component of every face: faces sharing an edge share one.
fn face_components(c: &CellComplex) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..c.faces.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for edge in &c.edges {
        let mut faces = edge.radial.iter().map(|f| c.fins[f.index()].face.index());
        let Some(first) = faces.next() else {
            continue;
        };
        for other in faces {
            let (a, b) = (find(&mut parent, first), find(&mut parent, other));
            parent[a] = b;
        }
    }
    (0..c.faces.len()).map(|f| find(&mut parent, f)).collect()
}

fn declined(reason: String) -> BoolError {
    BoolError::MergeFailed {
        stage: "Complex".to_string(),
        reason: format!("partition: {reason}"),
    }
}

/// Whether `p` lies inside `body`'s closed surface.
fn encloses(body: &BRep, p: &Point3, tolerance: &popkern_core::tolerance::Tolerance) -> bool {
    let boxes = crate::intersect::compute_face_aabbs(body);
    matches!(
        crate::classify::point_in_solid_final(p, body, &boxes, tolerance),
        Ok(PointClassification::Inside)
    )
}

impl Partition {
    /// Gather each region's shells (module doc). A complex of one face
    /// component is left as built.
    pub(super) fn nest_regions(&mut self) -> Result<(), BoolError> {
        let c = &self.built.complex.lift.complex;
        let component = face_components(c);
        let shell_component =
            |s: usize| c.shells[s].sides.first().map(|(f, _)| component[f.index()]);
        let mut distinct: Vec<usize> = component.clone();
        distinct.sort_unstable();
        distinct.dedup();
        if distinct.len() < 2 {
            return Ok(());
        }
        // Each region's shell surface, when it bounds the region from
        // outside; `None` for an inner shell.
        let regions = self.built.regions.len();
        let mut outer: Vec<Option<BRep>> = Vec::with_capacity(regions);
        for r in 0..regions {
            let s = self.built.region_shells[r][0];
            // A shell holding both sides of every face it has (a
            // lamina meeting nothing) encloses nothing.
            let two_sided = c.shells[s].sides.iter().all(|&(f, _)| {
                let sides = &c.faces[f.index()].sides;
                sides[0].shell == sides[1].shell
            });
            if two_sided {
                outer.push(None);
                continue;
            }
            let mut only = vec![false; c.shells.len()];
            only[s] = true;
            let selection = crate::complex::select::select_shells(
                &self.built.complex.lift,
                only,
                self.built.fragments.len(),
            );
            let surface = self.extract(&selection, NonManifoldPolicy::SplitToManifold, false)?;
            let inside = encloses(
                &surface,
                &self.built.regions[r].interior_point,
                &self.tolerance,
            );
            outer.push(inside.then_some(surface));
        }
        // The outer region each inner region joins.
        let mut host: Vec<Option<usize>> = vec![None; regions];
        for r in (0..regions).filter(|&r| outer[r].is_none()) {
            let shell = self.built.region_shells[r][0];
            let q = &self.built.regions[r].interior_point;
            let candidates: Vec<usize> = (0..regions)
                .filter(|&o| {
                    shell_component(self.built.region_shells[o][0]) != shell_component(shell)
                        && outer[o]
                            .as_ref()
                            .is_some_and(|b| encloses(b, q, &self.tolerance))
                })
                .collect();
            // The innermost: inside every other candidate.
            let innermost = candidates.iter().copied().find(|&o| {
                candidates.iter().all(|&o2| {
                    o2 == o
                        || outer[o2].as_ref().is_some_and(|b| {
                            encloses(b, &self.built.regions[o].interior_point, &self.tolerance)
                        })
                })
            });
            let Some(o) = innermost else {
                return Err(declined(format!(
                    "shell {shell} bounds a region from inside, but no other shell encloses it"
                )));
            };
            if self.built.regions[o].bodies != self.built.regions[r].bodies {
                return Err(declined(format!(
                    "shell {shell} lies in a region of other bodies ({:?}) than its own ({:?})",
                    self.built.regions[o].bodies, self.built.regions[r].bodies
                )));
            }
            host[r] = Some(o);
        }
        // Fold each inner region's shell into its host.
        let mut shells = std::mem::take(&mut self.built.region_shells);
        for r in 0..regions {
            if let Some(o) = host[r] {
                let moved = std::mem::take(&mut shells[r]);
                shells[o].extend(moved);
            }
        }
        let mut infos = std::mem::take(&mut self.built.regions).into_iter();
        for (r, region_shells) in shells.into_iter().enumerate() {
            let info = infos.next();
            if host[r].is_some() {
                continue;
            }
            if let Some(info) = info {
                self.built.regions.push(info);
                self.built.region_shells.push(region_shells);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "region_nesting_tests.rs"]
mod tests;
