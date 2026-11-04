//! Bodies that touch over a face, analysed bonded (PK-981).
//!
//! An assembly's bodies pass no load or heat between them by default:
//! each is its own part and must be supported. Under `bodies =
//! "bonded"` the placed bodies are partitioned together (one cell
//! complex, `popkern_bool::partition`). The faces two bodies' regions
//! share are their bonded interfaces, and each group of bodies joined
//! by interfaces is analysed as the one solid their fusion makes. One
//! material holds for the whole setup, so the fused solid is the bonded
//! model: material runs on across the interface, as a conformal mesh
//! of the bodies sharing the interface's nodes would make it.
//!
//! Bodies touching only along an edge or at a point share no area,
//! carry nothing across, and stay separate. Bodies that overlap are an
//! interference, not a bond: they decline, naming both. A partition or
//! fusion the kernel declines declines the run with the kernel's reason
//! rather than falling back to separate bodies, whose answer would be a
//! different physical model from the one asked for.

use std::collections::BTreeMap;

use popkern_bool::{BoolOptions, NonManifoldPolicy, SetExpr, partition};
use popkern_core::topology::BRep;

use popkern_analysis::SolidId;
use popkern_analysis::contact::ContactPatch;

use crate::step_model::PlacedBody;

/// Two bodies bonded over the faces they share.
#[derive(Debug, Clone, PartialEq)]
pub struct BondedInterface {
    /// The two bodies' names, in input order.
    pub bodies: (String, String),
    /// The shared area, mm².
    pub area_mm2: f64,
}

/// The bodies to analyse after bonding, and the interfaces bonded.
#[derive(Debug, Clone)]
pub struct Bonded {
    /// One body per touching group: a group of one is its body as
    /// given, a larger group the fusion of its bodies, named by its
    /// bodies' names joined with ` + `.
    pub bodies: Vec<PlacedBody>,
    /// Every pair of bodies bonded, in input order.
    pub interfaces: Vec<BondedInterface>,
    /// The shared patches as contacts ([`in_contact`]; empty when
    /// bonded).
    pub patches: Vec<ContactPatch>,
}

/// Why bodies could not be bonded.
#[derive(Debug, thiserror::Error)]
pub enum BondError {
    /// Two bodies overlap over a volume.
    #[error("bodies `{a}` and `{b}` overlap: an interference, not a bonded contact")]
    Interference {
        /// The first body.
        a: String,
        /// The second body.
        b: String,
    },
    /// The kernel declined the partition, a fusion or a measure.
    #[error("bonding declined: {0}")]
    Kernel(String),
}

fn kernel(e: impl std::fmt::Display) -> BondError {
    BondError::Kernel(e.to_string())
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

/// The bodies of `bodies`, each solid its own, bonded where they share
/// face area (module doc).
///
/// # Errors
///
/// [`BondError::Interference`] when two bodies overlap;
/// [`BondError::Kernel`] when the kernel declines.
pub fn bond_touching(bodies: Vec<PlacedBody>) -> Result<Bonded, BondError> {
    touching(bodies, false)
}

/// The bodies of `bodies`, each solid its own and kept apart, with every
/// patch two of them share as a contact (`popkern_analysis::contact`):
/// tessellated, wound with its normal out of the first body.
///
/// # Errors
///
/// As [`bond_touching`].
pub fn in_contact(bodies: Vec<PlacedBody>) -> Result<Bonded, BondError> {
    touching(bodies, true)
}

fn touching(bodies: Vec<PlacedBody>, contacts: bool) -> Result<Bonded, BondError> {
    let mut solids: Vec<PlacedBody> = Vec::new();
    for body in bodies {
        if body.brep.solids.len() <= 1 {
            solids.push(body);
            continue;
        }
        let pieces = popkern_bool::pieces::piece_bodies(&body.brep).map_err(kernel)?;
        for (k, brep) in pieces.into_iter().enumerate() {
            solids.push(PlacedBody {
                name: format!("{}#{}", body.name, k + 1),
                brep,
            });
        }
    }
    // Names the report can tell apart: a name several solids share (the
    // solids of one product) is numbered in order.
    let mut seen: BTreeMap<String, usize> = BTreeMap::new();
    for s in &solids {
        *seen.entry(s.name.clone()).or_default() += 1;
    }
    let mut numbered: BTreeMap<String, usize> = BTreeMap::new();
    for s in &mut solids {
        if seen.get(&s.name).is_some_and(|&n| n > 1) {
            let k = numbered.entry(s.name.clone()).or_default();
            *k += 1;
            s.name = format!("{}#{k}", s.name);
        }
    }
    if solids.len() < 2 {
        return Ok(Bonded {
            bodies: solids,
            interfaces: Vec::new(),
            patches: Vec::new(),
        });
    }
    // The loosest body tolerance: what every body's geometry is good to.
    let tolerance = solids
        .iter()
        .map(|s| s.brep.tolerance)
        .max_by(|a, b| a.linear().total_cmp(&b.linear()))
        .unwrap_or_default();
    let refs: Vec<&BRep> = solids.iter().map(|s| &s.brep).collect();
    let complex = partition(&refs, tolerance).map_err(kernel)?;
    for region in complex.regions() {
        let mut held = region.bodies.iter();
        if let (Some(a), Some(b)) = (held.next(), held.next()) {
            return Err(BondError::Interference {
                a: solids[a].name.clone(),
                b: solids[b].name.clone(),
            });
        }
    }
    let regions = complex.region_solids(&SetExpr::Any).map_err(kernel)?;
    let body_of = |i: usize| regions.solids.get(i).and_then(|r| r.bodies.iter().next());
    let mut parent: Vec<usize> = (0..solids.len()).collect();
    let mut shared: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    let mut patches = Vec::new();
    for interface in &regions.interfaces {
        let (Some(a), Some(b)) = (body_of(interface.solids.0), body_of(interface.solids.1)) else {
            continue;
        };
        if a == b {
            continue;
        }
        let mut area = 0.0;
        for face in interface.patch.faces.keys() {
            area +=
                popkern_measure::face_area(&interface.patch, face, &tolerance).map_err(kernel)?;
        }
        *shared.entry((a.min(b), a.max(b))).or_default() += area;
        if contacts {
            patches.push(contact_patch(&interface.patch, (a, b))?);
        }
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        parent[ra.max(rb)] = ra.min(rb);
    }
    let interfaces = shared
        .into_iter()
        .map(|((a, b), area_mm2)| BondedInterface {
            bodies: (solids[a].name.clone(), solids[b].name.clone()),
            area_mm2,
        })
        .collect();
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..solids.len() {
        let root = find(&mut parent, i);
        groups.entry(root).or_default().push(i);
    }
    if contacts {
        return Ok(Bonded {
            bodies: solids,
            interfaces,
            patches,
        });
    }
    let options = BoolOptions {
        non_manifold: NonManifoldPolicy::SplitToManifold,
    };
    let mut out = Vec::with_capacity(groups.len());
    for members in groups.values() {
        if let [only] = members.as_slice() {
            out.push(solids[*only].clone());
            continue;
        }
        let brep = complex
            .select_where(|held| held.iter().any(|i| members.contains(&i)), &options)
            .map_err(kernel)?;
        let name = members
            .iter()
            .map(|&i| solids[i].name.as_str())
            .collect::<Vec<_>>()
            .join(" + ");
        out.push(PlacedBody { name, brep });
    }
    Ok(Bonded {
        bodies: out,
        interfaces,
        patches: Vec::new(),
    })
}

/// `patch` (wound in the first body's sense) as the solver's contact
/// between bodies `a` and `b`: its triangles, each wound with the mesh
/// normals there (the face's outward sense, out of `a`).
fn contact_patch(patch: &BRep, (a, b): (usize, usize)) -> Result<ContactPatch, BondError> {
    let mesh = popkern_mesh::tessellate(patch, &popkern_mesh::TessellationOptions::display())
        .map_err(kernel)?;
    let at = |i: usize| mesh.vertices[i].position;
    let triangles = mesh
        .triangles
        .iter()
        .map(|t| {
            let (p, q, r) = (at(t[0]), at(t[1]), at(t[2]));
            let n = (q - p).cross(&(r - p));
            let outward = t.iter().fold(
                popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0),
                |acc, &i| acc + mesh.vertices[i].normal,
            );
            let corners = if n.dot(&outward) < 0.0 {
                [p, r, q]
            } else {
                [p, q, r]
            };
            corners.map(|c| [c.x, c.y, c.z])
        })
        .collect();
    Ok(ContactPatch {
        solids: (SolidId(a as u32), SolidId(b as u32)),
        triangles,
    })
}
