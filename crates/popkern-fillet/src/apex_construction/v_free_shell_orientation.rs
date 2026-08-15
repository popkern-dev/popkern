//! Shell winding reconciliation for the V-free consolidation (PK-153).
//!
//! The V-free repeated-apex rewrite authors its replacement shell
//! (strips, caps, face-family bridges) with loop winding mirrored
//! against the untouched body faces: the face-family bridge is built
//! as the twin cycle of the source support loop, and the post-bridge
//! mixed-winding replan then locks every downstream loop (cylinder
//! trims, cap boundary, support assemblies) to that mirror class.
//! Live twins stay antiparallel inside each class, so the model
//! validates — but every seam between the two classes is a coincident
//! PARALLEL dead pair that can never weld, and a face forced to
//! interface both classes contorts (slit loops, dead diagonal chords).
//!
//! Repair, in class terms rather than fixture terms:
//! 1. 2-color faces over seam parity constraints: a live twin pair is
//!    a same-class relation; a coincident parallel dead pair is
//!    opposite-class; a coincident antiparallel dead pair same-class.
//! 2. Faces that make the coloring infeasible AND carry structural
//!    contortion evidence (an edge traversed twice in one loop, or a
//!    partnerless dead side) are the broken interfaces: their
//!    constraints are excluded and their rings queued for rebuild.
//! 3. The mirror class is identified by the plane-winding invariant:
//!    on every default-construction body the outer-loop Newell normal
//!    of a plane face agrees in sign with the stored plane normal.
//!    The class whose plane faces unanimously violate that invariant
//!    reverses; ambiguity declines the stage (traced, no mutation).
//! 4. Reversal is twin-exchange per loop: interior seams swap halves
//!    between the two adjacent loops (staying antiparallel), dead
//!    sides flip in place, so the class boundary lands weldable.
//! 5. Queued interface faces get rings rebuilt by boundary
//!    completion: trusted live sides are kept and the ring closes by
//!    authoring twin halves against coincident dead sides on the
//!    face's own surface.
//!
//! The owning consolidation pass welds the reconciled seams, merges
//! seam vertices, heals pcurves, and validates-or-rolls-back.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::{face_sides, side_span, twin_face};

#[path = "v_free_shell_reversal.rs"]
mod reversal;
#[path = "v_free_shell_ring_rebuild.rs"]
mod ring_rebuild;
use reversal::reverse_class_loops;
use ring_rebuild::rebuild_face_ring;
pub(super) use ring_rebuild::rebuild_partnerless_chord_rings;

/// Structured counts for the winding reconciliation stage.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(super) struct WindingReconciliationOutcome {
    /// Faces whose loops were reversed (the mirror class).
    pub(super) reversed_faces: usize,
    /// Live cross-class seams severed ahead of ring rebuilds.
    pub(super) severed_seams: usize,
    /// Contorted interface faces whose rings were rebuilt.
    pub(super) rebuilt_rings: usize,
    /// Vertex handles that welds/rebuilds may have duplicated; the
    /// caller merges coincident records among these.
    pub(super) seam_vertices: Vec<VertexHandle>,
    /// Every face this reconciliation mutated (reversed class +
    /// contorted rebuild targets) — the caller's dirty set for
    /// region-scoped pcurve healing (PK-274).
    pub(super) touched_faces: Vec<FaceHandle>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SeamParity {
    Same,
    Opposite,
}

pub(super) fn reconcile_shell_winding(brep: &mut BRep) -> WindingReconciliationOutcome {
    let mut outcome = WindingReconciliationOutcome::default();
    let constraints = seam_parity_constraints(brep);
    if !constraints.iter().any(|c| c.parity == SeamParity::Opposite) {
        return outcome;
    }
    let contorted: BTreeSet<FaceHandle> = brep
        .faces
        .keys()
        .filter(|f| face_has_contorted_ring(brep, *f))
        .collect();
    let usable: Vec<&SeamConstraint> = constraints
        .iter()
        .filter(|c| !contorted.contains(&c.a) && !contorted.contains(&c.b))
        .collect();
    let mut forest = ParityForest::default();
    for c in &usable {
        if forest.union(c.a, c.b, c.parity).is_err() {
            trace(&format!(
                "winding reconciliation declined: parity conflict on {:?}~{:?} outside \
                 contorted set {:?}",
                c.a, c.b, contorted
            ));
            return outcome;
        }
    }
    let Some(reversal_faces) = select_reversal_class(brep, &mut forest, &usable) else {
        return outcome;
    };
    if reversal_faces.is_empty() {
        return outcome;
    }

    // Capture trusted sides of contorted faces, then sever their
    // remaining sides so the reversal sees clean class boundaries.
    let mut rebuild_plans: Vec<(FaceHandle, Vec<HalfEdgeHandle>)> = Vec::new();
    for &face in &contorted {
        let trusted = trusted_sides(brep, face, &reversal_faces, &contorted);
        if trusted.is_empty() {
            trace(&format!(
                "winding reconciliation declined: contorted face {face:?} has no trusted sides"
            ));
            return outcome;
        }
        rebuild_plans.push((face, trusted));
    }
    for (face, trusted) in &rebuild_plans {
        let keep: BTreeSet<HalfEdgeHandle> = trusted.iter().copied().collect();
        for side in face_sides(brep, *face) {
            if !keep.contains(&side) {
                remove_side(brep, side);
                outcome.severed_seams += 1;
            }
        }
    }
    // A live seam from the reversal class to an untouched face outside
    // both the class and the contorted set cannot be reversed safely.
    for &face in &reversal_faces {
        for side in face_sides(brep, face) {
            if let Some(partner) = twin_face(brep, side)
                && partner != face
                && !reversal_faces.contains(&partner)
            {
                trace(&format!(
                    "winding reconciliation declined: live seam from reversal face {face:?} \
                     to out-of-class face {partner:?}"
                ));
                return outcome;
            }
        }
    }

    outcome.reversed_faces = reverse_class_loops(brep, &reversal_faces);
    outcome.touched_faces = reversal_faces
        .iter()
        .copied()
        .chain(rebuild_plans.iter().map(|(face, _)| *face))
        .collect();
    for (face, trusted) in rebuild_plans {
        match rebuild_face_ring(brep, face, &trusted) {
            Some(mut seam) => {
                outcome.rebuilt_rings += 1;
                outcome.seam_vertices.append(&mut seam);
            }
            None => {
                trace(&format!(
                    "winding reconciliation: ring rebuild declined for {face:?}"
                ));
            }
        }
    }
    trace(&format!(
        "winding reconciliation: reversed {} face(s), severed {}, rebuilt {} ring(s)",
        outcome.reversed_faces, outcome.severed_seams, outcome.rebuilt_rings
    ));
    outcome
}

struct SeamConstraint {
    a: FaceHandle,
    b: FaceHandle,
    parity: SeamParity,
}

fn seam_parity_constraints(brep: &BRep) -> Vec<SeamConstraint> {
    let mut constraints = Vec::new();
    let mut seen_edges = BTreeSet::new();
    for (handle, data) in &brep.half_edges {
        if let Some(twin) = brep.half_edges.get(data.twin) {
            if !seen_edges.insert(data.edge) {
                continue;
            }
            let (Some(a), Some(b)) = (
                brep.loops.get(data.loop_).map(|l| l.face),
                brep.loops.get(twin.loop_).map(|l| l.face),
            ) else {
                continue;
            };
            if a != b {
                constraints.push(SeamConstraint {
                    a,
                    b,
                    parity: SeamParity::Same,
                });
            }
        } else {
            let _ = handle;
        }
    }
    let dead: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, data)| !brep.half_edges.contains_key(data.twin))
        .map(|(handle, _)| handle)
        .collect();
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    for i in 0..dead.len() {
        let Some((a_start, a_end)) = side_span(brep, dead[i]) else {
            continue;
        };
        for j in (i + 1)..dead.len() {
            let Some((b_start, b_end)) = side_span(brep, dead[j]) else {
                continue;
            };
            let antiparallel =
                a_start.distance_to(&b_end) <= tol && a_end.distance_to(&b_start) <= tol;
            let parallel = a_start.distance_to(&b_start) <= tol && a_end.distance_to(&b_end) <= tol;
            if !antiparallel && !parallel {
                continue;
            }
            let (Some(a), Some(b)) = (side_face(brep, dead[i]), side_face(brep, dead[j])) else {
                continue;
            };
            if a == b {
                continue;
            }
            constraints.push(SeamConstraint {
                a,
                b,
                parity: if parallel {
                    SeamParity::Opposite
                } else {
                    SeamParity::Same
                },
            });
        }
    }
    constraints
}

fn side_face(brep: &BRep, side: HalfEdgeHandle) -> Option<FaceHandle> {
    brep.loops
        .get(brep.half_edges.get(side)?.loop_)
        .map(|l| l.face)
}

/// Structural contortion evidence: an edge traversed twice within the
/// face's own loops (a slit — the shape a face takes when forced to
/// interface both winding classes). Partnerless dead sides are NOT
/// evidence — the victims of a contorted neighbor (e.g. strip rails
/// missing their support twin) carry those too.
fn face_has_contorted_ring(brep: &BRep, face: FaceHandle) -> bool {
    let sides = face_sides(brep, face);
    if sides.is_empty() {
        return false;
    }
    let mut edge_counts: BTreeMap<_, usize> = BTreeMap::new();
    for side in &sides {
        if let Some(data) = brep.half_edges.get(*side) {
            *edge_counts.entry(data.edge).or_default() += 1;
        }
    }
    edge_counts.values().any(|&count| count >= 2)
}

#[derive(Default)]
struct ParityForest {
    parent: BTreeMap<FaceHandle, (FaceHandle, bool)>,
}

impl ParityForest {
    fn find(&mut self, face: FaceHandle) -> (FaceHandle, bool) {
        let (parent, parity) = *self.parent.entry(face).or_insert((face, false));
        if parent == face {
            return (face, false);
        }
        let (root, root_parity) = self.find(parent);
        let total = parity ^ root_parity;
        self.parent.insert(face, (root, total));
        (root, total)
    }

    fn union(&mut self, a: FaceHandle, b: FaceHandle, parity: SeamParity) -> Result<(), ()> {
        let flip = parity == SeamParity::Opposite;
        let (ra, pa) = self.find(a);
        let (rb, pb) = self.find(b);
        if ra == rb {
            if (pa ^ pb) != flip {
                return Err(());
            }
            return Ok(());
        }
        self.parent.insert(ra, (rb, pa ^ pb ^ flip));
        Ok(())
    }
}

/// Pick the parity class to reverse via the plane-winding invariant.
/// Requires the two classes' plane faces to vote unanimously and
/// oppositely; anything else declines with a trace.
fn select_reversal_class(
    brep: &BRep,
    forest: &mut ParityForest,
    constraints: &[&SeamConstraint],
) -> Option<BTreeSet<FaceHandle>> {
    let mut faces = BTreeSet::new();
    for c in constraints {
        faces.insert(c.a);
        faces.insert(c.b);
    }
    // Group by root; only components that contain an Opposite relation
    // need reversal work.
    let mut components: BTreeMap<FaceHandle, Vec<(FaceHandle, bool)>> = BTreeMap::new();
    for face in faces {
        let (root, parity) = forest.find(face);
        components.entry(root).or_default().push((face, parity));
    }
    let mut split_roots = BTreeSet::new();
    for c in constraints {
        if c.parity == SeamParity::Opposite {
            split_roots.insert(forest.find(c.a).0);
        }
    }
    let mut reversal = BTreeSet::new();
    for root in split_roots {
        let members = components.remove(&root)?;
        // PK-260: calibrate with magnitude-gated votes. The parity
        // classes are exact topology, but a near-degenerate sliver's
        // Newell normal is noise — on the oblique base-corner class
        // three slivers voted "honest" at |dot| 0.0015..0.047 against
        // three genuinely mirror-wound production faces at 0.69..1.93
        // in the SAME parity class, splitting the calibration and
        // starving the reversal. Votes below 10% of the component's
        // strongest magnitude abstain.
        let ballots: Vec<(bool, bool, f64)> = members
            .iter()
            .filter_map(|(face, parity)| {
                let honest = plane_winding_honest(brep, *face)?;
                let magnitude = plane_winding_newell_magnitude(brep, *face)?.abs();
                Some((*parity, honest, magnitude))
            })
            .collect();
        let strongest = ballots
            .iter()
            .map(|(_, _, magnitude)| *magnitude)
            .fold(0.0_f64, f64::max);
        let mut votes: BTreeMap<bool, (usize, usize)> = BTreeMap::new();
        for (parity, honest, magnitude) in &ballots {
            if *magnitude < strongest * 0.1 {
                continue;
            }
            let entry = votes.entry(*parity).or_default();
            if *honest {
                entry.0 += 1;
            } else {
                entry.1 += 1;
            }
        }
        let dishonest_classes: Vec<bool> = votes
            .iter()
            .filter(|(_, (honest, dishonest))| *dishonest > 0 && *honest == 0)
            .map(|(class, _)| *class)
            .collect();
        let honest_classes: Vec<bool> = votes
            .iter()
            .filter(|(_, (honest, dishonest))| *honest > 0 && *dishonest == 0)
            .map(|(class, _)| *class)
            .collect();
        let (&[mirror], &[_keep]) = (dishonest_classes.as_slice(), honest_classes.as_slice())
        else {
            trace(&format!(
                "winding reconciliation declined: ambiguous plane calibration votes {votes:?}"
            ));
            // PK-260 calibration census: per-member verdicts with the
            // Newell magnitude — near-degenerate sliver rings vote
            // with noise-level normals and can split a class.
            if crate::direct_apex_trace::enabled() {
                for (face, parity) in &members {
                    let verdict = plane_winding_honest(brep, *face);
                    let magnitude = plane_winding_newell_magnitude(brep, *face);
                    trace(&format!(
                        "  calibration vote {face:?} parity={parity} honest={verdict:?} \
                         newell_dot={magnitude:?}"
                    ));
                }
            }
            return None;
        };
        for (face, parity) in members {
            if parity == mirror {
                reversal.insert(face);
            }
        }
    }
    Some(reversal)
}

/// The raw Newell·normal magnitude behind [`plane_winding_honest`] —
/// a signed-area scale that separates reliable votes from
/// near-degenerate sliver noise.
fn plane_winding_newell_magnitude(brep: &BRep, face: FaceHandle) -> Option<f64> {
    let sof = brep.faces.get(face)?.surface.as_ref()?;
    let SurfaceGeometry::Plane { normal, .. } = &sof.geometry else {
        return None;
    };
    let ring = brep.face_outer_loop_positions(face);
    if ring.len() < 3 {
        return None;
    }
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..ring.len() {
        let p = ring[i];
        let q = ring[(i + 1) % ring.len()];
        newell.x += (p.y - q.y) * (p.z + q.z);
        newell.y += (p.z - q.z) * (p.x + q.x);
        newell.z += (p.x - q.x) * (p.y + q.y);
    }
    Some(newell.dot(normal))
}

/// Plane-winding honesty: Newell normal of the outer ring agrees in
/// sign with the stored plane normal exactly when the orientation
/// flag says the surface normal agrees with the face normal.
fn plane_winding_honest(brep: &BRep, face: FaceHandle) -> Option<bool> {
    let sof = brep.faces.get(face)?.surface.as_ref()?;
    let SurfaceGeometry::Plane { normal, .. } = &sof.geometry else {
        return None;
    };
    let ring = brep.face_outer_loop_positions(face);
    if ring.len() < 3 {
        return None;
    }
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..ring.len() {
        let p = ring[i];
        let q = ring[(i + 1) % ring.len()];
        newell.x += (p.y - q.y) * (p.z + q.z);
        newell.y += (p.z - q.z) * (p.x + q.x);
        newell.z += (p.x - q.x) * (p.y + q.y);
    }
    let dot = newell.dot(normal);
    if dot.abs() < 1.0e-12 {
        return None;
    }
    Some((dot > 0.0) == sof.orientation)
}

fn trusted_sides(
    brep: &BRep,
    face: FaceHandle,
    reversal: &BTreeSet<FaceHandle>,
    contorted: &BTreeSet<FaceHandle>,
) -> Vec<HalfEdgeHandle> {
    face_sides(brep, face)
        .into_iter()
        .filter(|side| {
            twin_face(brep, *side).is_some_and(|partner| {
                partner != face && !reversal.contains(&partner) && !contorted.contains(&partner)
            })
        })
        .collect()
}

/// Remove one half-edge; a surviving twin becomes a dead boundary
/// side, a fully-unreferenced edge record is dropped.
fn remove_side(brep: &mut BRep, side: HalfEdgeHandle) {
    let Some(data) = brep.half_edges.remove(side) else {
        return;
    };
    if brep.half_edges.contains_key(data.twin) {
        brep.boundary_edges.insert(data.edge);
        if let Some(edge_data) = brep.edges.get_mut(data.edge)
            && edge_data.half_edge == side
        {
            edge_data.half_edge = data.twin;
        }
    } else if brep.edges.remove(data.edge).is_some() {
        brep.boundary_edges.remove(&data.edge);
    }
}

fn trace(message: &str) {
    if crate::direct_apex_trace::enabled() {
        eprintln!("[direct-apex-trace] {message}");
    }
}
