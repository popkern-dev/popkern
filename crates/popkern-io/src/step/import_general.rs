//! Importing a self-touching shell as one general body (PK-1108).
//!
//! A `CLOSED_SHELL` listing the faces of solids that touch (two blocks
//! sharing an edge or a corner, as OCC writes them) maps as manifold
//! pieces meeting at recorded contacts (`shell_pieces`). Under
//! [`StepNonManifold::General`] the pieces are glued at those contacts
//! (`BRep::glue_contacts`): one general body whose contact edge's ring
//! holds every piece's faces in radial order. A glue that declines, or a
//! glued body that does not validate, keeps the pieces and says why.

use popkern_core::topology::BRep;

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// How a STEP shell that touches itself imports.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StepNonManifold {
    /// Manifold pieces meeting at recorded contacts.
    #[default]
    Pieces,
    /// One general body: the pieces glued at their contacts.
    General,
}

/// Glue `brep` into one general body where its pieces touch: its
/// recorded contacts (a shell split into pieces), then whatever
/// coincides without being shared (`BRep::glue_coincident`: a writer's
/// cavities touching along an edge each with their own copy of it), then
/// the cells a face-by-face shell encloses (`BRep::fill_enclosed_cells`).
/// Kept when it validates; logged against STEP entity `id` either way.
pub(crate) fn glue_pieces(brep: &mut BRep, id: u64, log: &mut DiagnosticLog) {
    let mut glued = brep.clone();
    let mut steps = Vec::new();
    let outcome = (|| -> Result<(), String> {
        if !glued.topology_contacts.is_empty() {
            let report = glued.glue_contacts().map_err(|e| e.to_string())?;
            steps.push(format!(
                "{} contact(s) glued, {} kept recorded",
                report.glued,
                report.kept.len()
            ));
        }
        let coincident = glued.glue_coincident().map_err(|e| e.to_string())?;
        if coincident.vertices + coincident.edges > 0 {
            steps.push(format!(
                "{} coincident vertex(es) and {} edge(s) glued",
                coincident.vertices, coincident.edges
            ));
        }
        if steps.is_empty() {
            return Ok(());
        }
        // A shell imported face by face glues into sheets: the cells
        // they enclose become its solids.
        glued.fill_enclosed_cells();
        let issues = glued.validate_all();
        if issues.is_empty() {
            Ok(())
        } else {
            Err(format!("glued body invalid: {}", issues[0]))
        }
    })();
    match outcome {
        Ok(()) if steps.is_empty() => {}
        Ok(()) => {
            *brep = glued;
            log.warn(
                Phase::MapTopology,
                Some(id),
                format!(
                    "touching pieces imported as one general body: {}",
                    steps.join("; ")
                ),
            );
        }
        Err(why) => log.warn(
            Phase::MapTopology,
            Some(id),
            format!("touching pieces kept as they are: {why}"),
        ),
    }
}

/// Under [`StepNonManifold::General`], the items of one representation
/// that make one general body (a solid, the sheets, wires and acorns
/// `export_general` writes beside it) are put together again: items
/// touching where their vertices coincide are grouped, a group holding
/// no solid joins the group whose solid region contains it (a lamina,
/// wire or acorn floating in a solid), and each group of several is
/// appended into one body, glued where it coincides
/// (`BRep::glue_coincident`, a sheet glued on a solid's material side
/// becoming its lamina) and its cells filled (`BRep::fill_enclosed_cells`,
/// free shells joining the cells holding them). A group whose body does
/// not validate keeps its items, with a diagnostic. `None` when nothing
/// changed.
pub(crate) fn merge_touching_items(items: &[BRep], log: &mut DiagnosticLog) -> Option<Vec<BRep>> {
    let n = items.len();
    if n < 2 {
        return None;
    }
    let mut group: Vec<usize> = (0..n).collect();
    let find = |group: &mut Vec<usize>, mut i: usize| {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    };
    let union = |group: &mut Vec<usize>, a: usize, b: usize| {
        let (ra, rb) = (find(group, a), find(group, b));
        if ra != rb {
            group[ra.max(rb)] = ra.min(rb);
        }
    };
    // Coincident vertices, bucketed on a grid of the largest tolerance.
    let tolerance = items
        .iter()
        .map(|b| b.tolerance.linear())
        .fold(1e-9, f64::max);
    let cell = |p: popkern_core::geometry::Point3| {
        let k = |x: f64| (x / tolerance).floor() as i64;
        (k(p.x), k(p.y), k(p.z))
    };
    let mut grid: std::collections::HashMap<
        (i64, i64, i64),
        Vec<(usize, popkern_core::geometry::Point3)>,
    > = std::collections::HashMap::new();
    for (i, body) in items.iter().enumerate() {
        for v in body.vertices.values() {
            grid.entry(cell(v.position))
                .or_default()
                .push((i, v.position));
        }
    }
    let mut touching = Vec::new();
    for (&(x, y, z), here) in &grid {
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(near) = grid.get(&(x + dx, y + dy, z + dz)) else {
                        continue;
                    };
                    for &(i, p) in here {
                        for &(j, q) in near {
                            if i < j && p.distance_to(&q) <= tolerance {
                                touching.push((i, j));
                            }
                        }
                    }
                }
            }
        }
    }
    for (i, j) in touching {
        union(&mut group, i, j);
    }
    // Groups without a solid inside another group's solid region.
    let has_solid = |members: &[usize]| {
        members
            .iter()
            .any(|&i| !items[i].solid_regions().is_empty())
    };
    let mut members: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        members.entry(find(&mut group, i)).or_default().push(i);
    }
    let free: Vec<usize> = members
        .iter()
        .filter(|(_, m)| !has_solid(m))
        .map(|(&root, _)| root)
        .collect();
    for root in free {
        let Some(q) = items[root].vertices.values().next().map(|v| v.position) else {
            continue;
        };
        let host = (0..n).find(|&j| {
            find(&mut group, j) != find(&mut group, root)
                && items[j]
                    .solid_regions()
                    .into_iter()
                    .any(|r| items[j].region_winding(r, q) > 0.5)
        });
        if let Some(j) = host {
            union(&mut group, root, j);
        }
    }
    let mut members: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        members.entry(find(&mut group, i)).or_default().push(i);
    }
    if members.values().all(|m| m.len() == 1) {
        return None;
    }
    let mut out = Vec::new();
    for m in members.values() {
        if m.len() == 1 {
            out.push(items[m[0]].clone());
            continue;
        }
        let mut body = items[m[0]].clone();
        for &i in &m[1..] {
            body.append_body(&items[i]);
        }
        let glued = body.glue_coincident().map_err(|e| e.to_string());
        let outcome = glued.and_then(|glued| {
            body.fill_enclosed_cells();
            match body.validate_all().first() {
                None => Ok(glued),
                Some(issue) => Err(format!("merged body invalid: {issue}")),
            }
        });
        match outcome {
            Ok(glued) => {
                log.warn(
                    Phase::MapProduct,
                    None,
                    format!(
                        "{} representation items imported as one general body: {} coincident \
                         vertex(es) and {} edge(s) glued",
                        m.len(),
                        glued.vertices,
                        glued.edges
                    ),
                );
                out.push(body);
            }
            Err(why) => {
                log.warn(
                    Phase::MapProduct,
                    None,
                    format!(
                        "{} touching representation items kept apart: {why}",
                        m.len()
                    ),
                );
                out.extend(m.iter().map(|&i| items[i].clone()));
            }
        }
    }
    Some(out)
}
