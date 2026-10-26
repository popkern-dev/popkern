//! The design's pieces: which closed surfaces of the design body carry
//! load.
//!
//! A density design can come apart: material the solver left floating
//! (an island touching no support carries nothing), voids it closed in
//! (a cavity no printer can drain), or, under too few updates, a load
//! the material never joined to a support. The design surface's
//! edge-connected shells are sorted by the volume they enclose: an
//! inward shell is a cavity and is filled (its triangles dropped, the
//! enclosing shell solid there); an outward shell is kept when it holds
//! a triangle on a supported face and dropped otherwise. The design
//! declines when no kept shell reaches a support, or a loaded face's
//! material is only in dropped ones. Dropped and filled volumes are
//! reported, never silent.

use std::collections::HashMap;

/// What the sorting kept and why it dropped the rest.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Pieces {
    /// The kept triangles (indices into the design surface).
    pub triangles: Vec<usize>,
    /// Kept shells.
    pub kept: usize,
    /// Volume of the dropped floating shells, mm³.
    pub detached_volume: f64,
    /// Volume of the filled cavities, mm³.
    pub cavity_volume: f64,
}

/// Why the design does not connect its loads to its supports.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Disconnected {
    /// No material reaches a supported face.
    #[error("no material of the design reaches a support")]
    NoSupport,
    /// A loaded face's material is not joined to any support.
    #[error("the material at load {load} is not joined to a support")]
    LoadDetached {
        /// The load's index in the request.
        load: usize,
    },
}

fn signed_volume(vertices: &[[f64; 3]], triangles: &[[u32; 3]], shell: &[usize]) -> f64 {
    shell
        .iter()
        .map(|&t| {
            let [a, b, c] = triangles[t].map(|i| vertices[i as usize]);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

/// The edge-connected shells of a closed mesh, each a list of
/// triangles, in order of their first triangle.
fn shells(triangles: &[[u32; 3]]) -> Vec<Vec<usize>> {
    let mut by_edge: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    for (t, tri) in triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            by_edge.entry((a.min(b), a.max(b))).or_default().push(t);
        }
    }
    let mut shell_of = vec![usize::MAX; triangles.len()];
    let mut out: Vec<Vec<usize>> = Vec::new();
    for start in 0..triangles.len() {
        if shell_of[start] != usize::MAX {
            continue;
        }
        let s = out.len();
        shell_of[start] = s;
        let mut members = vec![start];
        let mut i = 0;
        while i < members.len() {
            let tri = triangles[members[i]];
            i += 1;
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                for &n in by_edge
                    .get(&(a.min(b), a.max(b)))
                    .map_or(&[][..], Vec::as_slice)
                {
                    if shell_of[n] == usize::MAX {
                        shell_of[n] = s;
                        members.push(n);
                    }
                }
            }
        }
        members.sort_unstable();
        out.push(members);
    }
    out
}

/// Sorts the design surface's shells: `supported(owner)` says a
/// triangle's owner is a supported face, `loads` lists each load's
/// owners.
pub(crate) fn sort(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    owners: &[u32],
    supported: &dyn Fn(u32) -> bool,
    loads: &[Vec<u32>],
) -> Result<Pieces, Disconnected> {
    let mut pieces = Pieces {
        triangles: Vec::new(),
        kept: 0,
        detached_volume: 0.0,
        cavity_volume: 0.0,
    };
    // Per load: whether any of its material is in a dropped shell, and
    // whether any is kept.
    let mut load_kept = vec![false; loads.len()];
    let mut load_dropped = vec![false; loads.len()];
    for shell in shells(triangles) {
        let volume = signed_volume(vertices, triangles, &shell);
        if volume < 0.0 {
            pieces.cavity_volume -= volume;
            continue;
        }
        let keep = shell.iter().any(|&t| supported(owners[t]));
        for (l, owners_of_load) in loads.iter().enumerate() {
            if shell.iter().any(|&t| owners_of_load.contains(&owners[t])) {
                if keep {
                    load_kept[l] = true;
                } else {
                    load_dropped[l] = true;
                }
            }
        }
        if keep {
            pieces.kept += 1;
            pieces.triangles.extend(shell);
        } else {
            pieces.detached_volume += volume;
        }
    }
    if pieces.kept == 0 {
        return Err(Disconnected::NoSupport);
    }
    if let Some(load) = (0..loads.len()).find(|&l| load_dropped[l] && !load_kept[l]) {
        return Err(Disconnected::LoadDetached { load });
    }
    pieces.triangles.sort_unstable();
    Ok(pieces)
}
