//! A design's history (PK-959): the density field at chosen updates of
//! its solve, kept so a host can show the design converging (a film).
//!
//! The intermediate designs are not bodies: early on most of the field
//! is grey, and its surface is no part anyone would make, so nothing
//! here facets or names them. They are display meshes: the iso-surface
//! at ½ of a recorded field ([`grid_surface`], the solver's own design
//! surface), or of two neighbouring records blended, so a film moves
//! smoothly between them. The design space's analysis mesh comes with
//! them, for a host to draw the space the design is carved from.
//!
//! The field changes fastest at the start, so the updates recorded are
//! every one of the first [`EVERY_UPDATE`], then each about 15% past
//! the last ([`next_recorded`]), and always the final design: 34 for a
//! 300-update solve. Each record is a byte per cell the part reaches
//! (the material quantised to 1/255): 15k cells for the gallery
//! bracket, about half a megabyte in all.

use std::ops::ControlFlow;

use popkern_analysis::model::AnalysisModel;
use popkern_analysis::structural::generative::surface::{MaterialGrid, grid_surface};
use popkern_analysis::structural::generative::{Design, Progress};

/// Updates recorded one by one from the first.
const EVERY_UPDATE: usize = 12;

/// The update recorded after `update`: the next one while under
/// [`EVERY_UPDATE`], else about 15% further on.
pub(crate) fn next_recorded(update: usize) -> usize {
    if update < EVERY_UPDATE {
        update + 1
    } else {
        update + (update * 3).div_ceil(20)
    }
}

/// One recorded design.
#[derive(Debug, Clone, PartialEq)]
pub struct HistoryFrame {
    /// Updates run when it was recorded (1 for the first).
    pub update: usize,
    /// Compliance of the design before that update, N·mm.
    pub compliance_nmm: f64,
    /// Its material volume over the design space's.
    pub volume_fraction: f64,
    /// Per cell of [`DesignHistory::active`]: its material, 0..=255.
    pub material: Vec<u8>,
}

/// A triangle mesh for display, in the design body's units.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct DisplayMesh {
    /// Vertex positions, xyz.
    pub positions: Vec<f32>,
    /// Unit vertex normals (area-weighted), xyz, pointing out of the
    /// material.
    pub normals: Vec<f32>,
    /// Triangles as vertex indices, counter-clockwise seen from outside.
    pub indices: Vec<u32>,
}

/// The fields a design's solve went through, and its design space.
#[derive(Debug, Clone, PartialEq)]
pub struct DesignHistory {
    /// Cell size, in the design body's units.
    pub cell: f64,
    /// The lower corner of cell (0, 0, 0), in the body's units.
    pub corner: [f64; 3],
    /// Cells along each axis (a zero layer around the design space).
    pub dims: [usize; 3],
    /// The cells the design space reaches (any fill), ascending: every
    /// other cell holds no material in any record.
    pub active: Vec<u32>,
    /// Updates the solve ran.
    pub updates: usize,
    /// The recorded designs, by update, the final one last.
    pub frames: Vec<HistoryFrame>,
    /// The design space's analysis mesh, in the body's units.
    pub space: DisplayMesh,
}

/// A field seen during a solve, before it is binned to cells.
struct Seen {
    update: usize,
    compliance: f64,
    volume_fraction: f64,
    densities: Vec<f32>,
}

/// Records a solve's fields at the updates [`next_recorded`] picks.
#[derive(Default)]
pub(crate) struct Recorder {
    seen: Vec<Seen>,
    next: usize,
}

impl Recorder {
    /// Wraps `observer`: every update goes on to it; the chosen ones
    /// are recorded first.
    pub fn observing<'a>(
        &'a mut self,
        observer: &'a mut dyn FnMut(Progress<'_>) -> ControlFlow<()>,
    ) -> impl FnMut(Progress<'_>) -> ControlFlow<()> + 'a {
        move |p: Progress<'_>| {
            self.see(&p);
            observer(p)
        }
    }

    fn see(&mut self, p: &Progress<'_>) {
        // A solve starting over (the shifted-grid fallback) starts the
        // record over.
        if p.count <= 1 {
            self.seen.clear();
            self.next = 1;
        }
        if p.count >= self.next {
            self.seen.push(Seen {
                update: p.count,
                compliance: p.iteration.compliance,
                volume_fraction: p.iteration.volume_fraction,
                densities: p.densities.iter().map(|&d| d as f32).collect(),
            });
            self.next = next_recorded(p.count);
        }
    }

    /// The history of `design`, the solve just recorded, on `model`
    /// (millimetres) in the body's units (`mm_per_unit`); `None` for a
    /// design without cells.
    pub fn finish(
        self,
        design: &Design,
        model: &AnalysisModel,
        mm_per_unit: f64,
    ) -> Option<DesignHistory> {
        let h = design.cell_size;
        let ones = vec![1.0; design.entries.len()];
        let fill = MaterialGrid::new(&design.entries, &ones, h)?;
        let active: Vec<u32> = (0..fill.material.len() as u32)
            .filter(|&c| fill.material[c as usize] > 0.0)
            .collect();
        let quantise = |densities: &[f64]| -> Option<Vec<u8>> {
            let grid = MaterialGrid::new(&design.entries, densities, h)?;
            Some(
                active
                    .iter()
                    .map(|&c| (grid.material[c as usize].clamp(0.0, 1.0) * 255.0).round() as u8)
                    .collect(),
            )
        };
        let mut frames = Vec::with_capacity(self.seen.len() + 1);
        for seen in &self.seen {
            let densities: Vec<f64> = seen.densities.iter().map(|&d| f64::from(d)).collect();
            frames.push(HistoryFrame {
                update: seen.update,
                compliance_nmm: seen.compliance,
                volume_fraction: seen.volume_fraction,
                material: quantise(&densities)?,
            });
        }
        // The final design, exactly (the last update seen may be
        // between records).
        let updates = design.history.len();
        frames.retain(|f| f.update < updates);
        frames.push(HistoryFrame {
            update: updates,
            compliance_nmm: design.compliance,
            volume_fraction: design.volume_fraction,
            material: quantise(&design.densities)?,
        });
        let scale = 1.0 / mm_per_unit;
        Some(DesignHistory {
            cell: h * scale,
            corner: fill.corner.map(|c| c * scale),
            dims: fill.dims,
            active,
            updates,
            frames,
            space: space_mesh(model, scale),
        })
    }
}

/// `model`'s triangles scaled by `scale`, every solid in one mesh.
fn space_mesh(model: &AnalysisModel, scale: f64) -> DisplayMesh {
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    for solid in &model.solids {
        let base = vertices.len() as u32;
        vertices.extend(solid.vertices.iter().map(|v| v.map(|x| x * scale)));
        triangles.extend(solid.triangles.iter().map(|t| t.map(|i| i + base)));
    }
    display_mesh(&vertices, &triangles)
}

/// `vertices` and `triangles` as display buffers with area-weighted
/// vertex normals.
fn display_mesh(vertices: &[[f64; 3]], triangles: &[[u32; 3]]) -> DisplayMesh {
    let mut normals = vec![[0.0f64; 3]; vertices.len()];
    for t in triangles {
        let [a, b, c] = t.map(|i| vertices[i as usize]);
        let u = [0, 1, 2].map(|d| b[d] - a[d]);
        let w = [0, 1, 2].map(|d| c[d] - a[d]);
        let n = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        for &i in t {
            for d in 0..3 {
                normals[i as usize][d] += n[d];
            }
        }
    }
    DisplayMesh {
        positions: vertices.iter().flatten().map(|&x| x as f32).collect(),
        normals: normals
            .iter()
            .flat_map(|n| {
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                let k = if len > 0.0 { 1.0 / len } else { 0.0 };
                n.map(|x| (x * k) as f32)
            })
            .collect(),
        indices: triangles.iter().flatten().copied().collect(),
    }
}

impl DesignHistory {
    /// The design at `at`, a position among the records (0 the first,
    /// `frames.len() - 1` the final design; clamped), as the
    /// iso-surface at ½ of the two records either side blended
    /// linearly.
    pub fn frame(&self, at: f64) -> DisplayMesh {
        let Some(last) = self.frames.len().checked_sub(1) else {
            return DisplayMesh::default();
        };
        let at = if at.is_finite() {
            at.clamp(0.0, last as f64)
        } else {
            0.0
        };
        let lo = (at.floor() as usize).min(last);
        let hi = (lo + 1).min(last);
        let w = at - lo as f64;
        let (a, b) = (&self.frames[lo].material, &self.frames[hi].material);
        let mut material = vec![0.0; self.dims[0] * self.dims[1] * self.dims[2]];
        for (k, &c) in self.active.iter().enumerate() {
            let (x, y) = (
                f64::from(a.get(k).copied().unwrap_or(0)),
                f64::from(b.get(k).copied().unwrap_or(0)),
            );
            if let Some(m) = material.get_mut(c as usize) {
                *m = (x + (y - x) * w) / 255.0;
            }
        }
        let surface = grid_surface(
            &MaterialGrid {
                cell_size: self.cell,
                dims: self.dims,
                corner: self.corner,
                material,
            },
            0.5,
        );
        display_mesh(&surface.vertices, &surface.triangles)
    }
}

#[cfg(test)]
#[path = "history_tests.rs"]
mod tests;
