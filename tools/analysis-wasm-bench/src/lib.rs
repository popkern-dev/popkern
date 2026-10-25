//! The solver built for wasm32, for gates DT1 (a solve gives the same
//! bits natively and in wasm) and E3 (browser solve time).
//!
//! A bench case is a small little-endian binary: the analysis model's
//! meshes and per-face mechanics, so the wasm side needs nothing but the
//! solver. The native tests encode the frozen STEP fixtures; the node
//! and Chrome runners (`run.sh`) decode and solve the same bytes.
//!
//! Layout: `PHB1`, u32 solids, u32 faces; per face u32 solid, f64 area,
//! u8 kind (0 free, 1 clamp, 2 traction, 3 roller), 3 × f64 traction;
//! per solid u32 vertices, u32 triangles, 3 × f64 per vertex, 3 × u32
//! per triangle, u32 face per triangle, f64 exact volume (NaN: none);
//! f64 Young's modulus, f64 Poisson's ratio.

use popkern_analysis::AnalysisOptions;
use popkern_analysis::element::vector_boundary::PieceSupport;
use popkern_analysis::model::{
    AnalysisModel, FaceId, FaceInfo, Provenance, RoleSource, SolidId, SolidMesh, Surface, Units,
};
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis::structural::{ElasticProblem, solve_elastic};

#[cfg(target_arch = "wasm32")]
mod web;

/// File magic.
const MAGIC: &[u8; 4] = b"PHB1";

/// A decoded bench case.
pub struct Case {
    /// The geometry.
    pub model: AnalysisModel,
    /// The elastic problem on it.
    pub problem: ElasticProblem,
}

/// What one solve produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Outcome {
    /// FNV-1a over the bits of every vertex displacement and von Mises
    /// stress, then the dof and iteration counts.
    pub hash: u64,
    /// Active cells.
    pub cells: usize,
    /// Degrees of freedom.
    pub dof: usize,
    /// CG iterations.
    pub iterations: usize,
}

/// Appends a case's bytes for `model` with `mechanics` per face.
pub fn encode(
    model: &AnalysisModel,
    mechanics: &[FaceMechanics],
    youngs: f64,
    poisson: f64,
) -> Vec<u8> {
    let mut b = Vec::new();
    b.extend_from_slice(MAGIC);
    let u32s = |b: &mut Vec<u8>, v: u32| b.extend_from_slice(&v.to_le_bytes());
    let f64s = |b: &mut Vec<u8>, v: f64| b.extend_from_slice(&v.to_le_bytes());
    u32s(&mut b, model.solids.len() as u32);
    u32s(&mut b, model.faces.len() as u32);
    for (f, m) in model.faces.iter().zip(mechanics) {
        u32s(&mut b, f.solid.0);
        f64s(&mut b, f.area);
        let (kind, t) = match m {
            FaceMechanics::Support(PieceSupport::Fixed { .. }) => (1u8, [0.0; 3]),
            FaceMechanics::Traction(t) => (2, *t),
            FaceMechanics::Support(PieceSupport::Roller) => (3, [0.0; 3]),
            _ => (0, [0.0; 3]),
        };
        b.push(kind);
        for c in t {
            f64s(&mut b, c);
        }
    }
    for s in &model.solids {
        u32s(&mut b, s.vertices.len() as u32);
        u32s(&mut b, s.triangles.len() as u32);
        for c in s.vertices.iter().flatten() {
            f64s(&mut b, *c);
        }
        for i in s.triangles.iter().flatten() {
            u32s(&mut b, *i);
        }
        for f in &s.triangle_face {
            u32s(&mut b, f.0);
        }
        f64s(&mut b, s.exact_volume.unwrap_or(f64::NAN));
    }
    f64s(&mut b, youngs);
    f64s(&mut b, poisson);
    b
}

/// Reads the little-endian fields of a case.
struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Reader<'_> {
    fn take<const N: usize>(&mut self) -> Result<[u8; N], String> {
        let end = self.at + N;
        let s = self
            .bytes
            .get(self.at..end)
            .ok_or_else(|| format!("case truncated at byte {}", self.at))?;
        self.at = end;
        let mut out = [0u8; N];
        out.copy_from_slice(s);
        Ok(out)
    }
    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(self.take()?))
    }
    fn f64(&mut self) -> Result<f64, String> {
        Ok(f64::from_le_bytes(self.take()?))
    }
}

/// Decodes a case.
pub fn decode(bytes: &[u8]) -> Result<Case, String> {
    let mut r = Reader { bytes, at: 0 };
    if &r.take::<4>()? != MAGIC {
        return Err("not a bench case".into());
    }
    let (ns, nf) = (r.u32()? as usize, r.u32()? as usize);
    let mut faces = Vec::with_capacity(nf);
    let mut mechanics = Vec::with_capacity(nf);
    for i in 0..nf {
        let solid = SolidId(r.u32()?);
        let area = r.f64()?;
        let [kind] = r.take::<1>()?;
        let t = [r.f64()?, r.f64()?, r.f64()?];
        mechanics.push(match kind {
            1 => FaceMechanics::Support(PieceSupport::Fixed {
                u0: [0.0; 3],
                gradient: [[0.0; 3]; 3],
            }),
            2 => FaceMechanics::Traction(t),
            3 => FaceMechanics::Support(PieceSupport::Roller),
            _ => FaceMechanics::Free,
        });
        faces.push(FaceInfo {
            id: FaceId(i as u32),
            solid,
            name: None,
            surface: Surface::Freeform,
            role: None,
            role_source: RoleSource::Inferred,
            area,
            // Not carried: the bench solves, it does not check models.
            exact_area: None,
            folded_area: None,
        });
    }
    let mut solids = Vec::with_capacity(ns);
    for s in 0..ns {
        let (nv, nt) = (r.u32()? as usize, r.u32()? as usize);
        let mut vertices = Vec::with_capacity(nv);
        for _ in 0..nv {
            vertices.push([r.f64()?, r.f64()?, r.f64()?]);
        }
        let mut triangles = Vec::with_capacity(nt);
        for _ in 0..nt {
            triangles.push([r.u32()?, r.u32()?, r.u32()?]);
        }
        let mut triangle_face = Vec::with_capacity(nt);
        for _ in 0..nt {
            triangle_face.push(FaceId(r.u32()?));
        }
        let v = r.f64()?;
        solids.push(SolidMesh {
            id: SolidId(s as u32),
            name: None,
            vertices,
            triangles,
            triangle_face,
            exact_volume: (!v.is_nan()).then_some(v),
        });
    }
    let (youngs, poisson) = (r.f64()?, r.f64()?);
    Ok(Case {
        model: AnalysisModel {
            units: Units::Millimetre,
            chord_tolerance: 0.01,
            solids,
            faces,
            concave_edges: Vec::new(),
            provenance: Provenance::Popkern {
                document: "wasm-bench".into(),
                revision: 0,
            },
        },
        problem: ElasticProblem {
            youngs,
            poisson,
            faces: mechanics,
        },
    })
}

fn fnv(hash: &mut u64, bytes: &[u8]) {
    for &b in bytes {
        *hash ^= u64::from(b);
        *hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
}

/// Solves a case under a budget of `max_active_cells`.
pub fn run(bytes: &[u8], max_active_cells: usize) -> Result<Outcome, String> {
    let case = decode(bytes)?;
    let options = AnalysisOptions {
        max_active_cells,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&case.model, &case.problem, &options).map_err(|e| e.to_string())?;
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for d in r.vertex_displacement.iter().flatten().flatten() {
        fnv(&mut hash, &d.to_bits().to_le_bytes());
    }
    for v in r.stresses.iter().flat_map(|s| &s.vertices) {
        fnv(&mut hash, &v.von_mises.to_bits().to_le_bytes());
    }
    fnv(&mut hash, &(r.report.dof as u64).to_le_bytes());
    fnv(&mut hash, &(r.report.cg_iterations as u64).to_le_bytes());
    Ok(Outcome {
        hash,
        cells: r.report.cells_active,
        dof: r.report.dof,
        iterations: r.report.cg_iterations,
    })
}
