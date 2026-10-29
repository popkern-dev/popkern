//! The scene a native film frame draws (cli-render-export.md §4, §5):
//! the session envelope's bodies placed as popkern live places them
//! (`instances.js` `drawList`), coloured as its film look colours them,
//! with the face-boundary edges `edges.js` extracts, the bounds the
//! camera and section fit, the explode directions, and the film camera
//! (`film-player.js` `frameCamera`).

use popkern_core::geometry::PortableMath;
use std::collections::HashMap;

use glam::{DMat4, DVec3};
use popkern_live::{AssemblyView, BodyMesh, Envelope};
use popkern_rhai::film::state::{CameraPose, ItemState};

/// `shading.js` PALETTE: a part's colour by its palette key.
const PALETTE: [u32; 8] = [
    0x6c8ebf, 0x8fb56c, 0xd9a25f, 0xb87ab8, 0x66b3b3, 0xc97a7a, 0x9c9c6c, 0x7a8fc9,
];
/// `viewer.js` CAD_FACE.
const CAD_FACE: u32 = 0xb4b8be;

fn rgb(hex: u32) -> [f32; 3] {
    [
        ((hex >> 16) & 255) as f32 / 255.0,
        ((hex >> 8) & 255) as f32 / 255.0,
        (hex & 255) as f32 / 255.0,
    ]
}

/// One draw: a body at an occurrence's placement.
#[derive(Debug, Clone)]
pub struct Draw {
    /// Index into the envelope's bodies.
    pub body: usize,
    /// The film selector (`"body N"` or the occurrence path).
    pub selector: String,
    /// Placement before any film motion.
    pub base: DMat4,
    /// Shown by the document.
    pub visible: bool,
    /// sRGB colour in the film's look.
    pub color: [f32; 3],
}

/// Everything static about a film's scene, and its per-body edges.
pub struct Scene {
    /// The draws.
    pub draws: Vec<Draw>,
    /// Face-boundary segments per body (two points each).
    pub edges: Vec<Vec<[f32; 3]>>,
    /// Unique selectors in draw order (the film state's paths).
    pub selectors: Vec<String>,
    /// The visible model's bounds.
    pub bounds: [DVec3; 2],
    /// Bounds per camera target (every path prefix, and `"all"`).
    pub targets: HashMap<String, [DVec3; 2]>,
}

/// A body's colour in a look: its appearance, else the toon palette or
/// the CAD grey (`viewer.js` `material`).
fn color(body: &BodyMesh, key: usize, toon: bool) -> [f32; 3] {
    match &body.appearance {
        Some(a) => a.color.map(|c| c as f32),
        None if toon => rgb(PALETTE[key % PALETTE.len()]),
        None => rgb(CAD_FACE),
    }
}

/// The draws of `bodies` placed by `assembly` (`drawList`).
pub fn draws(bodies: &[BodyMesh], assembly: Option<&AssemblyView>, toon: bool) -> Vec<Draw> {
    match assembly {
        None => bodies
            .iter()
            .enumerate()
            .map(|(i, b)| Draw {
                body: i,
                selector: format!("body {i}"),
                base: DMat4::IDENTITY,
                visible: true,
                color: color(b, i, toon),
            })
            .collect(),
        Some(assembly) => {
            let mut out = Vec::new();
            for instance in assembly.instances.iter().filter(|i| i.is_part) {
                let Some(definition) = assembly.definitions.get(instance.definition as usize)
                else {
                    continue;
                };
                let base = DMat4::from_cols_array(&instance.world.map(f64::from));
                for &body in &definition.bodies {
                    let Some(mesh) = bodies.get(body) else {
                        continue;
                    };
                    out.push(Draw {
                        body,
                        selector: instance.path.clone(),
                        base,
                        visible: instance.visible,
                        color: color(mesh, instance.definition as usize, toon),
                    });
                }
            }
            out
        }
    }
}

/// A welded mesh edge: one of its mesh edges (two vertex indices) and
/// how many times each face uses it.
type EdgeUse = (u32, u32, Vec<(u32, u32)>);

/// The face-boundary segments of `body` (`edges.js` `boundaryEdges`):
/// mesh edges one face uses an odd number of times, after welding at a
/// millionth of the body's diagonal.
pub fn boundary_edges(body: &BodyMesh) -> Vec<[f32; 3]> {
    let p = &body.positions;
    if body.indices.is_empty() || p.len() < 3 {
        return Vec::new();
    }
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for (i, v) in p.iter().enumerate() {
        lo[i % 3] = lo[i % 3].min(*v);
        hi[i % 3] = hi[i % 3].max(*v);
    }
    let diagonal =
        ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
    let tol = if diagonal > 0.0 {
        diagonal * 1e-6
    } else {
        1e-6
    };
    let mut cells: HashMap<(i64, i64, i64), u32> = HashMap::new();
    let weld: Vec<u32> = (0..p.len() / 3)
        .map(|v| {
            let key = (
                (p[3 * v] / tol).round() as i64,
                (p[3 * v + 1] / tol).round() as i64,
                (p[3 * v + 2] / tol).round() as i64,
            );
            let next = cells.len() as u32;
            *cells.entry(key).or_insert(next)
        })
        .collect();
    // Welded edge -> one mesh edge and its use count per face, in first
    // appearance order.
    let mut index: HashMap<(u32, u32), usize> = HashMap::new();
    let mut uses: Vec<EdgeUse> = Vec::new();
    for (t, tri) in body.indices.as_chunks::<3>().0.iter().enumerate() {
        let face = body.face_ids.get(t).copied().unwrap_or(0);
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let (wa, wb) = (weld[a as usize], weld[b as usize]);
            if wa == wb {
                continue;
            }
            let key = (wa.min(wb), wa.max(wb));
            let slot = *index.entry(key).or_insert_with(|| {
                uses.push((a, b, Vec::new()));
                uses.len() - 1
            });
            let faces = &mut uses[slot].2;
            match faces.iter_mut().find(|(f, _)| *f == face) {
                Some((_, n)) => *n += 1,
                None => faces.push((face, 1)),
            }
        }
    }
    let point = |v: u32| {
        let v = v as usize;
        [p[3 * v], p[3 * v + 1], p[3 * v + 2]]
    };
    uses.iter()
        .filter(|(_, _, faces)| faces.iter().any(|(_, n)| n % 2 == 1))
        .flat_map(|&(a, b, _)| [point(a), point(b)])
        .collect()
}

/// `body`'s bounds placed by `m` (the eight corners transformed).
fn placed_bounds(body: &BodyMesh, m: &DMat4) -> [DVec3; 2] {
    let [lo, hi] = body.bounds;
    let mut out = [DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY)];
    for c in 0..8 {
        let p = DVec3::new(
            if c & 1 != 0 { hi[0] } else { lo[0] },
            if c & 2 != 0 { hi[1] } else { lo[1] },
            if c & 4 != 0 { hi[2] } else { lo[2] },
        );
        let q = m.transform_point3(p);
        out[0] = out[0].min(q);
        out[1] = out[1].max(q);
    }
    out
}

fn union(a: Option<[DVec3; 2]>, b: [DVec3; 2]) -> [DVec3; 2] {
    match a {
        Some(a) => [a[0].min(b[0]), a[1].max(b[1])],
        None => b,
    }
}

impl Scene {
    /// The scene of `envelope` in a look (`toon` picks the palette).
    pub fn new(envelope: &Envelope, toon: bool) -> Option<Scene> {
        let draws = draws(&envelope.bodies, envelope.assembly.as_ref(), toon);
        let edges = envelope.bodies.iter().map(boundary_edges).collect();
        let mut selectors: Vec<String> = Vec::new();
        let mut bounds = None;
        let mut targets: HashMap<String, [DVec3; 2]> = HashMap::new();
        for d in &draws {
            if !selectors.contains(&d.selector) {
                selectors.push(d.selector.clone());
            }
            let b = placed_bounds(&envelope.bodies[d.body], &d.base);
            if d.visible {
                bounds = Some(union(bounds, b));
            }
            let mut names: Vec<String> = if d.selector.starts_with("body ") {
                vec![d.selector.clone()]
            } else {
                let parts: Vec<&str> = d.selector.split('/').collect();
                (1..=parts.len()).map(|k| parts[..k].join("/")).collect()
            };
            names.push("all".into());
            for name in names {
                let next = union(targets.get(&name).copied(), b);
                targets.insert(name, next);
            }
        }
        Some(Scene {
            draws,
            edges,
            selectors,
            bounds: bounds?,
            targets,
        })
    }

    /// Explode directions per top-level occurrence (`explodePlan`), from
    /// the draws' current placements.
    pub fn explode_directions(&self, envelope: &Envelope) -> HashMap<String, DVec3> {
        let mut boxes: Vec<(String, [DVec3; 2])> = Vec::new();
        for d in self
            .draws
            .iter()
            .filter(|d| !d.selector.starts_with("body "))
        {
            let top = d.selector.split('/').next().unwrap_or("").to_string();
            let b = placed_bounds(&envelope.bodies[d.body], &d.base);
            match boxes.iter_mut().find(|(t, _)| *t == top) {
                Some((_, existing)) => *existing = union(Some(*existing), b),
                None => boxes.push((top, b)),
            }
        }
        let mut out = HashMap::new();
        if boxes.len() < 2 {
            return out;
        }
        let whole = boxes.iter().fold(None, |acc, (_, b)| Some(union(acc, *b)));
        let Some([lo, hi]) = whole else {
            return out;
        };
        let centre = (lo + hi) * 0.5;
        let size = (hi - lo).length();
        let mut axial = Vec::new();
        for (top, [blo, bhi]) in &boxes {
            let c = (*blo + *bhi) * 0.5;
            let radial = DVec3::new(c.x - centre.x, c.y - centre.y, 0.0);
            if radial.length() > 1e-6 * size {
                let r = radial.normalize() * (size * 0.6);
                out.insert(top.clone(), DVec3::new(r.x, r.y, (c.z - centre.z) * 0.6));
            } else {
                axial.push((top.clone(), c.z));
            }
        }
        axial.sort_by(|a, b| a.1.total_cmp(&b.1));
        let n = axial.len();
        for (k, (top, _)) in axial.into_iter().enumerate() {
            let rank = if n > 1 {
                k as f64 / (n - 1) as f64 - 0.5
            } else {
                0.0
            };
            out.insert(top, DVec3::new(0.0, 0.0, rank * size * 1.2));
        }
        out
    }

    /// A draw's world matrix under a film item state and the explode
    /// factor: turns in order, then the offset (and explode), over its
    /// base placement.
    pub fn placement(
        draw: &Draw,
        item: &ItemState,
        explode: f64,
        directions: &HashMap<String, DVec3>,
    ) -> DMat4 {
        let mut m = DMat4::IDENTITY;
        for turn in &item.turns {
            let o = DVec3::from(turn.origin);
            let axis = DVec3::from(turn.direction).normalize_or_zero();
            if axis == DVec3::ZERO {
                continue;
            }
            let rot = DMat4::from_axis_angle(axis, turn.degrees.to_radians());
            m = DMat4::from_translation(o) * rot * DMat4::from_translation(-o) * m;
        }
        let top = draw.selector.split('/').next().unwrap_or("");
        let ex = directions.get(top).map_or(DVec3::ZERO, |d| *d * explode);
        m = DMat4::from_translation(DVec3::from(item.offset) + ex) * m;
        m * draw.base
    }

    /// The film camera for `pose` at `aspect` (`frameCamera`): the
    /// targets' bounding spheres blended, fitted in the narrower of the
    /// film's 30° vertical field and the matching horizontal one.
    pub fn camera(
        &self,
        pose: &CameraPose,
        aspect: f64,
    ) -> Option<popkern_viz::film_frame::FrameCamera> {
        let mut centre = DVec3::ZERO;
        let mut radius = 0.0;
        let mut weight = 0.0;
        for (name, w) in &pose.targets {
            let Some([lo, hi]) = self.targets.get(name).or_else(|| self.targets.get("all")) else {
                continue;
            };
            centre += (*lo + *hi) * 0.5 * *w;
            radius += w * (*hi - *lo).length() / 2.0;
            weight += w;
        }
        if weight <= 0.0 {
            return None;
        }
        centre /= weight;
        let radius = (radius / weight).max(1e-9);
        let base = 30f64.to_radians();
        let hfov = 2.0 * ((base / 2.0).ptan() * aspect).patan();
        let fit = base.min(hfov);
        let distance = radius / (fit / 2.0).psin() / pose.zoom;
        let (az, el) = (pose.azimuth.to_radians(), pose.elevation.to_radians());
        Some(popkern_viz::film_frame::FrameCamera {
            eye: centre
                + DVec3::new(
                    distance * el.pcos() * az.pcos(),
                    distance * el.pcos() * az.psin(),
                    distance * el.psin(),
                ),
            target: centre,
            up: DVec3::Z,
            vfov_deg: 30.0,
            near: (distance * 1e-3).max(distance - 4.0 * radius),
            far: distance + 6.0 * radius,
        })
    }
}
