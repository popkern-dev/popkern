//! The band rasterizer: opaque triangles into a G-buffer (eye depth,
//! item, world position, normal), shaded once per pixel; the ground
//! shadow on empty pixels; translucent triangles blended over in their
//! items' far-to-near order; then the edges, depth-tested against the
//! opaque surfaces with a small relative bias.
//!
//! Triangles are clipped to the near plane in clip space and sampled at
//! pixel centres with perspective-correct interpolation; the section
//! discards fragments on its removed side, and a back face seen through
//! the cut is drawn in the flat cap tone so a cut solid reads closed.

use popkern_core::geometry::PortableMath;
use std::ops::Range;

use glam::{Mat4, Vec3, Vec4};

use super::shade;
use super::shadow::ShadowMap;
use super::{FrameScene, SUPERSAMPLE, Tri};

/// What every band shares.
pub(crate) struct Context<'a, 'b> {
    pub scene: &'a FrameScene<'b>,
    pub vp_inverse: Mat4,
    pub width: usize,
    pub height: usize,
    pub shadow: Option<&'a ShadowMap>,
}

/// An edge segment ready to draw.
pub(crate) struct Segment {
    clip: [Vec4; 2],
    world: [Vec3; 2],
    alpha: f32,
}

/// A clip-space vertex with the attributes interpolated across it.
#[derive(Clone, Copy)]
struct Vertex {
    clip: Vec4,
    world: Vec3,
    normal: Vec3,
}

/// Clip a polygon to the near plane (`z ≥ −w`, GL clip space).
fn clip_near(poly: &[Vertex]) -> Vec<Vertex> {
    let inside = |v: &Vertex| v.clip.z + v.clip.w >= 0.0;
    let mut out = Vec::with_capacity(poly.len() + 1);
    for i in 0..poly.len() {
        let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
        if inside(&a) {
            out.push(a);
        }
        if inside(&a) != inside(&b) {
            let (da, db) = (a.clip.z + a.clip.w, b.clip.z + b.clip.w);
            let s = da / (da - db);
            out.push(Vertex {
                clip: a.clip.lerp(b.clip, s),
                world: a.world.lerp(b.world, s),
                normal: a.normal.lerp(b.normal, s),
            });
        }
    }
    out
}

/// A screen-space vertex: pixel position, 1/w, and attributes over w.
#[derive(Clone, Copy)]
struct Screen {
    x: f32,
    y: f32,
    inv_w: f32,
    world: Vec3,
    normal: Vec3,
}

fn to_screen(v: &Vertex, width: usize, height: usize) -> Screen {
    let inv_w = 1.0 / v.clip.w;
    Screen {
        x: (v.clip.x * inv_w * 0.5 + 0.5) * width as f32,
        y: (1.0 - (v.clip.y * inv_w * 0.5 + 0.5)) * height as f32,
        inv_w,
        world: v.world * inv_w,
        normal: v.normal * inv_w,
    }
}

/// A fragment: eye depth (w), world position, interpolated normal.
struct Fragment {
    w: f32,
    world: Vec3,
    normal: Vec3,
}

/// Rasterize `tri`'s fragments inside `rows`, calling `emit(x, y, f)`.
fn rasterize(
    tri: &Tri,
    width: usize,
    height: usize,
    rows: &Range<usize>,
    mut emit: impl FnMut(usize, usize, Fragment),
) {
    let poly: Vec<Vertex> = (0..3)
        .map(|i| Vertex {
            clip: tri.clip[i],
            world: tri.world[i],
            normal: tri.normal[i],
        })
        .collect();
    let poly = if poly.iter().all(|v| v.clip.z + v.clip.w >= 0.0) {
        poly
    } else {
        clip_near(&poly)
    };
    if poly.len() < 3 {
        return;
    }
    let screen: Vec<Screen> = poly.iter().map(|v| to_screen(v, width, height)).collect();
    for k in 1..screen.len() - 1 {
        let (a, b, c) = (screen[0], screen[k], screen[k + 1]);
        let area = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
        if !area.is_finite() || area.abs() < 1e-12 {
            continue;
        }
        let min_x = a.x.min(b.x).min(c.x).floor().max(0.0) as usize;
        let max_x = (a.x.max(b.x).max(c.x).ceil() as isize).min(width as isize - 1);
        let min_y = (a.y.min(b.y).min(c.y).floor().max(rows.start as f32)) as usize;
        let max_y = (a.y.max(b.y).max(c.y).ceil() as isize).min(rows.end as isize - 1);
        if max_x < 0 || max_y < min_y as isize {
            continue;
        }
        let inv_area = 1.0 / area;
        for y in min_y..=max_y as usize {
            let py = y as f32 + 0.5;
            for x in min_x..=max_x as usize {
                let px = x as f32 + 0.5;
                let w0 = ((b.x - px) * (c.y - py) - (c.x - px) * (b.y - py)) * inv_area;
                let w1 = ((c.x - px) * (a.y - py) - (a.x - px) * (c.y - py)) * inv_area;
                let w2 = 1.0 - w0 - w1;
                if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                    continue;
                }
                let inv_w = w0 * a.inv_w + w1 * b.inv_w + w2 * c.inv_w;
                if inv_w <= 0.0 {
                    continue;
                }
                let w = 1.0 / inv_w;
                emit(
                    x,
                    y,
                    Fragment {
                        w,
                        world: (a.world * w0 + b.world * w1 + c.world * w2) * w,
                        normal: (a.normal * w0 + b.normal * w1 + c.normal * w2) * w,
                    },
                );
            }
        }
    }
}

/// Whether the section removes `p`.
fn cut(scene: &FrameScene<'_>, p: Vec3) -> bool {
    scene
        .section
        .is_some_and(|(n, offset)| n.as_vec3().dot(p) > offset as f32)
}

/// The colour of a fragment of item `item` (cap tone for a back face
/// seen through the cut).
fn fragment_color(ctx: &Context<'_, '_>, item: usize, f: &Fragment) -> Vec3 {
    let scene = ctx.scene;
    let it = &scene.items[item];
    let albedo = shade::linear(it.color);
    let eye = scene.camera.eye.as_vec3();
    let view = (eye - f.world).normalize_or_zero();
    let mut n = f.normal.normalize_or_zero();
    if n.dot(view) < 0.0 {
        if scene.section.is_some() {
            return shade::cap(albedo);
        }
        n = -n;
    }
    shade::lit(scene.look, albedo, n, view, it.highlight)
}

/// The background at pixel `(x, y)`: the ground's shadow where the
/// pixel's ray meets the ground plane.
fn background(ctx: &Context<'_, '_>, x: usize, y: usize) -> Vec3 {
    let b = ctx.scene.background;
    let base = shade::linear([
        b[0] as f32 / 255.0,
        b[1] as f32 / 255.0,
        b[2] as f32 / 255.0,
    ]);
    let Some(shadow) = ctx.shadow else {
        return base;
    };
    let ndc_x = (x as f32 + 0.5) / ctx.width as f32 * 2.0 - 1.0;
    let ndc_y = 1.0 - (y as f32 + 0.5) / ctx.height as f32 * 2.0;
    let near = ctx.vp_inverse * Vec4::new(ndc_x, ndc_y, -1.0, 1.0);
    let far = ctx.vp_inverse * Vec4::new(ndc_x, ndc_y, 1.0, 1.0);
    let (near, far) = (near.truncate() / near.w, far.truncate() / far.w);
    let dir = far - near;
    if dir.z.abs() < 1e-12 {
        return base;
    }
    let s = (shadow.ground_z - near.z) / dir.z;
    if s <= 0.0 {
        return base;
    }
    let g = near + dir * s;
    base * (1.0 - 0.16 * shadow.darkness(g))
}

/// Render rows `rows` of the frame: linear colour, row-major.
pub(crate) fn render_band(
    ctx: &Context<'_, '_>,
    opaque: &[Tri],
    translucent: &[Tri],
    edges: &[Segment],
    rows: Range<usize>,
) -> Vec<[f32; 3]> {
    let (width, height) = (ctx.width, ctx.height);
    let n = rows.len() * width;
    let mut depth = vec![f32::INFINITY; n];
    let mut gbuffer: Vec<Option<(usize, Fragment)>> = (0..n).map(|_| None).collect();
    let at = |x: usize, y: usize| (y - rows.start) * width + x;
    for tri in opaque {
        rasterize(tri, width, height, &rows, |x, y, f| {
            let i = at(x, y);
            if f.w < depth[i] && !cut(ctx.scene, f.world) {
                depth[i] = f.w;
                gbuffer[i] = Some((tri.item as usize, f));
            }
        });
    }
    let mut color: Vec<Vec3> = Vec::with_capacity(n);
    for (i, g) in gbuffer.iter().enumerate() {
        let (x, y) = (i % width, rows.start + i / width);
        color.push(match g {
            Some((item, f)) => fragment_color(ctx, *item, f),
            None => background(ctx, x, y),
        });
    }
    for tri in translucent {
        let alpha = ctx.scene.items[tri.item as usize].opacity;
        rasterize(tri, width, height, &rows, |x, y, f| {
            let i = at(x, y);
            if f.w < depth[i] && !cut(ctx.scene, f.world) {
                let c = fragment_color(ctx, tri.item as usize, &f);
                color[i] = color[i] * (1.0 - alpha) + c * alpha;
            }
        });
    }
    draw_edges(ctx, edges, &rows, &depth, &mut color);
    color.into_iter().map(|c| [c.x, c.y, c.z]).collect()
}

/// Every item's edges as clip-space segments (near-clipped).
pub(crate) fn edge_segments(scene: &FrameScene<'_>, vp: Mat4) -> Vec<Segment> {
    let mut out = Vec::new();
    for item in &scene.items {
        if item.opacity <= 0.001 {
            continue;
        }
        let world = item.world.as_mat4();
        for pair in item.edges.as_chunks::<2>().0.iter() {
            let mut p =
                [Vec3::from(pair[0]), Vec3::from(pair[1])].map(|v| world.transform_point3(v));
            let mut c = p.map(|v| vp * v.extend(1.0));
            let d = [c[0].z + c[0].w, c[1].z + c[1].w];
            if d[0] < 0.0 && d[1] < 0.0 {
                continue;
            }
            if d[0] < 0.0 || d[1] < 0.0 {
                let s = d[0] / (d[0] - d[1]);
                let (cm, pm) = (c[0].lerp(c[1], s), p[0].lerp(p[1], s));
                let k = usize::from(d[1] < 0.0);
                c[k] = cm;
                p[k] = pm;
            }
            out.push(Segment {
                clip: c,
                world: p,
                alpha: item.opacity,
            });
        }
    }
    out
}

/// Draw `edges` over `color` where they are not behind the opaque
/// surfaces (`depth`, eye distance).
fn draw_edges(
    ctx: &Context<'_, '_>,
    edges: &[Segment],
    rows: &Range<usize>,
    depth: &[f32],
    color: &mut [Vec3],
) {
    let (width, height) = (ctx.width, ctx.height);
    let half = 0.5 * SUPERSAMPLE as f32;
    for seg in edges {
        let [a, b] = seg.clip.map(|c| {
            let iw = 1.0 / c.w;
            (
                (c.x * iw * 0.5 + 0.5) * width as f32,
                (1.0 - (c.y * iw * 0.5 + 0.5)) * height as f32,
                iw,
            )
        });
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let length = dx.phypot(dy);
        if !length.is_finite() {
            continue;
        }
        let lo = a.1.min(b.1) - half - 1.0;
        let hi = a.1.max(b.1) + half + 1.0;
        if hi < rows.start as f32 || lo > rows.end as f32 {
            continue;
        }
        let steps = (length.ceil() as usize).max(1);
        for s in 0..=steps {
            let t = s as f32 / steps as f32;
            let (x, y) = (a.0 + dx * t, a.1 + dy * t);
            // Perspective-correct depth and position along the segment.
            let inv_w = a.2 + (b.2 - a.2) * t;
            if inv_w <= 0.0 {
                continue;
            }
            let k = t * b.2 / inv_w;
            let w = 1.0 / inv_w;
            let p = seg.world[0].lerp(seg.world[1], k);
            if cut(ctx.scene, p) {
                continue;
            }
            let (x0, x1) = ((x - half).round() as isize, (x + half).round() as isize);
            let (y0, y1) = ((y - half).round() as isize, (y + half).round() as isize);
            for py in y0.max(rows.start as isize)..y1.min(rows.end as isize) {
                for px in x0.max(0)..x1.min(width as isize) {
                    let i = (py as usize - rows.start) * width + px as usize;
                    if w <= depth[i] * 1.002 + 1e-6 * w {
                        color[i] = color[i] * (1.0 - seg.alpha) + shade::EDGE * seg.alpha;
                    }
                }
            }
        }
    }
}
