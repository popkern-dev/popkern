//! The presentation renderer on a unit cube whose answers are known.

use super::*;

/// A unit cube centred on the origin: flat-shaded faces and its 12 edges.
struct Cube {
    positions: Vec<f32>,
    normals: Vec<f32>,
    indices: Vec<u32>,
    edges: Vec<[f32; 3]>,
}

fn cube() -> Cube {
    let mut c = Cube {
        positions: Vec::new(),
        normals: Vec::new(),
        indices: Vec::new(),
        edges: Vec::new(),
    };
    let faces: [([f32; 3], [f32; 3], [f32; 3]); 6] = [
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
        ([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
        ([0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        ([0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]),
    ];
    for (n, u, v) in faces {
        let base = (c.positions.len() / 3) as u32;
        let corners: Vec<[f32; 3]> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            .iter()
            .map(|&(a, b)| [0, 1, 2].map(|i| 0.5 * (n[i] + a * u[i] + b * v[i])))
            .collect();
        for p in &corners {
            c.positions.extend_from_slice(p);
            c.normals.extend_from_slice(&n);
        }
        c.indices
            .extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        for k in 0..4 {
            c.edges.push(corners[k]);
            c.edges.push(corners[(k + 1) % 4]);
        }
    }
    c
}

fn scene<'a>(cube: &'a Cube, look: Look) -> FrameScene<'a> {
    FrameScene {
        items: vec![FrameItem {
            positions: &cube.positions,
            normals: &cube.normals,
            indices: &cube.indices,
            edges: &cube.edges,
            world: DMat4::IDENTITY,
            color: [0.42, 0.56, 0.75],
            opacity: 1.0,
            highlight: 0.0,
        }],
        camera: FrameCamera {
            eye: DVec3::new(3.0, -4.0, 2.5),
            target: DVec3::ZERO,
            up: DVec3::Z,
            vfov_deg: 30.0,
            near: 0.1,
            far: 20.0,
        },
        look,
        background: [255, 255, 255],
        section: None,
        ground: Some([DVec3::splat(-0.5), DVec3::splat(0.5)]),
    }
}

fn pixel(frame: &[u8], w: u32, x: u32, y: u32) -> [u8; 3] {
    let o = ((y * w + x) * 4) as usize;
    [frame[o], frame[o + 1], frame[o + 2]]
}

const W: u32 = 160;
const H: u32 = 120;

#[test]
fn the_cube_covers_the_centre_on_a_white_background() {
    let cube = cube();
    let frame = render_frame(&scene(&cube, Look::Studio), W, H);
    assert_eq!(frame.len(), (W * H * 4) as usize);
    assert_ne!(
        pixel(&frame, W, W / 2, H / 2),
        [255, 255, 255],
        "centre is the cube"
    );
    assert_eq!(
        pixel(&frame, W, 2, 2),
        [255, 255, 255],
        "corner is the background"
    );
    assert!(
        frame.as_chunks::<4>().0.iter().all(|p| p[3] == 255),
        "opaque output"
    );
}

#[test]
fn renders_are_deterministic() {
    let cube = cube();
    let s = scene(&cube, Look::Toon);
    assert_eq!(render_frame(&s, W, H), render_frame(&s, W, H));
}

#[test]
fn edges_draw_black_and_the_ground_shadow_darkens() {
    let cube = cube();
    let frame = render_frame(&scene(&cube, Look::Studio), W, H);
    let darkest = frame
        .as_chunks::<4>()
        .0
        .iter()
        .map(|p| p[0] as u32 + p[1] as u32 + p[2] as u32)
        .min()
        .unwrap_or(765);
    assert!(darkest < 60, "edge pixels are near black: {darkest}");
    // The ground darkens background pixels (white without it), by at
    // most the shadow's 16 %.
    let mut without = scene(&cube, Look::Studio);
    without.ground = None;
    let plain = render_frame(&without, W, H);
    let mut darkened = 0;
    for (a, b) in frame
        .as_chunks::<4>()
        .0
        .iter()
        .zip(plain.as_chunks::<4>().0.iter())
    {
        if b[..3] == [255, 255, 255] && a[0] < 255 {
            darkened += 1;
            assert!(a[0] >= 205, "a shadow is soft and light: {}", a[0]);
        }
    }
    assert!(darkened > 10, "shadowed pixels: {darkened}");
}

#[test]
fn a_section_removes_the_high_side_and_shows_the_cut() {
    let cube = cube();
    let mut s = scene(&cube, Look::Studio);
    s.ground = None;
    let whole = render_frame(&s, W, H);
    // Remove everything above z = 0.
    s.section = Some((DVec3::Z, 0.0));
    let cut = render_frame(&s, W, H);
    let covered = |f: &[u8]| {
        f.as_chunks::<4>()
            .0
            .iter()
            .filter(|p| p[..3] != [255, 255, 255])
            .count()
    };
    assert!(covered(&cut) < covered(&whole), "less of the cube shows");
    assert!(covered(&cut) > covered(&whole) / 4, "the low half remains");
    // Removing nothing (plane above the cube) changes nothing.
    s.section = Some((DVec3::Z, 2.0));
    assert_eq!(render_frame(&s, W, H), whole);
}

#[test]
fn a_fade_blends_toward_the_background() {
    let cube = cube();
    let mut s = scene(&cube, Look::Studio);
    s.ground = None;
    let solid = pixel(&render_frame(&s, W, H), W, W / 2, H / 2);
    s.items[0].opacity = 0.4;
    let faded = pixel(&render_frame(&s, W, H), W, W / 2, H / 2);
    for c in 0..3 {
        assert!(faded[c] >= solid[c], "lighter toward white");
    }
    assert_ne!(faded, solid);
    s.items[0].opacity = 0.0;
    assert_eq!(
        pixel(&render_frame(&s, W, H), W, W / 2, H / 2),
        [255, 255, 255]
    );
}

#[test]
fn toon_shades_in_at_most_three_tones_per_face_colour() {
    let cube = cube();
    let mut s = scene(&cube, Look::Toon);
    s.ground = None;
    s.items[0].edges = &[];
    let frame = render_frame(&s, W, H);
    // Interior pixels (away from silhouettes, supersample blends) take
    // one of a handful of colours.
    let mut tones = std::collections::BTreeSet::new();
    for y in (H / 2 - 10)..(H / 2 + 10) {
        for x in (W / 2 - 10)..(W / 2 + 10) {
            tones.insert(pixel(&frame, W, x, y));
        }
    }
    assert!(tones.len() <= 6, "toon tones near the centre: {tones:?}");
}
