//! Native film rendering (cli-render-export.md §5): a film's frames
//! drawn by [`popkern_viz::film_frame`] from the session's envelope, the
//! state at each frame from [`popkern_rhai::film::state`] (the page's
//! `film-state.js`, ported), joints set through the session as the page
//! sets them, captions drawn on top, and the frames handed to an
//! encoder. No browser, no GPU.

use std::collections::HashMap;

use glam::DVec3;
use popkern_live::{Envelope, Session};
use popkern_rhai::film::state::{camera_at, section_at, state_at};
use popkern_rhai::film::{Channel, Film};
use popkern_viz::film_frame::{FrameItem, FrameScene, Look, render_frame};

use crate::film_caption::{Fonts, draw_caption};
use crate::film_scene::Scene;

/// The film-less turntable (cli-render-export.md §5): the document's
/// name over a first eighth of a turn, then the rest, 12 s in all.
pub fn turntable_script(title: &str) -> String {
    let title = serde_json::to_string(title).unwrap_or_else(|_| "\"\"".into());
    format!(
        "let film = storyboard(#{{ look: \"studio\", camera: #{{ azimuth: -35.0, elevation: 25.0, zoom: 1.0 }} }});\n\
         film.title({title}, \"\", 2.5, |s| {{ s.orbit(45.0); }});\n\
         film.step(\"\", \"\", 9.5, |s| {{ s.orbit(315.0); }});\n\
         film\n"
    )
}

/// A one-frame film holding a view (`popkern render`).
pub fn still_script(
    look: &str,
    background: &str,
    azimuth: f64,
    elevation: f64,
    zoom: f64,
) -> String {
    let (look, background) = (
        serde_json::to_string(look).unwrap_or_default(),
        serde_json::to_string(background).unwrap_or_default(),
    );
    format!(
        "let film = storyboard(#{{ look: {look}, background: {background}, \
         camera: #{{ azimuth: {azimuth:?}, elevation: {elevation:?}, zoom: {zoom:?} }} }});\n\
         film.hold(1.0);\n\
         film\n"
    )
}

/// A readable title from a file or folder stem (`boss_plate` → "Boss plate").
pub fn title_of(stem: &str) -> String {
    let spaced = stem.replace(['_', '-'], " ");
    let mut chars = spaced.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().chain(chars).collect(),
        None => String::new(),
    }
}

fn background_rgb(hex: &str) -> [u8; 3] {
    let n = u32::from_str_radix(hex.trim_start_matches('#'), 16).unwrap_or(0xffffff);
    [(n >> 16) as u8, (n >> 8) as u8, n as u8]
}

/// A film ready to draw frames of.
pub struct NativeFilm<'s> {
    session: &'s mut Session,
    envelope: Envelope,
    film: Film,
    scene: Scene,
    directions: HashMap<String, DVec3>,
    joints: Vec<(String, f64)>,
    look: Look,
    background: [u8; 3],
    fonts: Option<Fonts>,
}

impl<'s> NativeFilm<'s> {
    /// The film over the session's built document (`envelope` is its
    /// build). Fails for what the native player cannot show.
    pub fn new(
        session: &'s mut Session,
        envelope: Envelope,
        film: Film,
        fonts: Option<Fonts>,
    ) -> Result<NativeFilm<'s>, String> {
        if film
            .tracks
            .iter()
            .any(|t| matches!(t.channel, Channel::Evolve { .. }))
        {
            return Err(
                "this film evolves a generative design (s.evolve), which only the \
                        browser renders: use --browser"
                    .into(),
            );
        }
        let toon = film.look == "toon";
        let scene = Scene::new(&envelope, toon).ok_or("nothing to render: no visible body")?;
        let directions = scene.explode_directions(&envelope);
        Ok(NativeFilm {
            session,
            look: Look::named(&film.look),
            background: background_rgb(&film.background),
            envelope,
            film,
            scene,
            directions,
            joints: Vec::new(),
            fonts,
        })
    }

    /// The film.
    pub fn film(&self) -> &Film {
        &self.film
    }

    /// Set the joints a frame asks for (once per distinct value set) and
    /// re-read the placements they move.
    fn apply_joints(&mut self, joints: &[(String, f64)]) -> Result<(), String> {
        let same = joints.len() == self.joints.len()
            && joints
                .iter()
                .zip(&self.joints)
                .all(|(a, b)| a.0 == b.0 && (a.1 - b.1).abs() <= 1e-9);
        if joints.is_empty() || same {
            return Ok(());
        }
        let reply = self.session.set_joints(joints);
        if !reply.outcome.is_served() {
            return Err(format!(
                "setting joints failed: {}",
                serde_json::to_string(&reply.outcome).unwrap_or_default()
            ));
        }
        self.joints = joints.to_vec();
        if let Some(assembly) = &reply.assembly {
            let draws = crate::film_scene::draws(
                &self.envelope.bodies,
                Some(assembly),
                self.film.look == "toon",
            );
            if draws.len() == self.scene.draws.len() {
                self.scene.draws = draws;
                self.directions = self.scene.explode_directions(&self.envelope);
            }
        }
        Ok(())
    }

    /// The frame at `t`, `width × height` RGBA.
    pub fn frame(&mut self, t: f64, width: u32, height: u32) -> Result<Vec<u8>, String> {
        let t = t.clamp(0.0, self.film.duration);
        let state = state_at(&self.film, t, &self.scene.selectors);
        self.apply_joints(&state.joints)?;
        let items: HashMap<&str, &popkern_rhai::film::state::ItemState> =
            state.items.iter().map(|(p, i)| (p.as_str(), i)).collect();
        let mut frame_items = Vec::with_capacity(self.scene.draws.len());
        for draw in &self.scene.draws {
            let Some(item) = items.get(draw.selector.as_str()) else {
                continue;
            };
            if !(draw.visible && item.visible && item.opacity > 0.001) {
                continue;
            }
            let body = &self.envelope.bodies[draw.body];
            frame_items.push(FrameItem {
                positions: &body.positions,
                normals: &body.normals,
                indices: &body.indices,
                edges: &self.scene.edges[draw.body],
                world: Scene::placement(draw, item, state.explode, &self.directions),
                color: draw.color,
                opacity: item.opacity.clamp(0.0, 1.0) as f32,
                highlight: item.highlight as f32,
            });
        }
        let pose = camera_at(&self.film, t);
        let camera = self
            .scene
            .camera(&pose, width as f64 / height.max(1) as f64)
            .ok_or("the camera's target has no bounds")?;
        let [lo, hi] = self.scene.bounds;
        let section = state.section.as_ref().and_then(|s| {
            section_at(s, &[lo.to_array(), hi.to_array()])
                .map(|p| (DVec3::from(p.normal), p.offset))
        });
        let scene = FrameScene {
            items: frame_items,
            camera,
            look: self.look,
            background: self.background,
            section,
            ground: Some(self.scene.bounds),
        };
        let mut rgba = render_frame(&scene, width, height);
        if let (Some(caption), Some(fonts)) = (&state.caption, &self.fonts) {
            draw_caption(
                &mut rgba,
                width as usize,
                height as usize,
                caption,
                &self.film.background,
                fonts,
            );
        }
        Ok(rgba)
    }
}
