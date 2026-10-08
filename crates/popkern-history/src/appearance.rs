//! A part's appearance (six-axis-arm §4.5): one colour and finish for
//! every body of the document. Display and exchange data only — it
//! never changes geometry, so replay ignores it.

use serde::{Deserialize, Serialize};

/// Colour and finish.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Appearance {
    /// Linear-free sRGB components, each in `[0, 1]`.
    pub color: [f64; 3],
    /// Metalness in `[0, 1]` (renderers; STEP carries the colour only).
    #[serde(default)]
    pub metal: f64,
    /// Roughness in `[0, 1]`.
    #[serde(default = "default_rough")]
    pub rough: f64,
}

fn default_rough() -> f64 {
    0.5
}

/// Why an appearance was refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[error("appearance {field} = {value} is outside [0, 1]")]
pub struct AppearanceError {
    /// The component (`red`, `green`, `blue`, `metal`, `rough`).
    pub field: &'static str,
    /// Its value.
    pub value: f64,
}

impl Appearance {
    /// An appearance, every component checked to lie in `[0, 1]`.
    ///
    /// # Errors
    ///
    /// [`AppearanceError`] naming the first component out of range (or
    /// not finite).
    pub fn new(color: [f64; 3], metal: f64, rough: f64) -> Result<Appearance, AppearanceError> {
        for (field, value) in [
            ("red", color[0]),
            ("green", color[1]),
            ("blue", color[2]),
            ("metal", metal),
            ("rough", rough),
        ] {
            if !(0.0..=1.0).contains(&value) {
                return Err(AppearanceError { field, value });
            }
        }
        Ok(Appearance {
            color,
            metal,
            rough,
        })
    }

    /// The colour as `f32` RGB (the STEP document's precision).
    pub fn rgb_f32(&self) -> [f32; 3] {
        [
            self.color[0] as f32,
            self.color[1] as f32,
            self.color[2] as f32,
        ]
    }
}
