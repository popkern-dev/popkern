//! The film data model (docs/architecture/presentation-films.md §5): what
//! a film script produces and the page plays. Lengths are SI (document
//! units times the document's unit scale), angles degrees, times seconds.

use serde::Serialize;

/// A whole film: its frame, look, initial camera, captioned steps and
/// the tracks that animate the scene.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Film {
    /// The film's name (its file stem: `assembly_guide`).
    pub name: String,
    /// Frame size in pixels, `[width, height]`.
    pub size: [u32; 2],
    /// Frames per second.
    pub fps: u32,
    /// `"toon"`, `"studio"` or `"cad"`.
    pub look: String,
    /// Background colour, `#rrggbb`.
    pub background: String,
    /// Total length in seconds (the sum of the steps).
    pub duration: f64,
    /// The camera before any camera track.
    pub camera: CameraKey,
    /// Captioned steps, back to back from 0.
    pub steps: Vec<FilmStep>,
    /// Every animation track, in script order.
    pub tracks: Vec<Track>,
}

/// One step of the storyboard.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FilmStep {
    /// Caption title (empty for a hold).
    pub title: String,
    /// Caption subtitle.
    pub subtitle: String,
    /// Start time in seconds.
    pub start: f64,
    /// Length in seconds.
    pub duration: f64,
    /// `"step"`, `"title"` (drawn large) or `"hold"` (no caption).
    pub kind: &'static str,
    /// 1-based script line of the `step` call, when known.
    pub line: Option<usize>,
}

/// One animated change over `[start, start + duration]`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Track {
    /// Absolute start time in seconds.
    pub start: f64,
    /// Length in seconds (0: a jump at `start`).
    pub duration: f64,
    /// The easing curve.
    pub ease: Ease,
    /// 1-based script line of the action, when known.
    pub line: Option<usize>,
    /// What changes.
    pub channel: Channel,
}

/// An easing curve over a track's unit time.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Ease {
    /// Constant speed.
    Linear,
    /// Slow start.
    In,
    /// Slow end.
    Out,
    /// Slow start and end (smoothstep).
    InOut,
}

/// A camera pose: orbit angles about a target, and a zoom relative to
/// fitting the target in frame.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CameraKey {
    /// Degrees about +Z from +X.
    pub azimuth: f64,
    /// Degrees above the XY plane.
    pub elevation: f64,
    /// 1 fits the target; 2 is twice as close.
    pub zoom: f64,
    /// What the camera looks at and fits: an occurrence or body
    /// selector, or `"all"`.
    pub target: String,
}

/// A camera track's end pose: `None` keeps that field as it is when
/// the track starts.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CameraTo {
    /// Degrees about +Z from +X.
    pub azimuth: Option<f64>,
    /// Degrees above the XY plane.
    pub elevation: Option<f64>,
    /// Zoom relative to fitting the target.
    pub zoom: Option<f64>,
    /// The new target.
    pub target: Option<String>,
}

/// A section plane key: an axis preset at a 0..1 position across the
/// model's bounds, or an explicit plane (SI offset).
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SectionKey {
    /// `"x" | "y" | "z"` removing the high side at `at`.
    Axis {
        /// `"x"`, `"y"` or `"z"`.
        axis: String,
        /// 0 removes everything, 1 nothing.
        at: f64,
    },
    /// `normal · p = offset`, removing the side the normal points to.
    Plane {
        /// Unit normal.
        normal: [f64; 3],
        /// Offset in metres.
        offset: f64,
    },
    /// No section.
    Off,
}

/// What a track animates.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Channel {
    /// Turn the camera about its target by `degrees` over the track.
    Orbit {
        /// Degrees added to the azimuth.
        degrees: f64,
    },
    /// Move the camera to a pose.
    Camera {
        /// The pose at the track's end (unset fields keep the pose the
        /// camera has when the track starts).
        to: CameraTo,
    },
    /// Translate a selection by an offset that eases from `from` to `to`.
    Move {
        /// An occurrence path prefix, `"body N"`, or `"all"`.
        target: String,
        /// Offset at the start, metres.
        from: [f64; 3],
        /// Offset at the end, metres.
        to: [f64; 3],
    },
    /// Rotate a selection about an axis from `from` to `to` degrees.
    Turn {
        /// The selection.
        target: String,
        /// A point on the axis, metres.
        origin: [f64; 3],
        /// The axis direction (unit).
        direction: [f64; 3],
        /// Angle at the start.
        from: f64,
        /// Angle at the end.
        to: f64,
    },
    /// Sweep a placement-only assembly parameter.
    Joint {
        /// The parameter name.
        param: String,
        /// Value at the start, in the parameter's display unit.
        from: f64,
        /// Value at the end.
        to: f64,
    },
    /// The viewer's explode factor.
    Explode {
        /// Factor at the start.
        from: f64,
        /// Factor at the end.
        to: f64,
    },
    /// A section plane.
    Section {
        /// The plane at the start.
        from: SectionKey,
        /// The plane at the end.
        to: SectionKey,
        /// `"any"` or `"all"` (one plane: either).
        combine: String,
    },
    /// Opacity of a selection.
    Fade {
        /// The selection.
        target: String,
        /// Opacity at the start (0..1).
        from: f64,
        /// Opacity at the end.
        to: f64,
    },
    /// Show or hide a selection at the track's start.
    Show {
        /// The selection.
        target: String,
        /// Shown from `start` on.
        visible: bool,
    },
    /// Emphasise a selection (eases in over the first quarter and out
    /// over the last).
    Highlight {
        /// The selection.
        target: String,
    },
}
