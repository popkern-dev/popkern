//! Fastener recorder verbs (PK-460): `fastener` on [`Model`] (bolt /
//! screw / rod / nut bodies, built in the insertion frame) and
//! `tapped_hole` / `clearance_hole` on [`Body`].
//!
//! Threads record as plain length parameters (major diameter +
//! pitch): the Rhai skin's designation strings (`"M6"`, `"1/4-20"`)
//! resolve to millimetre-magnitude numbers at record time via
//! [`popkern_fastener::ThreadSpec::parse`] and bind in the model
//! unit, so the thread stays sliderable in the parameter table.
//! Author fastener models in a default-unit session with
//! millimetre-magnitude numbers: the kernel's robust scale envelope
//! for helical thread booleans is the ~1..100 range today (the
//! millimetre-in-SI scale campaign is a recorded follow-up).

use popkern_history::feature::Feature;
use popkern_history::feature_fastener::{ClearanceFitDesc, FastenerKind, ThreadDesc};

use crate::body::Body;
use crate::body_recipes::AxisSpec;
use crate::error::ScriptError;
use crate::model::Model;
use crate::record::record_feature;
use crate::scalar::{Scalar, bind_scalar};

/// A thread specification for the recorder: dimensions in the model
/// unit, plus winding and start count.
#[derive(Debug, Clone)]
pub struct ThreadSpecArg {
    /// Thread major (nominal) diameter.
    pub major: Scalar,
    /// Thread pitch.
    pub pitch: Scalar,
    /// Winding sense. `false` = right-hand.
    pub left_handed: bool,
    /// Number of thread starts (≥ 1).
    pub starts: u32,
}

impl ThreadSpecArg {
    /// A right-hand single-start thread.
    pub fn new(major: impl Into<Scalar>, pitch: impl Into<Scalar>) -> ThreadSpecArg {
        ThreadSpecArg {
            major: major.into(),
            pitch: pitch.into(),
            left_handed: false,
            starts: 1,
        }
    }
}

/// Options for [`Model::fastener`].
#[derive(Debug, Clone, Default)]
pub struct FastenerOptions {
    /// Thread only the last `thread_length` of the shank (`None` =
    /// full thread).
    pub thread_length: Option<Scalar>,
    /// Internal fit clearance (nut).
    pub clearance: Option<Scalar>,
    /// Suppress the 45° tip lead-in chamfer.
    pub no_tip_chamfer: bool,
}

/// Options for [`Body::tapped_hole`].
#[derive(Debug, Clone, Default)]
pub struct TapHoleOptions {
    /// Internal fit clearance (`None` = exact basic profile).
    pub clearance: Option<Scalar>,
    /// Suppress the 45° entry chamfer.
    pub no_entry_chamfer: bool,
}

/// Hole sizing for [`Body::clearance_hole`].
#[derive(Debug, Clone, Default)]
pub enum ClearanceFitArg {
    /// ISO 273 close series.
    Close,
    /// ISO 273 normal series.
    #[default]
    Normal,
    /// ISO 273 loose series.
    Loose,
    /// Explicit hole diameter.
    Custom(Scalar),
}

/// Options for [`Body::clearance_hole`].
#[derive(Debug, Clone, Default)]
pub struct ClearanceHoleArgs {
    /// Hole diameter series or explicit diameter.
    pub fit: ClearanceFitArg,
    /// Sink a 90° countersink with this extra radial room.
    pub countersink: Option<Scalar>,
    /// Sink a counterbore with this extra radial/axial room.
    pub counterbore: Option<Scalar>,
    /// Suppress the 45° entry chamfer.
    pub no_entry_chamfer: bool,
}

fn bind_thread(model: &mut Model, thread: &ThreadSpecArg) -> Result<ThreadDesc, ScriptError> {
    let unit = model.unit;
    if thread.starts == 0 {
        return Err(ScriptError::Argument {
            argument: "starts".into(),
            reason: "a thread needs at least one start".into(),
        });
    }
    Ok(ThreadDesc {
        major_diameter: bind_scalar(
            &mut model.tree,
            &mut model.pending_counter,
            "major_diameter",
            &thread.major,
            unit,
        )?,
        pitch: bind_scalar(
            &mut model.tree,
            &mut model.pending_counter,
            "pitch",
            &thread.pitch,
            unit,
        )?,
        left_handed: thread.left_handed,
        starts: thread.starts,
    })
}

impl Model {
    /// Record a fastener body (PK-460): a bolt, screw, threaded rod,
    /// or nut, built in the insertion frame (axis +Z, origin at the
    /// bearing plane) and placed afterwards with the transform verbs.
    /// For a nut, `length` overrides the standard thickness (pass a
    /// literal `0.0` to keep the nominal).
    pub fn fastener(
        &mut self,
        kind: FastenerKind,
        thread: ThreadSpecArg,
        length: impl Into<Scalar>,
        options: FastenerOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = self.unit;
        let thread = bind_thread(self, &thread)?;
        let length = bind_scalar(
            &mut self.tree,
            &mut self.pending_counter,
            "length",
            &length.into(),
            unit,
        )?;
        let thread_length = options
            .thread_length
            .map(|s| {
                bind_scalar(
                    &mut self.tree,
                    &mut self.pending_counter,
                    "thread_length",
                    &s,
                    unit,
                )
            })
            .transpose()?;
        let clearance = options
            .clearance
            .map(|s| {
                bind_scalar(
                    &mut self.tree,
                    &mut self.pending_counter,
                    "clearance",
                    &s,
                    unit,
                )
            })
            .transpose()?;
        let feature = Feature::Fastener {
            kind,
            thread,
            length,
            thread_length,
            clearance,
            tip_chamfer: !options.no_tip_chamfer,
        };
        let id = record_feature(
            &mut self.tree,
            self.tolerance,
            feature,
            label,
            None,
            Vec::new(),
            &mut self.geometry_current,
        )?;
        self.note_source(id);
        let body = self
            .tree
            .feature(id)
            .map(|n| n.body)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        self.active = Some(body);
        Ok(Body::new(body))
    }
}

impl Body {
    /// Cut a tapped (threaded) hole into the body: entered at the
    /// axis origin, threading into the material along the axis
    /// direction. `depth` is a blind depth; `None` cuts through.
    pub fn tapped_hole(
        self,
        model: &mut Model,
        axis: AxisSpec,
        thread: ThreadSpecArg,
        depth: Option<Scalar>,
        options: TapHoleOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let thread = bind_thread(model, &thread)?;
        let mut records = Vec::new();
        let axis = self.axis(model, axis, &mut records)?;
        let depth = depth
            .map(|s| Body::bind(model, "depth", &s, unit))
            .transpose()?;
        let clearance = Body::bind(
            model,
            "clearance",
            &options.clearance.unwrap_or(Scalar::from(0.0)),
            unit,
        )?;
        let feature = Feature::TappedHole {
            axis,
            thread,
            depth,
            clearance,
            entry_chamfer: !options.no_entry_chamfer,
        };
        self.record(model, feature, label, records)
    }

    /// Cut a plain clearance hole for a screw of the given thread,
    /// with optional countersink or counterbore. `depth` is a blind
    /// depth; `None` cuts through.
    pub fn clearance_hole(
        self,
        model: &mut Model,
        axis: AxisSpec,
        thread: ThreadSpecArg,
        depth: Option<Scalar>,
        options: ClearanceHoleArgs,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let thread = bind_thread(model, &thread)?;
        let mut records = Vec::new();
        let axis = self.axis(model, axis, &mut records)?;
        let depth = depth
            .map(|s| Body::bind(model, "depth", &s, unit))
            .transpose()?;
        let fit = match options.fit {
            ClearanceFitArg::Close => ClearanceFitDesc::Close,
            ClearanceFitArg::Normal => ClearanceFitDesc::Normal,
            ClearanceFitArg::Loose => ClearanceFitDesc::Loose,
            ClearanceFitArg::Custom(s) => {
                ClearanceFitDesc::Custom(Body::bind(model, "hole_diameter", &s, unit)?)
            }
        };
        let countersink = options
            .countersink
            .map(|s| Body::bind(model, "countersink", &s, unit))
            .transpose()?;
        let counterbore = options
            .counterbore
            .map(|s| Body::bind(model, "counterbore", &s, unit))
            .transpose()?;
        let feature = Feature::ClearanceHole {
            axis,
            thread,
            fit,
            depth,
            countersink,
            counterbore,
            entry_chamfer: !options.no_entry_chamfer,
        };
        self.record(model, feature, label, records)
    }
}
