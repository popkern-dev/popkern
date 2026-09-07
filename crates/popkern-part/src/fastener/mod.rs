//! Fastener recipes over the [`popkern_fastener`] kernel crate
//! (PK-460, design: `docs/decisions/pk460-fastener-subsystem.md`):
//! parametric threads, screws/bolts, nuts, and receiving holes as
//! [`Part`]-level builders.
//!
//! The geometry lives in `popkern-fastener` (operations tier, so the
//! parametric feature tree can replay it); this module wraps it in
//! the code-first facade's recipe idiom and re-exports the shared
//! vocabulary ([`ThreadSpec`], [`standards`], the hole options).
//!
//! ```no_run
//! # use popkern_part::prelude::*;
//! # use popkern_part::fastener::{ThreadSpec, HexBolt, HexNut};
//! # use popkern_core::tolerance::Tolerance;
//! let spec = ThreadSpec::parse("M6")?;              // ISO coarse: 6.0 × 1.0
//! let bolt = HexBolt::new(spec.clone(), 20.0).build(Tolerance::default())?;
//! let nut = HexNut::new(spec).clearance(0.2).build(Tolerance::default())?;
//! # Ok::<(), popkern_part::PartError>(())
//! ```
//!
//! Dimensions are plain numbers in the caller's model units
//! (standard tables are nominal millimetres); every generator builds
//! in the **insertion frame** — axis +Z, origin at the bearing
//! plane, threaded end pointing +Z — and is placed afterwards with
//! the part transforms. Receiving-hole cuts take an arbitrary axis
//! whose direction points into the material, so a fastener built on
//! the same axis mates by construction.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_fastener::builders::{
    self, CountersunkParams, HexBoltParams, NutParams, RodParams, SocketCapParams,
};

pub use popkern_fastener::{
    ClearanceFit, ClearanceHoleOptions, FastenerError, Handedness, HoleDepth, TapOptions,
    ThreadSpec, standards,
};

use crate::error::PartError;
use crate::part::Part;

fn wrap(
    stage: &'static str,
    result: Result<popkern_core::topology::BRep, FastenerError>,
    tolerance: Tolerance,
    summary: String,
) -> Result<Part, PartError> {
    let brep = result.map_err(|source| PartError::Fastener { stage, source })?;
    Part::from_validated(brep, tolerance, "fastener", summary)
}

/// A straight threaded rod: full-length external thread, both ends
/// free-running, optional 45° lead-in chamfer on the +Z tip.
#[derive(Debug, Clone)]
pub struct ThreadedRod {
    params: RodParams,
}

impl ThreadedRod {
    /// A rod of the given thread spec and length (z = 0..length).
    pub fn new(spec: ThreadSpec, length: f64) -> ThreadedRod {
        ThreadedRod {
            params: RodParams {
                spec,
                length,
                tip_chamfer: true,
            },
        }
    }

    /// Enable/disable the tip lead-in chamfer (default on).
    pub fn tip_chamfer(mut self, on: bool) -> ThreadedRod {
        self.params.tip_chamfer = on;
        self
    }

    /// Build the rod.
    ///
    /// # Errors
    ///
    /// [`PartError::Fastener`] carrying the structured
    /// [`FastenerError`] (spec problems, or kernel sweep / boolean /
    /// blend declines with the failing stage named).
    pub fn build(&self, tolerance: Tolerance) -> Result<Part, PartError> {
        wrap(
            "threaded rod",
            builders::threaded_rod(&self.params, tolerance),
            tolerance,
            format!(
                "threaded rod M{} × {}",
                self.params.spec.major_diameter(),
                self.params.length
            ),
        )
    }
}

/// A hex-head bolt (ISO 4017 proportions by default).
#[derive(Debug, Clone)]
pub struct HexBolt {
    params: HexBoltParams,
}

impl HexBolt {
    /// A bolt with under-head length `length` (head at z ≤ 0, thread
    /// tip at z = length).
    pub fn new(spec: ThreadSpec, length: f64) -> HexBolt {
        HexBolt {
            params: HexBoltParams {
                spec,
                length,
                thread_length: None,
                across_flats: None,
                head_height: None,
                head_chamfer: true,
                tip_chamfer: true,
            },
        }
    }

    /// Thread only the last `length` of the shank (default: full
    /// thread up to the head runout).
    pub fn thread_length(mut self, length: f64) -> HexBolt {
        self.params.thread_length = Some(length);
        self
    }

    /// Override the hex width across flats.
    pub fn across_flats(mut self, s: f64) -> HexBolt {
        self.params.across_flats = Some(s);
        self
    }

    /// Override the head height.
    pub fn head_height(mut self, k: f64) -> HexBolt {
        self.params.head_height = Some(k);
        self
    }

    /// Enable/disable the 30° hex corner chamfer (default on).
    pub fn head_chamfer(mut self, on: bool) -> HexBolt {
        self.params.head_chamfer = on;
        self
    }

    /// Enable/disable the tip lead-in chamfer (default on).
    pub fn tip_chamfer(mut self, on: bool) -> HexBolt {
        self.params.tip_chamfer = on;
        self
    }

    /// Build the bolt.
    ///
    /// # Errors
    ///
    /// As [`ThreadedRod::build`].
    pub fn build(&self, tolerance: Tolerance) -> Result<Part, PartError> {
        wrap(
            "hex bolt",
            builders::hex_bolt(&self.params, tolerance),
            tolerance,
            format!(
                "hex bolt M{} × {}",
                self.params.spec.major_diameter(),
                self.params.length
            ),
        )
    }
}

/// A socket-head cap screw (ISO 4762 proportions by default).
#[derive(Debug, Clone)]
pub struct SocketCapScrew {
    params: SocketCapParams,
}

impl SocketCapScrew {
    /// A cap screw with under-head length `length`.
    pub fn new(spec: ThreadSpec, length: f64) -> SocketCapScrew {
        SocketCapScrew {
            params: SocketCapParams {
                spec,
                length,
                thread_length: None,
                tip_chamfer: true,
            },
        }
    }

    /// Thread only the last `length` of the shank.
    pub fn thread_length(mut self, length: f64) -> SocketCapScrew {
        self.params.thread_length = Some(length);
        self
    }

    /// Enable/disable the tip lead-in chamfer (default on).
    pub fn tip_chamfer(mut self, on: bool) -> SocketCapScrew {
        self.params.tip_chamfer = on;
        self
    }

    /// Build the screw.
    ///
    /// # Errors
    ///
    /// As [`ThreadedRod::build`].
    pub fn build(&self, tolerance: Tolerance) -> Result<Part, PartError> {
        wrap(
            "socket cap screw",
            builders::socket_cap_screw(&self.params, tolerance),
            tolerance,
            format!(
                "socket cap screw M{} × {}",
                self.params.spec.major_diameter(),
                self.params.length
            ),
        )
    }
}

/// A countersunk (flat-head) screw: 90° head sinking into the
/// material, top face flush at z = 0. `length` is overall (flush
/// plane to tip).
#[derive(Debug, Clone)]
pub struct CountersunkScrew {
    params: CountersunkParams,
}

impl CountersunkScrew {
    /// A countersunk screw of overall length `length`.
    pub fn new(spec: ThreadSpec, length: f64) -> CountersunkScrew {
        CountersunkScrew {
            params: CountersunkParams {
                spec,
                length,
                thread_length: None,
                tip_chamfer: true,
            },
        }
    }

    /// Thread only the last `length` of the shank.
    pub fn thread_length(mut self, length: f64) -> CountersunkScrew {
        self.params.thread_length = Some(length);
        self
    }

    /// Enable/disable the tip lead-in chamfer (default on).
    pub fn tip_chamfer(mut self, on: bool) -> CountersunkScrew {
        self.params.tip_chamfer = on;
        self
    }

    /// Build the screw.
    ///
    /// # Errors
    ///
    /// As [`ThreadedRod::build`].
    pub fn build(&self, tolerance: Tolerance) -> Result<Part, PartError> {
        wrap(
            "countersunk screw",
            builders::countersunk_screw(&self.params, tolerance),
            tolerance,
            format!(
                "countersunk screw M{} × {}",
                self.params.spec.major_diameter(),
                self.params.length
            ),
        )
    }
}

/// A hex nut (ISO 4032 proportions by default) with a real modeled
/// internal thread, fit clearance, and lead-in chamfers on both
/// faces.
#[derive(Debug, Clone)]
pub struct HexNut {
    params: NutParams,
}

impl HexNut {
    /// A nut for the given thread (z = 0..thickness).
    pub fn new(spec: ThreadSpec) -> HexNut {
        HexNut {
            params: NutParams {
                spec,
                thickness: None,
                across_flats: None,
                clearance: 0.0,
                corner_chamfer: true,
                mouth_chamfer: true,
            },
        }
    }

    /// Radial/flank fit clearance to the nominal external thread
    /// (0 = exact basic profile; 0.15–0.3 typical for FDM printing).
    pub fn clearance(mut self, clearance: f64) -> HexNut {
        self.params.clearance = clearance;
        self
    }

    /// Override the nut thickness.
    pub fn thickness(mut self, m: f64) -> HexNut {
        self.params.thickness = Some(m);
        self
    }

    /// Override the hex width across flats.
    pub fn across_flats(mut self, s: f64) -> HexNut {
        self.params.across_flats = Some(s);
        self
    }

    /// Enable/disable the 30° hex corner chamfers (default on).
    pub fn corner_chamfer(mut self, on: bool) -> HexNut {
        self.params.corner_chamfer = on;
        self
    }

    /// Enable/disable the 45° thread-mouth lead-in chamfers on both
    /// faces (default on).
    pub fn mouth_chamfer(mut self, on: bool) -> HexNut {
        self.params.mouth_chamfer = on;
        self
    }

    /// Build the nut.
    ///
    /// # Errors
    ///
    /// As [`ThreadedRod::build`]; additionally
    /// [`FastenerError::ProfileOverlap`] when the clearance consumes
    /// the thread form.
    pub fn build(&self, tolerance: Tolerance) -> Result<Part, PartError> {
        wrap(
            "hex nut",
            builders::hex_nut(&self.params, tolerance),
            tolerance,
            format!("hex nut M{}", self.params.spec.major_diameter()),
        )
    }
}

impl Part {
    /// Cut a tapped (threaded) hole: a bore at the internal minor
    /// diameter plus a real modeled internal thread with the
    /// requested fit clearance, entered at `origin` along
    /// `direction` (pointing into the material). See
    /// [`popkern_fastener::holes::tapped_hole`].
    ///
    /// # Errors
    ///
    /// [`PartError::Fastener`] for spec/clearance problems; kernel
    /// boolean diagnostics if a cut declines (the part is unchanged).
    pub fn tapped_hole(
        &self,
        spec: &ThreadSpec,
        origin: Point3,
        direction: Vec3,
        depth: HoleDepth,
        options: &TapOptions,
    ) -> Result<Part, PartError> {
        let mut part = wrap(
            "tapped hole",
            popkern_fastener::holes::tapped_hole(
                self.brep(),
                spec,
                origin,
                direction,
                depth,
                options,
                self.tolerance(),
            ),
            self.tolerance(),
            format!("tapped hole M{}", spec.major_diameter()),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Cut a plain clearance hole for a screw of `spec`, with
    /// optional countersink (flush 90° head) or counterbore (sunk
    /// socket head). See [`popkern_fastener::holes::clearance_hole`].
    ///
    /// # Errors
    ///
    /// As [`Part::tapped_hole`].
    pub fn clearance_hole(
        &self,
        spec: &ThreadSpec,
        origin: Point3,
        direction: Vec3,
        depth: HoleDepth,
        options: &ClearanceHoleOptions,
    ) -> Result<Part, PartError> {
        let mut part = wrap(
            "clearance hole",
            popkern_fastener::holes::clearance_hole(
                self.brep(),
                spec,
                origin,
                direction,
                depth,
                options,
                self.tolerance(),
            ),
            self.tolerance(),
            format!("clearance hole M{}", spec.major_diameter()),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }

    /// Cut an external thread onto an existing (cylindrical) shank
    /// region between `z_from` and `z_to` measured along `direction`
    /// from `origin`. See
    /// [`popkern_fastener::holes::external_thread`].
    ///
    /// # Errors
    ///
    /// As [`Part::tapped_hole`].
    pub fn external_thread(
        &self,
        spec: &ThreadSpec,
        origin: Point3,
        direction: Vec3,
        z_from: f64,
        z_to: f64,
    ) -> Result<Part, PartError> {
        let mut part = wrap(
            "external thread",
            popkern_fastener::holes::external_thread(
                self.brep(),
                spec,
                origin,
                direction,
                z_from,
                z_to,
                self.tolerance(),
            ),
            self.tolerance(),
            format!("external thread M{}", spec.major_diameter()),
        )?;
        part.merge_history_from(self);
        Ok(part)
    }
}
