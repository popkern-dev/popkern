//! [`ThreadSpec`]: the parametric description of one thread —
//! standard designations parse into it, custom sizes construct it
//! directly, and everything downstream (tools, generators, holes)
//! reads only this.

use super::FastenerError;
use super::standards;

/// Thread winding direction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Handedness {
    /// Right-hand thread (tightens clockwise viewed from the head;
    /// the default everywhere).
    Right,
    /// Left-hand thread.
    Left,
}

/// A thread specification: major diameter, pitch, flank angle,
/// handedness, and start count.
///
/// Standard sizes come from [`ThreadSpec::parse`] (`"M6"`,
/// `"M6x0.75"`, `"1/4-20"`, `"#8-32"` — dimensions in nominal
/// millimetres); non-standard sizes from [`ThreadSpec::custom`]. Both
/// produce complete, first-class specs — the standards tables refine
/// defaults, they never gate.
///
/// The profile model is the ISO 68-1 basic form generalized to any
/// flank angle: fundamental triangle height `H = P / (2·tan(α))`
/// (α = half the flank angle), thread depth `5H/8`, crest flat `P/8`,
/// root flat `P/4`.
#[derive(Debug, Clone, PartialEq)]
pub struct ThreadSpec {
    major_diameter: f64,
    pitch: f64,
    flank_angle: f64,
    handedness: Handedness,
    starts: u32,
}

impl ThreadSpec {
    /// The ISO metric / unified flank angle: 60°.
    pub const ISO_FLANK_ANGLE: f64 = std::f64::consts::FRAC_PI_3;

    /// A custom thread: any major diameter and pitch (same length
    /// unit), 60° flank angle, right-handed, single-start. Adjust
    /// with the builder setters.
    pub fn custom(major_diameter: f64, pitch: f64) -> ThreadSpec {
        ThreadSpec {
            major_diameter,
            pitch,
            flank_angle: Self::ISO_FLANK_ANGLE,
            handedness: Handedness::Right,
            starts: 1,
        }
    }

    /// Parse a standard designation. Accepted grammars
    /// (case-insensitive, dimensions produced in millimetres):
    ///
    /// - `"M6"` — ISO metric, coarse pitch from the ISO 261 table
    /// - `"M6x0.75"` — ISO metric with explicit pitch
    /// - `"1/4-20"` — unified: fractional-inch major, threads per inch
    /// - `"#8-32"` — unified numbered size (major = 0.060 + 0.013·N inch)
    ///
    /// # Errors
    ///
    /// [`FastenerError::UnknownDesignation`] naming the accepted
    /// grammars; for a bare metric size missing from the coarse-pitch
    /// table, the error says to give the pitch explicitly.
    pub fn parse(designation: &str) -> Result<ThreadSpec, FastenerError> {
        let text = designation.trim().to_ascii_uppercase();
        let fail = |detail: String| FastenerError::UnknownDesignation {
            designation: designation.to_string(),
            detail,
        };
        let number = |s: &str, what: &str| -> Result<f64, FastenerError> {
            s.parse::<f64>()
                .ok()
                .filter(|v| v.is_finite() && *v > 0.0)
                .ok_or_else(|| fail(format!("'{s}' is not a positive number ({what})")))
        };
        if let Some(rest) = text.strip_prefix('M') {
            let (major, pitch) = match rest.split_once('X') {
                Some((d, p)) => (
                    number(d, "major diameter in mm")?,
                    Some(number(p, "pitch in mm")?),
                ),
                None => (number(rest, "major diameter in mm")?, None),
            };
            let pitch = match pitch {
                Some(p) => p,
                None => standards::iso_coarse_pitch(major).ok_or_else(|| {
                    fail(format!(
                        "M{major} is not in the ISO coarse-pitch table; \
                         give the pitch explicitly as M{major}x<pitch>"
                    ))
                })?,
            };
            return Ok(ThreadSpec::custom(major, pitch));
        }
        if let Some(rest) = text.strip_prefix('#') {
            let (size, tpi) = rest
                .split_once('-')
                .ok_or_else(|| fail("expected '#<size>-<threads per inch>'".into()))?;
            let size = number(size, "unified number size")?;
            let tpi = number(tpi, "threads per inch")?;
            let major_inch = 0.060 + 0.013 * size;
            return Ok(ThreadSpec::custom(
                major_inch * standards::MM_PER_INCH,
                standards::MM_PER_INCH / tpi,
            ));
        }
        if let Some((frac, tpi)) = text.split_once('-')
            && let Some((num, den)) = frac.split_once('/')
        {
            let num = number(num, "inch fraction numerator")?;
            let den = number(den, "inch fraction denominator")?;
            let tpi = number(tpi, "threads per inch")?;
            return Ok(ThreadSpec::custom(
                num / den * standards::MM_PER_INCH,
                standards::MM_PER_INCH / tpi,
            ));
        }
        Err(fail(
            "accepted grammars: 'M<d>' (ISO coarse), 'M<d>x<pitch>', \
             '<a>/<b>-<tpi>' (unified fractional), '#<n>-<tpi>' (unified numbered)"
                .into(),
        ))
    }

    /// Make the thread left-handed.
    pub fn left_handed(mut self) -> ThreadSpec {
        self.handedness = Handedness::Left;
        self
    }

    /// Set the winding direction.
    pub fn with_handedness(mut self, handedness: Handedness) -> ThreadSpec {
        self.handedness = handedness;
        self
    }

    /// Set the number of thread starts (≥ 1). The lead becomes
    /// `pitch × starts`; starts are phased `360°/starts` apart.
    pub fn with_starts(mut self, starts: u32) -> ThreadSpec {
        self.starts = starts;
        self
    }

    /// Set the full included flank angle in radians (default 60°).
    pub fn with_flank_angle(mut self, radians: f64) -> ThreadSpec {
        self.flank_angle = radians;
        self
    }

    /// Major (nominal) diameter.
    pub fn major_diameter(&self) -> f64 {
        self.major_diameter
    }

    /// Pitch: axial crest-to-crest distance of adjacent thread turns
    /// (of any start).
    pub fn pitch(&self) -> f64 {
        self.pitch
    }

    /// Full included flank angle in radians.
    pub fn flank_angle(&self) -> f64 {
        self.flank_angle
    }

    /// Winding direction.
    pub fn handedness(&self) -> Handedness {
        self.handedness
    }

    /// Number of thread starts.
    pub fn starts(&self) -> u32 {
        self.starts
    }

    /// Axial advance per revolution: `pitch × starts`.
    pub fn lead(&self) -> f64 {
        self.pitch * self.starts as f64
    }

    /// Fundamental triangle height `H = P / (2·tan(flank/2))`
    /// (`≈ 0.866·P` at 60°).
    pub fn fundamental_height(&self) -> f64 {
        self.pitch / (2.0 * (self.flank_angle / 2.0).tan())
    }

    /// Basic thread depth `5H/8`: radial distance between the major
    /// and minor surfaces.
    pub fn thread_depth(&self) -> f64 {
        self.fundamental_height() * 5.0 / 8.0
    }

    /// Minor diameter `major − 2·(5H/8)` (the internal thread's bore
    /// diameter at zero clearance).
    pub fn minor_diameter(&self) -> f64 {
        self.major_diameter - 2.0 * self.thread_depth()
    }

    /// Check the numbers describe a realizable thread.
    ///
    /// # Errors
    ///
    /// [`FastenerError::InvalidSpec`] for non-positive or non-finite
    /// fields; [`FastenerError::ProfileOverlap`] when the thread
    /// depth reaches the axis.
    pub fn validate(&self) -> Result<(), FastenerError> {
        let positive = |field: &'static str, value: f64| -> Result<(), FastenerError> {
            if value.is_finite() && value > 0.0 {
                Ok(())
            } else {
                Err(FastenerError::InvalidSpec {
                    field,
                    value,
                    requirement: "finite and > 0",
                })
            }
        };
        positive("major_diameter", self.major_diameter)?;
        positive("pitch", self.pitch)?;
        positive("flank_angle", self.flank_angle)?;
        if self.flank_angle >= std::f64::consts::PI {
            return Err(FastenerError::InvalidSpec {
                field: "flank_angle",
                value: self.flank_angle,
                requirement: "< 180°",
            });
        }
        if self.starts == 0 {
            return Err(FastenerError::InvalidSpec {
                field: "starts",
                value: 0.0,
                requirement: ">= 1",
            });
        }
        let minor = self.minor_diameter();
        if minor <= 0.0 {
            return Err(FastenerError::ProfileOverlap {
                detail: format!(
                    "thread depth {:.6} reaches the axis (minor diameter {minor:.6}); \
                     increase the flank angle or major diameter, or decrease the pitch",
                    self.thread_depth()
                ),
                gap: minor,
                pitch: self.pitch,
                lead: self.lead(),
            });
        }
        Ok(())
    }
}
