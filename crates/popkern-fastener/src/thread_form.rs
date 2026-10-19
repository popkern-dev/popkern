//! Thread-form profile math: the 2D cross-sections of the helical
//! cut tools, in absolute `(radial, axial)` coordinates.
//!
//! Model: ISO 68-1 basic profile generalized to any flank angle.
//! With pitch `P`, half flank angle `α`, and fundamental height
//! `H = P / (2·tan α)`: crest flat `P/8`, root flat `P/4`, thread
//! depth `5H/8`. These are mutually consistent for every `α` — the
//! flank slope `dy/dr = tan α` carries the groove half-width from
//! `7P/16` at the major surface to `P/8` at the minor surface over
//! exactly `5H/8` of depth.
//!
//! Both tools extend past the wall they cut by an **embed margin**
//! (`P/16`, flanks continued straight) so every boolean crossing is
//! strictly transversal — no tangent or coincident pairs by
//! construction. Printer-fit clearance is an exact offset of the
//! ridge profile: flanks and tip each move `clearance` along their
//! outward normals and are re-intersected, so the printed gap is the
//! requested clearance *everywhere*, not just radially.

use super::FastenerError;
use super::spec::ThreadSpec;
use popkern_math::PortableMath;

/// A closed convex tool cross-section in absolute
/// `(radial, axial)` coordinates (radial > 0; axial measured from
/// the groove/ridge centreline).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ToolProfile {
    /// Polygon corners, consistently wound, first point not repeated.
    pub points: Vec<(f64, f64)>,
}

#[cfg_attr(not(test), allow(dead_code))] // test-oracle accessors
impl ToolProfile {
    /// Largest axial reach of the profile from its centre. A headed
    /// screw floors its thread tool by this much above the head plane:
    /// floored by its centre only, the groove notched the head's
    /// underside, and every full-thread screw's cut fell back from the
    /// cellular boolean to the graph-merge retry (PK-762: M5x12 5.3 s,
    /// 2.7 s of it the retry).
    pub fn axial_half_width(&self) -> f64 {
        self.points.iter().fold(0.0_f64, |m, p| m.max(p.1.abs()))
    }

    /// Largest radial coordinate of the profile.
    pub fn radial_max(&self) -> f64 {
        self.points.iter().fold(f64::MIN, |m, p| m.max(p.0))
    }

    /// Smallest radial coordinate of the profile.
    pub fn radial_min(&self) -> f64 {
        self.points.iter().fold(f64::MAX, |m, p| m.min(p.0))
    }
}

/// Embed margin: how far a tool reaches past the wall it cuts.
pub(crate) fn embed_margin(spec: &ThreadSpec) -> f64 {
    spec.pitch() / 16.0
}

/// How far below a bearing head (socket, hex) a screw's thread stays
/// at full major diameter: the groove tool's axial half width plus a
/// tenth of a lead of clip clearance (`headed_screw`'s runout). A
/// tapped hole's entry cone must clear that stretch for the screw to
/// seat flush (PK-849): the two builders share this one rule.
pub(crate) fn head_side_full_major(spec: &ThreadSpec) -> Result<f64, FastenerError> {
    Ok(external_groove(spec)?.axial_half_width() + 0.1 * spec.lead())
}

/// The external-thread groove tool: the thread *space* trapezoid,
/// widened past the major surface by the embed margin. Cutting it
/// from a major-diameter shank leaves the ISO basic external thread.
///
/// # Errors
///
/// [`FastenerError::ProfileOverlap`] when adjacent turns of the tool
/// would touch (`lead − axial width at the embed radius ≤ 0`).
pub(crate) fn external_groove(spec: &ThreadSpec) -> Result<ToolProfile, FastenerError> {
    spec.validate()?;
    let p = spec.pitch();
    let tan_a = (spec.flank_angle() / 2.0).ptan();
    let r_major = spec.major_diameter() / 2.0;
    let r_minor = spec.minor_diameter() / 2.0;
    let embed = embed_margin(spec);
    let r_out = r_major + embed;
    let w_major = 7.0 * p / 16.0;
    let w_minor = p / 8.0;
    // The flanks continue past the wall only briefly, then run
    // parallel (constant width) for the rest of the embed reach: the
    // crest flat left between adjacent grooves is P/8, and unlimited
    // flank continuation outside the material would eat the tool's
    // own inter-turn gap at steep flank angles.
    let flare = flank_flare(embed, p / 8.0, tan_a);
    let r_corner = r_major + flare;
    let w_corner = w_major + flare * tan_a;
    self_adjacency_gap(spec, "external groove tool", w_corner)?;
    Ok(ToolProfile {
        points: vec![
            (r_out, w_corner),
            (r_corner, w_corner),
            (r_minor, w_minor),
            (r_minor, -w_minor),
            (r_corner, -w_corner),
            (r_out, -w_corner),
        ],
    })
}

/// How far a cut tool's flanks keep their slope past the wall they
/// cross before turning parallel: half the embed margin, capped so
/// the flare consumes at most half the remaining `flat` between
/// adjacent turns. Strictly positive whenever `flat` is — the corner
/// sits off the wall, so the flank face always crosses the wall in
/// its interior, transversally.
fn flank_flare(embed: f64, flat: f64, tan_a: f64) -> f64 {
    (0.5 * embed).min(0.25 * flat / tan_a.max(1e-9))
}

/// The internal-thread ridge tool: the external thread's *material*
/// trapezoid grown by `clearance` (exact flank-normal offset) and
/// extended inward past the bore wall by the embed margin. Cutting
/// it outward from a bore of [`internal_bore_radius`] leaves the
/// internal thread with the requested clearance to the nominal
/// external thread everywhere.
///
/// # Errors
///
/// [`FastenerError::ProfileOverlap`] when the clearance consumes the
/// ridge tip or the internal crest flat, or adjacent tool turns
/// would touch; [`FastenerError::InvalidSpec`] for a negative
/// clearance or a tool reaching the axis.
pub(crate) fn internal_ridge(
    spec: &ThreadSpec,
    clearance: f64,
) -> Result<ToolProfile, FastenerError> {
    spec.validate()?;
    if !clearance.is_finite() || clearance < 0.0 {
        return Err(FastenerError::InvalidSpec {
            field: "clearance",
            value: clearance,
            requirement: "finite and >= 0",
        });
    }
    let p = spec.pitch();
    let half = spec.flank_angle() / 2.0;
    let (sin_a, cos_a) = half.psin_cos();
    let tan_a = half.ptan();
    let r_major = spec.major_diameter() / 2.0;
    let r_minor = spec.minor_diameter() / 2.0;
    let embed = embed_margin(spec);
    // Grown +axial flank: through the offset image of the nominal
    // base corner (r_minor, 3P/8), sloping down-axial as r grows.
    let flank_r0 = r_minor + clearance * sin_a;
    let flank_w0 = 3.0 * p / 8.0 + clearance * cos_a;
    let flank_w = |r: f64| flank_w0 - tan_a * (r - flank_r0);
    let r_tip = r_major + clearance;
    let w_tip = flank_w(r_tip);
    if w_tip <= 0.0 {
        return Err(FastenerError::ProfileOverlap {
            detail: format!(
                "clearance {clearance:.6} consumes the thread ridge tip \
                 (grown tip half-width {w_tip:.6})"
            ),
            gap: w_tip,
            pitch: p,
            lead: spec.lead(),
        });
    }
    let r_bore = internal_bore_radius(spec, clearance);
    let crest_flat = p - 2.0 * flank_w(r_bore);
    if crest_flat <= 0.0 {
        return Err(FastenerError::ProfileOverlap {
            detail: format!(
                "clearance {clearance:.6} consumes the internal crest flat at \
                 the bore wall (remaining flat {crest_flat:.6})"
            ),
            gap: crest_flat,
            pitch: p,
            lead: spec.lead(),
        });
    }
    let r_base = r_bore - embed;
    if r_base <= 0.0 {
        return Err(FastenerError::InvalidSpec {
            field: "bore radius minus embed margin",
            value: r_base,
            requirement: "> 0 (thread too small for its pitch)",
        });
    }
    // Flanks keep their slope just under the bore wall, then run
    // parallel to the base: continuing them the whole embed reach
    // would let a near-limit clearance self-overlap the tool inside
    // the bore void even while the in-material crest flat survives.
    let flare = flank_flare(embed, crest_flat, tan_a);
    let r_corner = r_bore - flare;
    let w_corner = flank_w(r_corner);
    self_adjacency_gap(spec, "internal ridge tool", w_corner)?;
    Ok(ToolProfile {
        points: vec![
            (r_base, w_corner),
            (r_corner, w_corner),
            (r_tip, w_tip),
            (r_tip, -w_tip),
            (r_corner, -w_corner),
            (r_base, -w_corner),
        ],
    })
}

/// MATING PHASE CONVENTION (PK-476): an internal thread's sweep
/// start shifts DOWN by half a lead relative to the external
/// convention, so an internal and an external thread of one spec,
/// cut with coincident axis placements (or axial offsets of whole
/// leads), MATE with the design clearance. Without it the phase was
/// an accidental byproduct of start margins, and a nominal rod
/// placed one lead into a tapped hole clashed crest-on-ridge
/// (measured: 862 penetrating mesh vertices; the fit window sat at
/// a half-lead offset). Shifting the START only widens the swept
/// margin — turn planning ceils coverage.
pub(crate) fn internal_phase_start(z_start: f64, spec: &ThreadSpec) -> f64 {
    z_start - spec.lead() / 2.0
}

/// Bore radius of an internal thread: the minor radius grown by the
/// fit clearance.
pub(crate) fn internal_bore_radius(spec: &ThreadSpec, clearance: f64) -> f64 {
    spec.minor_diameter() / 2.0 + clearance
}

/// Adjacent turns of one tool body are `lead` apart axially; the
/// widest half-width must leave a positive gap or the swept tool
/// self-intersects.
fn self_adjacency_gap(
    spec: &ThreadSpec,
    tool: &str,
    widest_half_width: f64,
) -> Result<(), FastenerError> {
    let gap = spec.lead() - 2.0 * widest_half_width;
    if gap <= 0.0 {
        return Err(FastenerError::ProfileOverlap {
            detail: format!("adjacent turns of the {tool} would merge"),
            gap,
            pitch: spec.pitch(),
            lead: spec.lead(),
        });
    }
    Ok(())
}
