//! Standard fastener dimension tables (nominal millimetres) with
//! documented derivation formulas for sizes outside the tables — a
//! custom spec always gets a complete fastener; the tables refine,
//! they never gate.
//!
//! Values are the nominal product dimensions of the common standards
//! (ISO 261 coarse pitches, ISO 4017 hex heads, ISO 4762 socket
//! heads, ISO 10642 countersunk heads, ISO 4032 nuts, ISO 273
//! clearance holes). They size printable, correctly-proportioned
//! parts; they are not a substitute for a toleranced standards
//! document.

/// Millimetres per inch (exact).
pub const MM_PER_INCH: f64 = 25.4;

/// ISO 261 coarse pitches by major diameter (mm).
const ISO_COARSE: &[(f64, f64)] = &[
    (1.6, 0.35),
    (2.0, 0.4),
    (2.5, 0.45),
    (3.0, 0.5),
    (3.5, 0.6),
    (4.0, 0.7),
    (5.0, 0.8),
    (6.0, 1.0),
    (8.0, 1.25),
    (10.0, 1.5),
    (12.0, 1.75),
    (14.0, 2.0),
    (16.0, 2.0),
    (18.0, 2.5),
    (20.0, 2.5),
    (22.0, 2.5),
    (24.0, 3.0),
    (27.0, 3.0),
    (30.0, 3.5),
    (36.0, 4.0),
];

/// (major, across flats, head height) — ISO 4017 hex heads.
const HEX_HEAD: &[(f64, f64, f64)] = &[
    (3.0, 5.5, 2.0),
    (4.0, 7.0, 2.8),
    (5.0, 8.0, 3.5),
    (6.0, 10.0, 4.0),
    (8.0, 13.0, 5.3),
    (10.0, 16.0, 6.4),
    (12.0, 18.0, 7.5),
    (14.0, 21.0, 8.8),
    (16.0, 24.0, 10.0),
    (20.0, 30.0, 12.5),
    (24.0, 36.0, 15.0),
];

/// (major, head diameter, head height, socket across flats) — ISO 4762.
const SOCKET_HEAD: &[(f64, f64, f64, f64)] = &[
    (3.0, 5.5, 3.0, 2.5),
    (4.0, 7.0, 4.0, 3.0),
    (5.0, 8.5, 5.0, 4.0),
    (6.0, 10.0, 6.0, 5.0),
    (8.0, 13.0, 8.0, 6.0),
    (10.0, 16.0, 10.0, 8.0),
    (12.0, 18.0, 12.0, 10.0),
    (16.0, 24.0, 16.0, 14.0),
    (20.0, 30.0, 20.0, 17.0),
    (24.0, 36.0, 24.0, 19.0),
];

/// (major, head diameter) — ISO 10642 countersunk heads (90°).
const COUNTERSUNK_HEAD: &[(f64, f64)] = &[
    (3.0, 6.7),
    (4.0, 9.0),
    (5.0, 11.2),
    (6.0, 13.4),
    (8.0, 17.9),
    (10.0, 22.4),
    (12.0, 26.9),
];

/// (major, nut thickness) — ISO 4032 style-1 hex nuts (across flats
/// as [`hex_head`]).
const HEX_NUT: &[(f64, f64)] = &[
    (3.0, 2.4),
    (4.0, 3.2),
    (5.0, 4.7),
    (6.0, 5.2),
    (8.0, 6.8),
    (10.0, 8.4),
    (12.0, 10.8),
    (14.0, 12.8),
    (16.0, 14.8),
    (20.0, 18.0),
    (24.0, 21.5),
];

/// (major, close, normal, loose) — ISO 273 clearance hole diameters.
const CLEARANCE_HOLE: &[(f64, f64, f64, f64)] = &[
    (1.6, 1.7, 1.8, 2.0),
    (2.0, 2.2, 2.4, 2.6),
    (2.5, 2.7, 2.9, 3.1),
    (3.0, 3.2, 3.4, 3.6),
    (4.0, 4.3, 4.5, 4.8),
    (5.0, 5.3, 5.5, 5.8),
    (6.0, 6.4, 6.6, 7.0),
    (8.0, 8.4, 9.0, 10.0),
    (10.0, 10.5, 11.0, 12.0),
    (12.0, 13.0, 13.5, 14.5),
    (14.0, 15.0, 15.5, 16.5),
    (16.0, 17.0, 17.5, 18.5),
    (20.0, 21.0, 22.0, 24.0),
    (24.0, 25.0, 26.0, 28.0),
];

/// How much larger than the screw a plain clearance hole is drilled.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClearanceFit {
    /// ISO 273 close series (≈ 1.06·d off-table).
    Close,
    /// ISO 273 normal (medium) series (≈ 1.10·d off-table).
    Normal,
    /// ISO 273 loose (coarse) series (≈ 1.20·d off-table).
    Loose,
    /// An explicit hole diameter in model units.
    Custom(f64),
}

/// Table row match: majors are catalogue values, compared snapped.
fn row<T: Copy>(table: &[(f64, T)], major: f64) -> Option<T> {
    table
        .iter()
        .find(|(d, _)| (d - major).abs() <= 1e-6)
        .map(|&(_, v)| v)
}

/// ISO coarse pitch for a metric major diameter, when tabulated.
pub fn iso_coarse_pitch(major: f64) -> Option<f64> {
    row(ISO_COARSE, major)
}

/// Hex head proportions (bolt heads and nuts share across-flats).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HexHeadDims {
    /// Width across flats of the hexagon.
    pub across_flats: f64,
    /// Head height (axial).
    pub height: f64,
}

/// Hex head dimensions: ISO 4017 nominal, else derived as
/// `across_flats = 1.6·d`, `height = 0.65·d`.
pub fn hex_head(major: f64) -> HexHeadDims {
    match HEX_HEAD.iter().find(|(d, _, _)| (d - major).abs() <= 1e-6) {
        Some(&(_, s, k)) => HexHeadDims {
            across_flats: s,
            height: k,
        },
        None => HexHeadDims {
            across_flats: 1.6 * major,
            height: 0.65 * major,
        },
    }
}

/// Socket (cap) head proportions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SocketHeadDims {
    /// Cylindrical head diameter.
    pub diameter: f64,
    /// Head height (axial).
    pub height: f64,
    /// Hex socket width across flats.
    pub socket_across_flats: f64,
}

/// Socket head dimensions: ISO 4762 nominal, else derived as
/// `diameter = 1.5·d`, `height = d`, `socket = 0.75·d`.
pub fn socket_head(major: f64) -> SocketHeadDims {
    match SOCKET_HEAD
        .iter()
        .find(|(d, _, _, _)| (d - major).abs() <= 1e-6)
    {
        Some(&(_, dk, k, s)) => SocketHeadDims {
            diameter: dk,
            height: k,
            socket_across_flats: s,
        },
        None => SocketHeadDims {
            diameter: 1.5 * major,
            height: major,
            socket_across_flats: 0.75 * major,
        },
    }
}

/// Countersunk head diameter (90° included angle): ISO 10642
/// nominal, else derived as `2.24·d`.
pub fn countersunk_head_diameter(major: f64) -> f64 {
    row(COUNTERSUNK_HEAD, major).unwrap_or(2.24 * major)
}

/// Hex nut thickness: ISO 4032 nominal, else derived as `0.85·d`.
pub fn hex_nut_thickness(major: f64) -> f64 {
    row(HEX_NUT, major).unwrap_or(0.85 * major)
}

/// Clearance hole diameter for a screw of the given major diameter:
/// ISO 273 nominal, else derived per the [`ClearanceFit`] ratios.
pub fn clearance_diameter(major: f64, fit: ClearanceFit) -> f64 {
    let tabulated = CLEARANCE_HOLE
        .iter()
        .find(|(d, _, _, _)| (d - major).abs() <= 1e-6);
    match (fit, tabulated) {
        (ClearanceFit::Custom(d), _) => d,
        (ClearanceFit::Close, Some(&(_, c, _, _))) => c,
        (ClearanceFit::Normal, Some(&(_, _, n, _))) => n,
        (ClearanceFit::Loose, Some(&(_, _, _, l))) => l,
        (ClearanceFit::Close, None) => 1.06 * major,
        (ClearanceFit::Normal, None) => 1.10 * major,
        (ClearanceFit::Loose, None) => 1.20 * major,
    }
}
