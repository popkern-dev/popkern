//! ISO size-table lookups for a [`ThreadSpec`] in its own length unit
//! (PK-462). The tables and their derived proportions are millimetres:
//! a spec whose lengths are SI metres (a mm-unit document) is looked
//! up at `major / millimetre` and its dimensions scaled back, so an M6
//! gets the ISO 4017 head (s = 10 mm) at either scale instead of
//! missing every table row and falling back to derived proportions.

use crate::spec::ThreadSpec;
use crate::standards::{self, ClearanceFit, HexHeadDims, SocketHeadDims};

impl ThreadSpec {
    /// The major diameter in millimetres (the table key).
    fn major_mm(&self) -> f64 {
        self.major_diameter() / self.millimetre()
    }

    /// ISO 4017 hex head (or its derived proportions), in this
    /// spec's unit.
    pub fn iso_hex_head(&self) -> HexHeadDims {
        let mm = standards::hex_head(self.major_mm());
        HexHeadDims {
            across_flats: mm.across_flats * self.millimetre(),
            height: mm.height * self.millimetre(),
        }
    }

    /// ISO 4762 socket head (or its derived proportions), in this
    /// spec's unit.
    pub fn iso_socket_head(&self) -> SocketHeadDims {
        let mm = standards::socket_head(self.major_mm());
        SocketHeadDims {
            diameter: mm.diameter * self.millimetre(),
            height: mm.height * self.millimetre(),
            socket_across_flats: mm.socket_across_flats * self.millimetre(),
        }
    }

    /// ISO 10642 countersunk head diameter, in this spec's unit.
    pub fn iso_countersunk_head_diameter(&self) -> f64 {
        standards::countersunk_head_diameter(self.major_mm()) * self.millimetre()
    }

    /// ISO 4032 hex nut thickness, in this spec's unit.
    pub fn iso_hex_nut_thickness(&self) -> f64 {
        standards::hex_nut_thickness(self.major_mm()) * self.millimetre()
    }

    /// ISO 273 clearance hole diameter, in this spec's unit. A custom
    /// fit's diameter is already in this unit and passes through.
    pub fn iso_clearance_diameter(&self, fit: ClearanceFit) -> f64 {
        match fit {
            ClearanceFit::Custom(d) => d,
            fit => standards::clearance_diameter(self.major_mm(), fit) * self.millimetre(),
        }
    }
}

#[cfg(test)]
#[path = "spec_standards_tests.rs"]
mod tests;
