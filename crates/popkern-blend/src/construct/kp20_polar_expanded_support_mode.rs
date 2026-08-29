//! Expanded KP20 support-lane reauthoring modes.

#[derive(Debug, Clone, Copy)]
pub(in crate::construct) enum Kp20ExpandedSupportMode {
    SupportNormalRelief,
    ChordRelief,
    InterfaceContract,
    CoupledLaneSpread,
}

impl Kp20ExpandedSupportMode {
    pub(in crate::construct) fn all() -> [Self; 4] {
        [
            Self::SupportNormalRelief,
            Self::ChordRelief,
            Self::InterfaceContract,
            Self::CoupledLaneSpread,
        ]
    }

    pub(in crate::construct) fn diagnostic_label(self) -> &'static str {
        match self {
            Self::SupportNormalRelief => "support_normal_relief",
            Self::ChordRelief => "chord_relief",
            Self::InterfaceContract => "interface_contract",
            Self::CoupledLaneSpread => "coupled_lane_spread",
        }
    }
}
