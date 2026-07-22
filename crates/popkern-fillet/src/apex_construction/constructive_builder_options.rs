//! Execution options for constructive apex topology application.

use super::cap_loop_construction_options::CapLoopConstructionOptions;

/// How the V-free repeated-apex contract was requested.
///
/// V-free topology is the production DEFAULT for the N >= 4
/// repeated-apex partition class (PK-124): the Stage-D seam-ownership
/// contract holds by construction there and the pyramid4 zoo diagnostic
/// set is clean under it (0 slivers, 0 overlaps, 0 tessellation seams).
/// The N = 3 class is the production default since PK-153: the
/// CONSTRUCTION class closed under the strict env gate (PK-151
/// pyramid3 spoke anchoring, PK-129 pocket-junction settlement,
/// PK-128 mixed-winding replan; strict matrix 8/8), and the SERVING
/// class closed with the vestigial-shell consolidation (flap/sliver
/// retirement, shell winding reconciliation, seam welds, periodic
/// pcurve branch repair) plus the rim-trim short-sweep convention —
/// all 13 previously-broken workspace serving gates are green, with
/// the mixed cube serving a watertight 12-face shell at 0.25% area
/// parity and an empty fidelity-blocker census. PK-167 tracks the
/// constructive replacement for the winding repair (bridges as
/// in-place support-loop rewrites; Reauthored cycle primary).
/// `POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX=1` forces V-free
/// at every valence (diagnostic gates); `=0` restores the
/// apex-retaining path everywhere as a containment escape hatch.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum VFreeRepeatedApexMode {
    /// Promoted production default: V-free for N >= 3.
    #[default]
    DefaultOn,
    /// Environment-requested strict mode at every valence.
    ExplicitOn,
    /// Environment-requested apex-retaining path everywhere.
    Off,
}

impl VFreeRepeatedApexMode {
    pub(crate) fn from_environment_adapter() -> Self {
        match std::env::var("POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX") {
            Ok(value) if matches!(value.trim(), "0" | "false" | "off") => Self::Off,
            Ok(_) => Self::ExplicitOn,
            Err(_) => Self::DefaultOn,
        }
    }

    fn applies_to(self, n_fillets: usize) -> bool {
        match self {
            Self::ExplicitOn => true,
            Self::Off => false,
            Self::DefaultOn => n_fillets >= 3,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct ConstructiveApexExecutionOptions {
    cap_loop: CapLoopConstructionOptions,
    v_free_mode: VFreeRepeatedApexMode,
}

impl ConstructiveApexExecutionOptions {
    pub(crate) fn from_environment_adapter() -> Self {
        Self {
            cap_loop: CapLoopConstructionOptions::from_environment_adapter(),
            v_free_mode: VFreeRepeatedApexMode::from_environment_adapter(),
        }
    }

    /// Cap-loop options with the V-free flag resolved for this apex's
    /// fillet count.
    pub(super) fn cap_loop_for(self, n_fillets: usize) -> CapLoopConstructionOptions {
        self.cap_loop
            .with_v_free_repeated_apex(self.v_free_mode.applies_to(n_fillets))
    }

    /// Copy of these options with the V-free repeated-apex contract off
    /// (the apex-retaining path). The dispatch retries with this variant
    /// when the promoted V-free default fails native acceptance.
    pub(crate) fn with_v_free_repeated_apex_disabled(mut self) -> Self {
        self.v_free_mode = VFreeRepeatedApexMode::Off;
        self
    }

    pub(crate) fn require_v_free_repeated_apex_for(self, n_fillets: usize) -> bool {
        self.v_free_mode.applies_to(n_fillets)
    }
}
