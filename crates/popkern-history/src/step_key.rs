//! Content addresses of timeline steps (the step cache's keys).
//!
//! A step's key chains the previous step's key with everything the
//! step's execution reads that is not already part of an earlier
//! step: the feature definition and identity, its body and operation
//! id, the label its report row carries, the evaluated value, unit and
//! expression of every parameter it uses, its resolved suppression,
//! whether a pattern in the replayed range seeds from it (seeds do not
//! auto-merge),
//! its pattern instance operation ids, the build tolerance and the
//! continuity-certification policy. Two replays whose keys agree up to
//! step *i* executed identical inputs through step *i*, so a state
//! cached under the key is the state a fresh replay would produce —
//! correct by construction, with no invalidation bookkeeping: any edit
//! (a parameter value, a configuration switch, an inserted or
//! reordered feature) changes the key of the first step it affects and
//! of every step after it.

use popkern_core::tolerance::Tolerance;
use slotmap::Key;

use crate::FeatureTree;
use crate::config_eval;
use crate::feature_params::feature_param_ids;
use crate::measure_query::MeasureContext;

/// The content address of one timeline step (128-bit FNV-1a chain).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct StepKey(u128);

const FNV_OFFSET: u128 = 0x6c62_272e_07bb_0142_62b8_2175_6295_c58d;
const FNV_PRIME: u128 = 0x0000_0000_0100_0000_0000_0000_0000_013b;

struct Fnv(u128);

impl Fnv {
    fn bytes(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 ^= u128::from(b);
            self.0 = self.0.wrapping_mul(FNV_PRIME);
        }
        // Length-delimit every field so adjacent fields cannot alias.
        for b in (bytes.len() as u64).to_le_bytes() {
            self.0 ^= u128::from(b);
            self.0 = self.0.wrapping_mul(FNV_PRIME);
        }
    }

    fn text(&mut self, text: &str) {
        self.bytes(text.as_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }
}

impl FeatureTree {
    /// The keys of timeline steps `0..end` under the current parameter
    /// values (call after evaluating parameters, as replay does). A
    /// step whose definition cannot be serialized gets no key, and
    /// neither does any step after it: those steps are never cached.
    pub(crate) fn step_keys(&self, end: usize, tolerance: &Tolerance) -> Vec<Option<StepKey>> {
        let mut keys = Vec::with_capacity(end);
        let mut chain = Some({
            let mut h = Fnv(FNV_OFFSET);
            h.text(&format!("{tolerance:?}"));
            h.text(&format!("{:?}", self.deferred.policy));
            h.0
        });
        let seeds = self.pattern_seed_ids(end);
        for &id in &self.timeline[..end.min(self.timeline.len())] {
            chain = chain.and_then(|prev| self.step_hash(prev, id, seeds.contains(&id)));
            keys.push(chain.map(StepKey));
        }
        keys
    }

    fn step_hash(&self, prev: u128, id: crate::node::FeatureId, seed: bool) -> Option<u128> {
        let node = self.features.get(id)?;
        let mut h = Fnv(prev);
        h.u64(id.data().as_ffi());
        h.u64(u64::from(seed));
        h.bytes(&serde_json::to_vec(&node.feature).ok()?);
        h.u64(u64::from(node.body.0));
        h.u64(node.operation_id.0);
        h.text(&node.label);
        let config = config_eval::resolve_suppression(id, self.active_config, &self.configs);
        h.text(&format!("{config:?}/{}", node.suppressed));
        if let Some(expr) = &node.suppress_if {
            let value =
                crate::expression::evaluate(expr, &self.parameter_table, &MeasureContext::empty());
            h.text(&format!("{expr:?}={value:?}"));
        }
        for (role, pid) in feature_param_ids(&node.feature) {
            h.text(&role);
            match self.parameter_table.get(pid) {
                Some(p) => {
                    h.u64(p.value.to_bits());
                    h.text(&format!("{:?}/{:?}", p.unit, p.expression));
                }
                None => h.text("missing"),
            }
        }
        let mut instances: Vec<_> = node
            .instance_op_ids
            .iter()
            .map(|(slot, op)| (*slot, op.0))
            .collect();
        instances.sort();
        h.text(&format!("{instances:?}"));
        Some(h.0)
    }
}
