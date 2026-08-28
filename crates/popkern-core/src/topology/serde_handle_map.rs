//! Serde representation for maps keyed by arena handles.
//!
//! `serde_json` cannot serialize maps whose keys are not strings, so
//! any `HashMap<HandleType, _>` field silently poisons whole-model
//! JSON serialization the first time an operator populates it
//! (measured: the vertex-blend pyramid apex serve stores 16
//! `surface_boundary_hints` entries and every `serde_json::to_vec`
//! of the body fails with "key must be a string"). These helpers
//! serialize such maps as a `Vec<(K, V)>` SORTED BY KEY — sorted so
//! that serialization is deterministic regardless of `HashMap`
//! iteration order, which bit-identical output fingerprints rely on.

use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::collections::HashMap;

/// Serialize a handle-keyed map as a key-sorted `Vec<(K, V)>`.
pub fn serialize<K, V, S>(map: &HashMap<K, V>, serializer: S) -> Result<S::Ok, S::Error>
where
    K: Ord + Copy + Serialize + std::hash::Hash,
    V: Serialize,
    S: Serializer,
{
    let mut pairs: Vec<(&K, &V)> = map.iter().collect();
    pairs.sort_by_key(|(k, _)| **k);
    let mut seq = serializer.serialize_seq(Some(pairs.len()))?;
    for (k, v) in pairs {
        seq.serialize_element(&(k, v))?;
    }
    seq.end()
}

/// Deserialize a `Vec<(K, V)>` back into a handle-keyed map.
pub fn deserialize<'de, K, V, D>(deserializer: D) -> Result<HashMap<K, V>, D::Error>
where
    K: Eq + std::hash::Hash + Deserialize<'de>,
    V: Deserialize<'de>,
    D: Deserializer<'de>,
{
    let pairs: Vec<(K, V)> = Vec::deserialize(deserializer)?;
    Ok(pairs.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// A body carrying handle-keyed auxiliary maps must round-trip
    /// through JSON (regression: "key must be a string").
    #[test]
    fn brep_with_handle_keyed_maps_roundtrips_json() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let he = brep.half_edges.keys().next().unwrap();
        brep.surface_boundary_hints
            .insert(he, crate::geometry::edge_triple::BoundarySide::VMin);
        brep.nonmanifold_radial.insert(he, he);
        let bytes = serde_json::to_vec(&brep).expect("hint-carrying body must serialize");
        let back: BRep = serde_json::from_slice(&bytes).expect("roundtrip");
        assert_eq!(back.surface_boundary_hints.len(), 1);
        assert_eq!(back.nonmanifold_radial.len(), 1);
        assert_eq!(
            back.surface_boundary_hints.get(&he),
            brep.surface_boundary_hints.get(&he)
        );
    }

    /// Serialization must be byte-deterministic across identical
    /// maps built in different insertion orders (fingerprint
    /// contract).
    #[test]
    fn handle_map_serialization_is_order_independent() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let handles: Vec<_> = brep.half_edges.keys().take(4).collect();
        for &h in &handles {
            brep.surface_boundary_hints
                .insert(h, crate::geometry::edge_triple::BoundarySide::VMin);
        }
        let forward = serde_json::to_vec(&brep).unwrap();
        let mut reversed = brep.clone();
        reversed.surface_boundary_hints.clear();
        for &h in handles.iter().rev() {
            reversed
                .surface_boundary_hints
                .insert(h, crate::geometry::edge_triple::BoundarySide::VMin);
        }
        let backward = serde_json::to_vec(&reversed).unwrap();
        assert_eq!(
            forward, backward,
            "sorted-pair serialization must not depend on insertion order"
        );
    }
}
