//! A fast, deterministic hasher for the kernel's internal maps
//! (PK-1223).
//!
//! The standard library's maps hash with SipHash-1-3 under a random
//! per-process key: collision-resistant against adversarial input,
//! which the kernel's keys (handles, index pairs, quantized points)
//! never are, and slow on them: hashing held 25-30% of the mesher's
//! angular refinement and 5-10% of boolean and blend profiles. This is
//! the Fx hash (the multiply-rotate word hash the Rust compiler uses for
//! its own tables): a few instructions per word, and deterministic, so a
//! map's iteration order no longer varies from run to run.

use std::hash::{BuildHasherDefault, Hasher};

/// The Fx multiply constant (the golden ratio scaled to 64 bits).
const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

/// The Fx word hasher.
#[derive(Debug, Default, Clone, Copy)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    #[inline]
    fn add(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(5) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let (words, rest) = bytes.as_chunks::<8>();
        for word in words {
            self.add(u64::from_le_bytes(*word));
        }
        if !rest.is_empty() {
            let mut word = [0u8; 8];
            word[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(word));
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(u64::from(i));
    }

    #[inline]
    fn write_u16(&mut self, i: u16) {
        self.add(u64::from(i));
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(u64::from(i));
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

/// The builder of [`FxHasher`]s.
pub type FxBuildHasher = BuildHasherDefault<FxHasher>;

/// A hash map on the Fx hasher: build with `default()` or
/// `with_capacity_and_hasher(n, Default::default())`.
pub type FastHashMap<K, V> = std::collections::HashMap<K, V, FxBuildHasher>;

/// A hash set on the Fx hasher.
pub type FastHashSet<T> = std::collections::HashSet<T, FxBuildHasher>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::hash::BuildHasher;

    #[test]
    fn equal_keys_hash_equal_and_distinct_pairs_spread() {
        let b = FxBuildHasher::default();
        let h = |k: &(usize, usize)| b.hash_one(k);
        assert_eq!(h(&(3, 7)), h(&(3, 7)));
        let mut seen = FastHashSet::default();
        for a in 0..64 {
            for c in 0..64 {
                seen.insert(h(&(a, c)));
            }
        }
        assert_eq!(seen.len(), 64 * 64, "index pairs collide");
    }

    #[test]
    fn a_map_iterates_in_one_order_every_time() {
        let fill = || {
            let mut m: FastHashMap<u64, u64> = FastHashMap::default();
            for k in 0..1000u64 {
                m.insert(k.wrapping_mul(2654435761), k);
            }
            m.keys().copied().collect::<Vec<_>>()
        };
        assert_eq!(fill(), fill());
    }
}
