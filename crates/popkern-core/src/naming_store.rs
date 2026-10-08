//! Copy-on-write storage for [`crate::naming::NamePool`] (PK-731).
//!
//! A pool only grows during modelling: every operation appends its
//! names, and a body's pool carries every name of its history. Each
//! BRep clone used to deep-copy all of it — every name, its ancestor
//! list, its lineage values — so a feature recorded late in a long
//! history cost in proportion to the whole history (a 64-slot disc:
//! the 64th slot cut 20× slower than the first). [`ChunkedVec`] keeps
//! the entries in fixed-size chunks shared ([`Arc`]) between clones: a
//! clone copies one pointer per chunk, and a write copies only the
//! chunk it lands in. Shared chunks also make "this pool extends that
//! one" a pointer comparison ([`ChunkedVec::extends`]).

use std::sync::Arc;

use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Entries per chunk.
const CHUNK: usize = 256;

/// A vector stored as shared fixed-size chunks (see the module docs).
#[derive(Debug, Clone)]
pub(crate) struct ChunkedVec<T> {
    chunks: Vec<Arc<Vec<T>>>,
    len: usize,
}

impl<T> Default for ChunkedVec<T> {
    fn default() -> Self {
        ChunkedVec {
            chunks: Vec::new(),
            len: 0,
        }
    }
}

impl<T: Clone> ChunkedVec<T> {
    /// The number of entries.
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// Entry `i`.
    pub(crate) fn get(&self, i: usize) -> Option<&T> {
        if i >= self.len {
            return None;
        }
        self.chunks.get(i / CHUNK)?.get(i % CHUNK)
    }

    /// Entry `i`, mutable: copies its chunk first when shared.
    pub(crate) fn get_mut(&mut self, i: usize) -> Option<&mut T> {
        if i >= self.len {
            return None;
        }
        Arc::make_mut(self.chunks.get_mut(i / CHUNK)?).get_mut(i % CHUNK)
    }

    /// Append `value`.
    pub(crate) fn push(&mut self, value: T) {
        match self.chunks.last_mut() {
            Some(last) if last.len() < CHUNK => Arc::make_mut(last).push(value),
            _ => {
                let mut chunk = Vec::with_capacity(CHUNK);
                chunk.push(value);
                self.chunks.push(Arc::new(chunk));
            }
        }
        self.len += 1;
    }

    /// The entries in order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = &T> {
        self.chunks.iter().flat_map(|chunk| chunk.iter())
    }

    /// Visit every entry mutably. Chunks are copied only where `edit`
    /// reports a change it wants to make (`wants` returns true).
    pub(crate) fn edit_where(
        &mut self,
        mut wants: impl FnMut(usize, &T) -> bool,
        mut edit: impl FnMut(usize, &mut T),
    ) {
        for (c, chunk) in self.chunks.iter_mut().enumerate() {
            let base = c * CHUNK;
            let hits: Vec<usize> = chunk
                .iter()
                .enumerate()
                .filter(|(k, v)| wants(base + k, v))
                .map(|(k, _)| k)
                .collect();
            if hits.is_empty() {
                continue;
            }
            let owned = Arc::make_mut(chunk);
            for k in hits {
                if let Some(v) = owned.get_mut(k) {
                    edit(base + k, v);
                }
            }
        }
    }

    /// Build from owned entries.
    pub(crate) fn from_vec(values: Vec<T>) -> Self {
        let mut out = ChunkedVec::default();
        for value in values {
            out.push(value);
        }
        out
    }
}

impl<T: Clone + PartialEq> ChunkedVec<T> {
    /// Whether `self` begins with every entry of `prefix`: whole chunks
    /// compare by pointer (shared, so equal), falling back to values.
    pub(crate) fn extends(&self, prefix: &ChunkedVec<T>) -> bool {
        if prefix.len > self.len {
            return false;
        }
        prefix
            .chunks
            .iter()
            .zip(&self.chunks)
            .all(|(p, s)| Arc::ptr_eq(p, s) || s.get(..p.len()) == Some(p.as_slice()))
    }
}

impl<T: Serialize> Serialize for ChunkedVec<T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.collect_seq(self.chunks.iter().flat_map(|chunk| chunk.iter()))
    }
}

impl<'de, T: Deserialize<'de> + Clone> Deserialize<'de> for ChunkedVec<T> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        Vec::<T>::deserialize(deserializer).map(ChunkedVec::from_vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clones_share_chunks_until_written() {
        let mut a = ChunkedVec::from_vec((0..600).collect::<Vec<u32>>());
        let mut b = a.clone();
        assert!(b.extends(&a) && a.extends(&b));
        b.push(600);
        assert!(b.extends(&a), "an append keeps the prefix");
        assert!(!a.extends(&b));
        if let Some(v) = b.get_mut(3) {
            *v = 99;
        }
        assert!(!b.extends(&a), "a write breaks it");
        assert_eq!(
            a.get(3),
            Some(&3),
            "the clone's write never reaches the original"
        );
        a.push(7);
        assert_eq!(a.get(600), Some(&7));
        assert_eq!(a.iter().count(), 601);
    }

    #[test]
    fn serializes_as_a_flat_sequence() {
        let v = ChunkedVec::from_vec((0..300).collect::<Vec<u32>>());
        let json = serde_json::to_string(&v).unwrap_or_default();
        let flat = serde_json::to_string(&(0..300).collect::<Vec<u32>>()).unwrap_or_default();
        assert_eq!(json, flat);
        let back: ChunkedVec<u32> = serde_json::from_str(&json).unwrap_or_default();
        assert_eq!(back.len(), 300);
        assert!(back.extends(&v));
    }
}
