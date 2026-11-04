//! A set of operand indices: which operands contain a region.
//!
//! General fuse may take more than 64 operands, so the set is a
//! growable bitset rather than a `u64` (the cellular `Cell` mask).

/// Set of operand indices (bit `k` set ⇔ operand `k` is in the set).
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct OperandSet {
    words: Vec<u64>,
}

impl OperandSet {
    /// The empty set.
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Add operand `k`.
    pub(crate) fn insert(&mut self, k: usize) {
        let (w, b) = (k / 64, k % 64);
        if self.words.len() <= w {
            self.words.resize(w + 1, 0);
        }
        self.words[w] |= 1u64 << b;
    }

    /// Whether operand `k` is in the set.
    pub(crate) fn contains(&self, k: usize) -> bool {
        let (w, b) = (k / 64, k % 64);
        self.words
            .get(w)
            .is_some_and(|word| word & (1u64 << b) != 0)
    }

    /// Number of operands in the set.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    /// Whether the set is empty.
    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }

    /// Operand indices in increasing order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.words.iter().enumerate().flat_map(|(wi, &word)| {
            (0..64)
                .filter(move |b| word & (1u64 << b) != 0)
                .map(move |b| wi * 64 + b)
        })
    }
}
