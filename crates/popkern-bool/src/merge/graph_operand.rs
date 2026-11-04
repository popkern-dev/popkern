//! Which operand a piece of graph topology comes from.

/// Identifies which boolean operand a piece of topology comes from: its
/// index among the operands (`A` = 0, `B` = 1 for a two-operand
/// boolean). An index, not two cases, so the graph can hold N operands
/// (design: operations on the cell complex, 7.1).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Operand(u8);

/// `A` and `B` as the two-operand logs and artifacts have always read;
/// further operands by index.
impl std::fmt::Debug for Operand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.0 {
            0 => f.write_str("A"),
            1 => f.write_str("B"),
            n => write!(f, "Operand({n})"),
        }
    }
}

impl Operand {
    /// The first (left-hand) operand.
    pub(crate) const A: Operand = Operand(0);
    /// The second (right-hand) operand.
    pub(crate) const B: Operand = Operand(1);

    /// The operand at `index` among the operands.
    pub(crate) fn from_index(index: usize) -> Self {
        Operand(u8::try_from(index).unwrap_or(u8::MAX))
    }

    /// Its index among the operands.
    pub(crate) fn index(self) -> usize {
        usize::from(self.0)
    }
}
