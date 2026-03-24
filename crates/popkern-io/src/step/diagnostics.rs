//! Structured diagnostic log for import/export operations.

/// Severity level for a diagnostic entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// Non-fatal issue — import continues.
    Warning,
    /// Fatal issue for the affected entity — entity is skipped.
    Error,
}

/// Phase of the import/export pipeline where the diagnostic occurred.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// Part 21 parsing.
    Parse,
    /// Entity reference resolution.
    Resolve,
    /// Geometry mapping.
    MapGeometry,
    /// Topology mapping.
    MapTopology,
    /// Product structure mapping.
    MapProduct,
    /// Presentation (color, layer) mapping.
    MapPresentation,
    /// BRep validation.
    Validate,
    /// Healing passes.
    Heal,
}

/// A single diagnostic entry.
#[derive(Debug, Clone)]
pub struct DiagnosticEntry {
    /// Severity level.
    pub severity: Severity,
    /// Pipeline phase.
    pub phase: Phase,
    /// The STEP entity ID involved, if any.
    pub entity_id: Option<u64>,
    /// Human-readable description.
    pub message: String,
}

/// Collects diagnostic entries during import/export.
#[derive(Debug, Clone, Default)]
pub struct DiagnosticLog {
    entries: Vec<DiagnosticEntry>,
}

impl DiagnosticLog {
    /// Create an empty log.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a diagnostic entry.
    pub fn push(&mut self, entry: DiagnosticEntry) {
        self.entries.push(entry);
    }

    /// Log a warning.
    pub fn warn(&mut self, phase: Phase, entity_id: Option<u64>, message: impl Into<String>) {
        self.push(DiagnosticEntry {
            severity: Severity::Warning,
            phase,
            entity_id,
            message: message.into(),
        });
    }

    /// Log an error.
    pub fn error(&mut self, phase: Phase, entity_id: Option<u64>, message: impl Into<String>) {
        self.push(DiagnosticEntry {
            severity: Severity::Error,
            phase,
            entity_id,
            message: message.into(),
        });
    }

    /// Iterate over warnings.
    pub fn warnings(&self) -> impl Iterator<Item = &DiagnosticEntry> {
        self.entries
            .iter()
            .filter(|e| e.severity == Severity::Warning)
    }

    /// Iterate over errors.
    pub fn errors(&self) -> impl Iterator<Item = &DiagnosticEntry> {
        self.entries
            .iter()
            .filter(|e| e.severity == Severity::Error)
    }

    /// True if any error-level entries exist.
    pub fn has_errors(&self) -> bool {
        self.entries.iter().any(|e| e.severity == Severity::Error)
    }

    /// Total number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True if no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// All entries as a slice.
    pub fn entries(&self) -> &[DiagnosticEntry] {
        &self.entries
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_log_has_no_errors() {
        let log = DiagnosticLog::new();
        assert!(!log.has_errors());
        assert!(log.is_empty());
    }

    #[test]
    fn warn_and_error() {
        let mut log = DiagnosticLog::new();
        log.warn(Phase::Parse, None, "minor issue");
        log.error(Phase::Parse, Some(42), "bad entity");
        assert!(log.has_errors());
        assert_eq!(log.len(), 2);
        assert_eq!(log.warnings().count(), 1);
        assert_eq!(log.errors().count(), 1);
    }
}
