//! A pose of several joints in one request (PK-853): a film frame moves
//! the six-axis arm's six joints and its grip together. One mate solve
//! and one placement envelope serve the frame; seven `set_param` calls
//! solved seven times and sent seven 270 kB envelopes (63 ms native,
//! about 600 ms in wasm per frame).

use popkern_assembly::Document;
use popkern_rhai::from_assembly_error;

use crate::envelope::Envelope;
use crate::session::Session;

impl Session {
    /// Set every joint in `values` (`name`, value in its unit) and
    /// re-place once. Joints are placement-only assembly parameters
    /// (they drive mates, not parts); when any name is not one, each is
    /// set in turn through [`Session::set_param`] and the last envelope
    /// (or the first failure) is the reply.
    pub fn set_joints(&mut self, values: &[(String, f64)]) -> Envelope {
        if values.is_empty() || !values.iter().all(|(n, _)| self.moves_placements_only(n)) {
            let mut last = None;
            for (name, value) in values {
                let envelope = self.set_param(name, *value);
                let failed = !envelope.outcome.is_served();
                last = Some(envelope);
                if failed {
                    break;
                }
            }
            return last.unwrap_or_else(|| self.describe_placements());
        }
        let Some(Document::Assembly(assembly)) = self.document.as_mut() else {
            return self.describe_placements();
        };
        let pose: Vec<(&str, f64)> = values.iter().map(|(n, v)| (n.as_str(), *v)).collect();
        // `set_param`'s order: the drives' solve, then the mates' solve.
        let result = assembly.set_mate_params(&pose).and_then(|()| {
            if !assembly.mates().is_empty() {
                assembly.solve_mates()?;
            }
            Ok(())
        });
        match result {
            Ok(()) => self.describe_placements(),
            Err(error) => {
                let names: Vec<&str> = pose.iter().map(|(n, _)| *n).collect();
                let outcome =
                    from_assembly_error(&error, &format!("set_joints {}", names.join(", ")));
                self.describe(outcome, Vec::new(), None)
            }
        }
    }
}
