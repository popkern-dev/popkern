//! Thermal contacts in the steady system (`crate::contact`): a contact
//! conductance h_c on the temperature jump across each shared patch.

use crate::error::AnalysisError;

/// The reduced contact couplings of `prepared` (`crate::contact`):
/// conductance h_c on the temperature jump, W/(mm²·K).
pub(super) fn contact_triplets(
    prepared: &super::system::Prepared,
    offsets: &[usize],
) -> Result<Vec<(usize, usize, f64)>, AnalysisError> {
    let mut out = Vec::new();
    for contact in &prepared.contacts {
        let (a, b) = (
            contact.patch.solids.0.0 as usize,
            contact.patch.solids.1.0 as usize,
        );
        let solid = |s: usize| {
            prepared
                .solids
                .get(s)
                .zip(offsets.get(s))
                .ok_or(AnalysisError::Unsupported {
                    what: "a contact names a solid the model does not have",
                })
        };
        let ((_, ma, ga), oa) = solid(a)?;
        let ((_, mb, gb), ob) = solid(b)?;
        let ea = ga
            .extension
            .to_csr()
            .map_err(|e| super::internal("contact", e))?;
        let eb = gb
            .extension
            .to_csr()
            .map_err(|e| super::internal("contact", e))?;
        let sa = crate::contact::Side {
            moments: ma,
            extension: &ga.extension,
            reduce: &ea,
            offset: *oa,
        };
        let sb = crate::contact::Side {
            moments: mb,
            extension: &gb.extension,
            reduce: &eb,
            offset: *ob,
        };
        let coefficient = 1e-6 * contact.conductance;
        let tables = &prepared.tables;
        out.extend(crate::contact::coupling_triplets(
            &prepared.grid,
            &contact.patch,
            (&sa, &sb),
            1,
            |s, _| {
                tables
                    .facet_mass(s)
                    .iter()
                    .map(|row| row.iter().map(|v| coefficient * v).collect())
                    .collect()
            },
        )?);
    }
    Ok(out)
}
