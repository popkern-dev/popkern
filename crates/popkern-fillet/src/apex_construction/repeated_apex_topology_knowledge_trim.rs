//! Compatibility shim for repeated-apex trim endpoint knowledge.

pub(super) use popkern_core::topology::{
    BrepTrimEndpointAgreementFacts as TrimEndpointAgreementFacts,
    collect_brep_trim_endpoint_agreements as collect_trim_endpoint_agreements,
};
