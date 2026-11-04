//! What a touching pair's contact is (PK-979, plan 4.12 item 1): a
//! point, a curve or an area, and for an area how much.
//!
//! The cell complex of the two world bodies decides area: the faces
//! their regions share are the contact patch, measured exactly, and a
//! patch of positive area is an area contact. Otherwise the kernel's
//! contact report tells a curve from a point. Its per-pair extent is
//! not the authority for area: at a patch corner it names the face
//! pair touching along the patch's edge (measured: a block offset on
//! another's face read as a curve). A touch the report finds no contact
//! for at the kernel tolerance (a gap under the interference tolerance)
//! has no extent. A kernel decline is reported with the pair, never
//! turned into an extent.

use popkern_bool::{SetExpr, TouchExtent, contact_report, partition};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::Serialize;

/// How far a touch extends.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ContactExtent {
    /// The bodies meet at isolated points.
    Point,
    /// They meet along a curve (an edge on a face, two edges).
    Curve,
    /// They share face area.
    Area,
}

/// A touching pair's contact.
#[derive(Debug, Clone, PartialEq, Default)]
pub(crate) struct ContactDetail {
    /// The widest extent the bodies touch over.
    pub extent: Option<ContactExtent>,
    /// The shared face area for an area contact.
    pub area: Option<f64>,
    /// Why the kernel could not say.
    pub decline: Option<String>,
}

/// The contact between world bodies `a` and `b`, at `tolerance`.
pub(crate) fn contact_detail(a: &BRep, b: &BRep, tolerance: f64) -> ContactDetail {
    let tol = Tolerance::default().with_linear(tolerance);
    let contacts = match contact_report(a, b, &tol) {
        Ok(c) => c,
        Err(e) => {
            return ContactDetail {
                decline: Some(format!("contact report: {e}")),
                ..ContactDetail::default()
            };
        }
    };
    let rank = |e: TouchExtent| match e {
        TouchExtent::Point => ContactExtent::Point,
        TouchExtent::Curve => ContactExtent::Curve,
        TouchExtent::Area => ContactExtent::Area,
    };
    let touched = contacts.iter().map(|c| rank(c.extent)).max();
    let Some(touched) = touched else {
        return ContactDetail::default();
    };
    let area = partition(&[a, b], tol)
        .and_then(|p| p.region_solids(&SetExpr::Any))
        .map_err(|e| e.to_string())
        .and_then(|regions| {
            let body = |i: usize| regions.solids.get(i).and_then(|r| r.bodies.iter().next());
            let mut total = 0.0;
            for interface in &regions.interfaces {
                let (Some(x), Some(y)) = (body(interface.solids.0), body(interface.solids.1))
                else {
                    continue;
                };
                if x == y {
                    continue;
                }
                for face in interface.patch.faces.keys() {
                    total += popkern_measure::face_area(&interface.patch, face, &tol)
                        .map_err(|e| e.to_string())?;
                }
            }
            Ok(total)
        });
    let below_area = |e: ContactExtent| e.min(ContactExtent::Curve);
    match area {
        Ok(area) if area > tolerance * tolerance => ContactDetail {
            extent: Some(ContactExtent::Area),
            area: Some(area),
            decline: None,
        },
        Ok(_) => ContactDetail {
            extent: Some(below_area(touched)),
            ..ContactDetail::default()
        },
        Err(e) => ContactDetail {
            extent: Some(touched),
            area: None,
            decline: Some(format!("contact area: {e}")),
        },
    }
}
