//! Mates (design `docs/architecture/assemblies.md` §8.2, PK-631):
//! geometric relations between entities of instances, solved to the
//! placements of the root instances that hold them.
//!
//! An entity is a face of the instance's part (by its stable name or a
//! label the part's script gave it) or one of the part's own origin
//! datums (`xy`/`yz`/`zx` planes, `x`/`y`/`z` axes). A face reduces to
//! its geometric primitive: a plane (point, outward normal), an axis
//! (cylinder, cone, torus: point, direction) or a point (sphere
//! centre). Free-form faces carry no primitive and are refused.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::names::InstancePath;
use crate::placement::Placement;

/// An origin datum of a part.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Datum {
    /// The XY plane (normal +z).
    Xy,
    /// The YZ plane (normal +x).
    Yz,
    /// The ZX plane (normal +y).
    Zx,
    /// The X axis.
    X,
    /// The Y axis.
    Y,
    /// The Z axis.
    Z,
}

impl Datum {
    /// Parse `xy`, `yz`, `zx` (or `xz`), `x`, `y`, `z`.
    pub fn parse(text: &str) -> Option<Datum> {
        Some(match text.to_ascii_lowercase().as_str() {
            "xy" => Datum::Xy,
            "yz" => Datum::Yz,
            "zx" | "xz" => Datum::Zx,
            "x" => Datum::X,
            "y" => Datum::Y,
            "z" => Datum::Z,
            _ => return None,
        })
    }

    fn local(self) -> Primitive {
        let o = Point3::new(0.0, 0.0, 0.0);
        match self {
            Datum::Xy => Primitive::Plane(o, Vec3::new(0.0, 0.0, 1.0)),
            Datum::Yz => Primitive::Plane(o, Vec3::new(1.0, 0.0, 0.0)),
            Datum::Zx => Primitive::Plane(o, Vec3::new(0.0, 1.0, 0.0)),
            Datum::X => Primitive::Axis(o, Vec3::new(1.0, 0.0, 0.0)),
            Datum::Y => Primitive::Axis(o, Vec3::new(0.0, 1.0, 0.0)),
            Datum::Z => Primitive::Axis(o, Vec3::new(0.0, 0.0, 1.0)),
        }
    }
}

/// What a mate refers to on one instance.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum MateEntity {
    /// A face of the part, by its spelled stable name.
    Face {
        /// The spelled entity name (`extrude-1/Lateral/0`).
        name: String,
    },
    /// An origin datum of the part.
    Datum {
        /// Which datum.
        datum: Datum,
    },
}

/// An entity of one occurrence.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MateRef {
    /// The occurrence (`planet_1`, `stage/planet_1`).
    pub instance: String,
    /// The entity on it.
    pub entity: MateEntity,
}

/// The relation a mate asks for.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum MateKind {
    /// Planes coplanar (normals opposed: faces touch; `aligned`: same
    /// facing), axes collinear, points coincident, or a point on a
    /// plane or an axis.
    Coincident {
        /// Normals point the same way instead of facing each other.
        #[serde(default)]
        aligned: bool,
    },
    /// Axes collinear (or a point on an axis).
    Concentric,
    /// Planes parallel at `distance` along the first plane's normal
    /// (or a point at `distance` from a plane).
    Distance {
        /// The signed offset.
        distance: f64,
        /// Normals point the same way instead of facing each other.
        #[serde(default)]
        aligned: bool,
    },
    /// The angle between two directions (plane normals or axes); signed
    /// about `about` when given, else unsigned.
    Angle {
        /// The angle in degrees.
        degrees: f64,
        /// The axis the angle is measured about (right-handed).
        #[serde(default)]
        about: Option<MateRef>,
    },
}

/// One mate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Mate {
    /// The relation.
    #[serde(flatten)]
    pub kind: MateKind,
    /// The first entity.
    pub a: MateRef,
    /// The second entity.
    pub b: MateRef,
}

impl Mate {
    /// "concentric planet_1:z with carrier_1:pin_0" — the mate in words.
    pub fn describe(&self) -> String {
        let kind = match &self.kind {
            MateKind::Coincident { aligned: false } => "coincident".to_string(),
            MateKind::Coincident { aligned: true } => "coincident (aligned)".to_string(),
            MateKind::Concentric => "concentric".to_string(),
            MateKind::Distance { distance, .. } => format!("distance {distance}"),
            MateKind::Angle { degrees, .. } => format!("angle {degrees}°"),
        };
        format!("{kind} {} with {}", spell(&self.a), spell(&self.b))
    }
}

fn spell(r: &MateRef) -> String {
    match &r.entity {
        MateEntity::Face { name } => format!("{}:{name}", r.instance),
        MateEntity::Datum { datum } => format!("{}:{datum:?}", r.instance).to_lowercase(),
    }
}

/// A resolved geometric primitive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Primitive {
    /// A point and an outward unit normal.
    Plane(Point3, Vec3),
    /// A point on the axis and a unit direction.
    Axis(Point3, Vec3),
    /// A point.
    Point(Point3),
}

impl Primitive {
    /// Move by `p`.
    pub fn placed(&self, p: &Placement) -> Primitive {
        let t = p.transform();
        match *self {
            Primitive::Plane(o, n) => Primitive::Plane(t.apply_point(&o), t.apply_direction(&n)),
            Primitive::Axis(o, d) => Primitive::Axis(t.apply_point(&o), t.apply_direction(&d)),
            Primitive::Point(o) => Primitive::Point(t.apply_point(&o)),
        }
    }

    /// The direction (plane normal or axis), when there is one.
    pub fn direction(&self) -> Option<Vec3> {
        match *self {
            Primitive::Plane(_, n) | Primitive::Axis(_, n) => Some(n),
            Primitive::Point(_) => None,
        }
    }
}

/// A mate entity resolved against the document: its primitive in the
/// root instance's frame and the root instance it rides on.
#[derive(Debug, Clone)]
pub struct Resolved {
    /// Index of the root instance carrying it.
    pub root: usize,
    /// The primitive in the root instance's local frame.
    pub local: Primitive,
}

impl Assembly {
    /// Resolve `r` to its root instance and its primitive in that
    /// instance's frame.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownInstance`]; [`AssemblyError::Mate`] for a
    /// face the part does not have or a free-form face.
    pub fn resolve_mate_ref(&mut self, r: &MateRef) -> Result<Resolved, AssemblyError> {
        let path = InstancePath::parse(&r.instance);
        let unknown = || AssemblyError::UnknownInstance {
            path: r.instance.clone(),
        };
        let root_name = path.0.first().cloned().ok_or_else(unknown)?;
        let root = self
            .root
            .iter()
            .position(|i| i.name == root_name)
            .ok_or_else(unknown)?;
        let rows = self.tree();
        let row = rows
            .iter()
            .find(|row| row.path == path)
            .ok_or_else(unknown)?;
        let root_world = self.root[root].placement;
        // The occurrence relative to its root instance (fixed while the
        // root moves).
        let inner = row.world.then(&root_world.inverse());
        let definition = row.definition;
        let local = match &r.entity {
            MateEntity::Datum { datum } => datum.local(),
            MateEntity::Face { name } => {
                let model = self
                    .part_model_mut(definition)
                    .map_err(|_| AssemblyError::Mate {
                        mate: r.instance.clone(),
                        detail: "a face mate needs a part occurrence, not a sub-assembly".into(),
                    })?;
                face_primitive(model, name).map_err(|detail| AssemblyError::Mate {
                    mate: format!("{}:{name}", r.instance),
                    detail,
                })?
            }
        };
        Ok(Resolved {
            root,
            local: local.placed(&inner),
        })
    }
}

/// The primitive of the face labelled `name` (a label the part's
/// script bound) or spelled `name` (a stable entity name), resolved by
/// the part's own selection rules so a face a later feature split or
/// renamed is still found; exactly one face must match.
fn face_primitive(model: &mut popkern_script::Model, name: &str) -> Result<Primitive, String> {
    let selector = if model.entity_by_label(name).is_ok() {
        popkern_script::selection::faces().labeled(name)
    } else {
        popkern_script::selection::faces().spelled(name)
    };
    let mut found = Vec::new();
    for body in model.alive_bodies() {
        let brep = model.brep_of(body).map_err(|e| e.to_string())?;
        let snapshot = model
            .resolve_faces_in(&brep, &selector)
            .map_err(|e| e.to_string())?;
        for entity in snapshot.names {
            for face in brep.faces.values() {
                if face.name.and_then(|n| brep.name_pool.get(n)) == Some(&entity) {
                    found.push(face.surface.clone());
                }
            }
        }
    }
    let surface = match found.as_slice() {
        [Some(surface)] => surface.clone(),
        [None] => return Err("the face has no surface".into()),
        [] => return Err(format!("the part has no face '{name}'")),
        more => {
            return Err(format!(
                "'{name}' names {} faces; a mate needs exactly one",
                more.len()
            ));
        }
    };
    let sign = if surface.orientation { 1.0 } else { -1.0 };
    match &surface.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => Ok(Primitive::Plane(
            *origin,
            normal.normalized().unwrap_or(*normal) * sign,
        )),
        SurfaceGeometry::Cylinder { origin, axis, .. }
        | SurfaceGeometry::Torus {
            center: origin,
            axis,
            ..
        } => Ok(Primitive::Axis(*origin, axis.normalized().unwrap_or(*axis))),
        SurfaceGeometry::Cone { apex, axis, .. } => {
            Ok(Primitive::Axis(*apex, axis.normalized().unwrap_or(*axis)))
        }
        SurfaceGeometry::Sphere { center, .. } => Ok(Primitive::Point(*center)),
        other => Err(format!(
            "a {} face has no plane, axis or centre to mate",
            other.kind_name()
        )),
    }
}
