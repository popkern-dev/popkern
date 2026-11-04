//! Parts split by tools, imprinted, and partitioned into regions
//! (design: operations on the cell complex, 5.3).

use popkern_bool::{BodySet, BoolOptions, FaceClass, SetExpr};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::FaceHandle;

use crate::error::PartError;
use crate::part::Part;

/// One piece of a split part.
#[derive(Debug, Clone)]
pub struct SplitPart {
    /// The piece.
    pub part: Part,
    /// The tools (by index into the tools given) the piece lies inside.
    pub inside_tools: BodySet,
}

/// A part imprinted by a tool, each face classified against it.
#[derive(Debug, Clone)]
pub struct ImprintedPart {
    /// The part with its faces split where the tool crosses them.
    pub part: Part,
    classes: Vec<(FaceHandle, FaceClass)>,
}

impl ImprintedPart {
    /// Where `face` (of [`ImprintedPart::part`]) lies relative to the
    /// tool, when decided.
    pub fn class_of(&self, face: FaceHandle) -> Option<FaceClass> {
        self.classes
            .iter()
            .find(|(f, _)| *f == face)
            .map(|(_, c)| *c)
    }

    /// The part's faces of `class`.
    pub fn faces_classified(&self, class: FaceClass) -> impl Iterator<Item = FaceHandle> + '_ {
        self.classes
            .iter()
            .filter(move |(_, c)| *c == class)
            .map(|(f, _)| *f)
    }
}

/// One region of a partition of parts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartRegion {
    /// The parts (by index) holding it.
    pub parts: BodySet,
    /// A point strictly inside it.
    pub interior_point: Point3,
    /// Its volume.
    pub volume: f64,
}

/// The regions N parts cut space into.
pub struct PartPartition {
    parts: Vec<Part>,
    tolerance: Tolerance,
    partition: popkern_bool::Partition,
}

impl std::fmt::Debug for PartPartition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartPartition")
            .field("parts", &self.parts.len())
            .field("partition", &self.partition)
            .finish()
    }
}

fn boolean_error(operation: &'static str) -> impl Fn(popkern_bool::BoolError) -> PartError {
    move |source| PartError::Boolean { operation, source }
}

impl PartPartition {
    /// Every region, with its parts, an interior point and its volume
    /// (measured on the region as a solid).
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] when a region cannot be extracted;
    /// [`PartError::Measure`] when its volume cannot be measured.
    pub fn regions(&self) -> Result<Vec<PartRegion>, PartError> {
        let cells = self
            .partition
            .region_solids(&SetExpr::Any)
            .map_err(boolean_error("partition"))?;
        let infos = self.partition.regions();
        cells
            .solids
            .iter()
            .map(|s| {
                let volume =
                    popkern_measure::body_volume(&s.body, &self.tolerance).map_err(|source| {
                        PartError::Measure {
                            what: "region volume",
                            source,
                        }
                    })?;
                Ok(PartRegion {
                    parts: s.bodies,
                    interior_point: infos[s.region].interior_point,
                    volume,
                })
            })
            .collect()
    }

    /// The regions `expr` keeps, as one part.
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn select(&self, expr: &SetExpr, options: &BoolOptions) -> Result<Part, PartError> {
        let body = self
            .partition
            .select(expr, options)
            .map_err(boolean_error("partition select"))?;
        self.part_of(body, "partition select")
    }

    /// Every region `expr` keeps as its own part.
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn region_parts(&self, expr: &SetExpr) -> Result<Vec<Part>, PartError> {
        let cells = self
            .partition
            .region_solids(expr)
            .map_err(boolean_error("partition regions"))?;
        cells
            .solids
            .into_iter()
            .map(|s| self.part_of(s.body, "partition region"))
            .collect()
    }

    fn part_of(
        &self,
        body: popkern_core::topology::BRep,
        name: &'static str,
    ) -> Result<Part, PartError> {
        let refs: Vec<&Part> = self.parts.iter().collect();
        Part::derived(body, self.tolerance, name, &refs)
    }
}

impl Part {
    /// This part cut into its pieces by `tools`, every piece kept with
    /// the tools it lies inside (see [`popkern_bool::split`]).
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn split(&self, tools: &[&Part]) -> Result<Vec<SplitPart>, PartError> {
        let tolerance = tools.iter().fold(self.tolerance(), |t, p| {
            crate::part::finer(t, p.tolerance())
        });
        let breps: Vec<&popkern_core::topology::BRep> = tools.iter().map(|p| p.brep()).collect();
        let cut =
            popkern_bool::split(self.brep(), &breps, tolerance).map_err(boolean_error("split"))?;
        let mut inputs = vec![self];
        inputs.extend_from_slice(tools);
        cut.pieces
            .into_iter()
            .map(|p| {
                Ok(SplitPart {
                    part: Part::derived(p.body, tolerance, "split", &inputs)?,
                    inside_tools: p.inside_tools,
                })
            })
            .collect()
    }

    /// This part with its faces split where `tool` crosses them, each
    /// classified against the tool (see [`popkern_bool::ops::imprint`]).
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn imprint(&self, tool: &Part) -> Result<ImprintedPart, PartError> {
        let tolerance = crate::part::finer(self.tolerance(), tool.tolerance());
        let r = popkern_bool::ops::imprint::imprint(self.brep(), tool.brep(), &tolerance)
            .map_err(boolean_error("imprint"))?;
        Ok(ImprintedPart {
            part: Part::derived(r.imprinted_target, tolerance, "imprint", &[self, tool])?,
            classes: r.classes,
        })
    }

    /// The regions `parts` cut space into (see [`popkern_bool::partition`]).
    ///
    /// # Errors
    ///
    /// [`PartError::Boolean`] with the kernel diagnostics.
    pub fn partition(parts: &[&Part]) -> Result<PartPartition, PartError> {
        let tolerance = parts
            .iter()
            .map(|p| p.tolerance())
            .reduce(crate::part::finer)
            .unwrap_or_default();
        let breps: Vec<&popkern_core::topology::BRep> = parts.iter().map(|p| p.brep()).collect();
        let partition =
            popkern_bool::partition(&breps, tolerance).map_err(boolean_error("partition"))?;
        Ok(PartPartition {
            parts: parts.iter().map(|p| (*p).clone()).collect(),
            tolerance,
            partition,
        })
    }
}

#[cfg(test)]
#[path = "part_partition_tests.rs"]
mod tests;
