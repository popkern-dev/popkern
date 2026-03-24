//! Frame computation for orienting profiles along sweep paths.

pub mod auxiliary;
pub mod bishop;
pub mod fixed;
pub mod frenet;

use crate::arc_length::ArcLengthTable;
use crate::error::{SweepError, SweepWarning};
use crate::types::FrameMethod;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

/// An orthonormal frame at a point along a path.
#[derive(Debug, Clone)]
pub struct PathFrame {
    /// Position on the curve.
    pub position: Point3,
    /// Unit tangent vector.
    pub tangent: Vec3,
    /// Unit normal vector.
    pub normal: Vec3,
    /// Unit binormal vector (tangent x normal).
    pub binormal: Vec3,
}

/// Sequence of frames along a path.
#[derive(Debug)]
pub struct FrameSequence {
    /// Ordered frames along the path.
    pub frames: Vec<PathFrame>,
    /// Arc-length values corresponding to each frame.
    pub s_values: Vec<f64>,
}

/// Compute frames along a path using the specified method.
pub fn compute_frames(
    curve: &NurbsCurve,
    table: &ArcLengthTable,
    method: &FrameMethod,
    n_samples: usize,
    twist: f64,
    is_closed: bool,
) -> Result<(FrameSequence, Vec<SweepWarning>), SweepError> {
    match method {
        FrameMethod::Bishop => {
            let seq = bishop::bishop_frames(curve, table, n_samples, twist, is_closed)?;
            Ok((seq, vec![]))
        }
        FrameMethod::Frenet => frenet::frenet_frames(curve, table, n_samples, twist),
        FrameMethod::Fixed(dir) => {
            let seq = fixed::fixed_frames(curve, table, *dir, n_samples, twist)?;
            Ok((seq, vec![]))
        }
        FrameMethod::Auxiliary(aux) => {
            let seq = auxiliary::auxiliary_frames(curve, aux, table, n_samples, twist)?;
            Ok((seq, vec![]))
        }
    }
}
