//! Core scene state for the debug viewer.
//!
//! [`VizServer`] holds the loaded BRep model, its tessellated mesh,
//! camera, selection state, and display options.

use glam::Vec3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use popkern_mesh::TriMesh;

use crate::camera::Camera;
use crate::error::VizError;

/// Selection state: what entity (if any) is currently highlighted.
#[derive(Debug, Clone)]
pub enum Selection {
    /// Nothing is selected.
    None,
    /// A BRep face is selected.
    Face(FaceHandle),
    /// A BRep edge is selected.
    Edge(EdgeHandle),
    /// A BRep vertex is selected.
    Vertex(VertexHandle),
}

/// How faces are colored in the viewport.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorMode {
    /// Color each face by its underlying surface type (plane, cylinder, etc.).
    Surface,
    /// Flat single gray for all faces.
    Flat,
    /// No face coloring — wireframe only.
    None,
}

/// Display toggles controlling how the model is rendered.
#[derive(Debug, Clone)]
pub struct DisplayOptions {
    /// Whether to draw a wireframe overlay on top of shaded faces.
    pub show_wireframe: bool,
    /// Whether to draw surface normals as debug lines.
    pub show_normals: bool,
    /// Face coloring strategy.
    pub color_mode: ColorMode,
}

impl Default for DisplayOptions {
    fn default() -> Self {
        Self {
            show_wireframe: true,
            show_normals: false,
            color_mode: ColorMode::Surface,
        }
    }
}

/// Core scene state for the debug viewer.
///
/// Holds the loaded BRep model and its tessellated mesh, the camera,
/// the current selection, and display options. Created with
/// [`VizServer::new`] and populated via [`VizServer::load_brep`].
pub struct VizServer {
    /// The loaded BRep model, if any.
    pub brep: Option<BRep>,
    /// Tessellated triangle mesh derived from the BRep.
    pub mesh: Option<TriMesh>,
    /// Camera controlling the viewport.
    pub camera: Camera,
    /// Current selection state.
    pub selection: Selection,
    /// Display toggles.
    pub display: DisplayOptions,
}

impl VizServer {
    /// Create a new scene with no model loaded and a default isometric camera.
    pub fn new() -> Self {
        Self {
            brep: Option::None,
            mesh: Option::None,
            camera: Camera::default(),
            selection: Selection::None,
            display: DisplayOptions::default(),
        }
    }

    /// Load a BRep model, tessellating it into a triangle mesh.
    ///
    /// Replaces any previously loaded model and resets the selection.
    /// The camera is fitted to the bounding box of the resulting mesh.
    pub fn load_brep(&mut self, brep: BRep) -> Result<(), VizError> {
        let options = popkern_mesh::TessellationOptions::display();
        let mesh =
            popkern_mesh::tessellate(&brep, &options).map_err(|e| VizError::Load(e.to_string()))?;

        // Fit camera to model bounds
        if let Some((min, max)) = Self::mesh_bounds(&mesh) {
            self.camera.fit_to_bounds(min, max);
        }

        self.brep = Some(brep);
        self.mesh = Some(mesh);
        self.selection = Selection::None;
        Ok(())
    }

    /// Compute the axis-aligned bounding box of the loaded mesh.
    ///
    /// Returns `None` if no mesh is loaded or the mesh has no vertices.
    pub fn model_bounds(&self) -> Option<(Vec3, Vec3)> {
        self.mesh.as_ref().and_then(Self::mesh_bounds)
    }

    /// Update the current selection.
    pub fn select(&mut self, sel: Selection) {
        self.selection = sel;
    }

    /// Compute the AABB of a [`TriMesh`], returning `None` if empty.
    fn mesh_bounds(mesh: &TriMesh) -> Option<(Vec3, Vec3)> {
        let first = mesh.vertices.first()?;
        let mut min = Vec3::new(
            first.position.x as f32,
            first.position.y as f32,
            first.position.z as f32,
        );
        let mut max = min;

        for v in &mesh.vertices {
            let p = Vec3::new(
                v.position.x as f32,
                v.position.y as f32,
                v.position.z as f32,
            );
            min = min.min(p);
            max = max.max(p);
        }

        Some((min, max))
    }
}

impl Default for VizServer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_has_no_model() {
        let server = VizServer::new();
        assert!(server.brep.is_none());
        assert!(server.mesh.is_none());
    }

    #[test]
    fn default_display_options() {
        let opts = DisplayOptions::default();
        assert!(opts.show_wireframe);
        assert!(!opts.show_normals);
        assert_eq!(opts.color_mode, ColorMode::Surface);
    }

    #[test]
    fn default_is_same_as_new() {
        let server = VizServer::default();
        assert!(server.brep.is_none());
        assert!(server.mesh.is_none());
    }

    #[test]
    fn model_bounds_none_when_no_mesh() {
        let server = VizServer::new();
        assert!(server.model_bounds().is_none());
    }

    #[test]
    fn select_updates_selection() {
        let mut server = VizServer::new();
        server.select(Selection::Face(FaceHandle::default()));
        assert!(matches!(server.selection, Selection::Face(_)));
        server.select(Selection::Edge(EdgeHandle::default()));
        assert!(matches!(server.selection, Selection::Edge(_)));
        server.select(Selection::Vertex(VertexHandle::default()));
        assert!(matches!(server.selection, Selection::Vertex(_)));
        server.select(Selection::None);
        assert!(matches!(server.selection, Selection::None));
    }

    #[test]
    fn model_bounds_with_mesh() {
        use popkern_core::geometry::{Point3 as P3, Vec3 as V3};
        use popkern_mesh::{MeshVertex, TriMesh};

        let mesh = TriMesh {
            zipped_t_junctions: 0,
            vertices: vec![
                MeshVertex {
                    position: P3::new(-1.0, 0.0, 2.0),
                    normal: V3::new(0.0, 0.0, 1.0),
                    uv: (0.0, 0.0),
                },
                MeshVertex {
                    position: P3::new(3.0, -2.0, 0.5),
                    normal: V3::new(0.0, 0.0, 1.0),
                    uv: (1.0, 0.0),
                },
                MeshVertex {
                    position: P3::new(1.0, 4.0, -1.0),
                    normal: V3::new(0.0, 0.0, 1.0),
                    uv: (0.0, 1.0),
                },
            ],
            triangles: vec![[0, 1, 2]],
            face_map: vec![FaceHandle::default()],
            face_meshes: vec![],
        };

        let mut server = VizServer::new();
        server.mesh = Some(mesh);

        let (min, max) = server.model_bounds().expect("should have bounds");
        assert!((min.x - (-1.0)).abs() < f32::EPSILON);
        assert!((min.y - (-2.0)).abs() < f32::EPSILON);
        assert!((min.z - (-1.0)).abs() < f32::EPSILON);
        assert!((max.x - 3.0).abs() < f32::EPSILON);
        assert!((max.y - 4.0).abs() < f32::EPSILON);
        assert!((max.z - 2.0).abs() < f32::EPSILON);
    }

    #[test]
    fn mesh_bounds_empty_mesh_is_none() {
        let mesh = TriMesh::default();
        assert!(VizServer::mesh_bounds(&mesh).is_none());
    }

    #[test]
    fn load_brep_tessellates_unit_cube() {
        use popkern_core::tolerance::Tolerance;
        let mut server = VizServer::new();
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube should build");
        let result = server.load_brep(brep);
        assert!(result.is_ok());
        assert!(server.brep.is_some());
        assert!(server.mesh.is_some());
        assert!(server.model_bounds().is_some());
    }

    #[test]
    fn color_mode_variants_are_distinct() {
        assert_ne!(ColorMode::Surface, ColorMode::Flat);
        assert_ne!(ColorMode::Flat, ColorMode::None);
        assert_ne!(ColorMode::Surface, ColorMode::None);
    }
}

#[cfg(test)]
#[path = "scene_tests.rs"]
mod scene_tests;
