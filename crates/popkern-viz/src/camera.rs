//! Camera module for the debug viewer.
//!
//! Provides a perspective camera with orbit, pan, zoom controls and
//! named view presets (front, back, top, bottom, left, right, iso).

use glam::{Mat4, Quat, Vec3};

/// Named orthographic view directions plus isometric.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewPreset {
    /// Looking along −Z (camera on +Z side).
    Front,
    /// Looking along +Z (camera on −Z side).
    Back,
    /// Looking along −Y (camera above).
    Top,
    /// Looking along +Y (camera below).
    Bottom,
    /// Looking along +X (camera on −X side).
    Left,
    /// Looking along −X (camera on +X side).
    Right,
    /// Isometric view from (+1, +1, +1) direction.
    Iso,
}

/// Perspective camera with orbit / pan / zoom controls.
#[derive(Debug, Clone)]
pub struct Camera {
    /// Camera position in world space.
    pub eye: Vec3,
    /// The point the camera looks at.
    pub target: Vec3,
    /// Up direction (not necessarily unit length after mutations).
    pub up: Vec3,
    /// Vertical field-of-view in **radians**.
    pub fov_y: f32,
    /// Viewport aspect ratio (width / height).
    pub aspect: f32,
    /// Near clipping plane distance.
    pub near: f32,
    /// Far clipping plane distance.
    pub far: f32,
}

/// Minimum allowed distance between eye and target during zoom.
const MIN_DISTANCE: f32 = 0.01;

/// Default distance used by presets (matches isometric default magnitude).
const PRESET_DISTANCE: f32 = 5.196;

impl Default for Camera {
    fn default() -> Self {
        Self::default_iso()
    }
}

impl Camera {
    /// Create the default isometric camera.
    ///
    /// * eye = (3, 3, 3)
    /// * target = (0, 0, 0)
    /// * up = (0, 1, 0)
    /// * fov_y = 45 degrees
    pub fn default_iso() -> Self {
        Self {
            eye: Vec3::new(3.0, 3.0, 3.0),
            target: Vec3::ZERO,
            up: Vec3::Y,
            fov_y: std::f32::consts::FRAC_PI_4,
            aspect: 16.0 / 9.0,
            near: 0.1,
            far: 1000.0,
        }
    }

    /// Compute the right-handed view (world → camera) matrix.
    pub fn view_matrix(&self) -> Mat4 {
        Mat4::look_at_rh(self.eye, self.target, self.up)
    }

    /// Compute the right-handed perspective projection matrix.
    pub fn projection_matrix(&self) -> Mat4 {
        Mat4::perspective_rh(self.fov_y, self.aspect, self.near, self.far)
    }

    /// Combined view-projection matrix (projection * view).
    pub fn view_projection(&self) -> Mat4 {
        self.projection_matrix() * self.view_matrix()
    }

    /// Orbit the camera around the target point.
    ///
    /// `dx` rotates around the world up axis (yaw).
    /// `dy` rotates around the camera right axis (pitch).
    /// Both values are in radians.
    pub fn orbit(&mut self, dx: f32, dy: f32) {
        let offset = self.eye - self.target;
        let distance = offset.length();
        if distance < f32::EPSILON {
            return;
        }

        // Yaw: rotate around world Y
        let yaw = Quat::from_rotation_y(-dx);
        let offset = yaw * offset;

        // Pitch: rotate around the camera-local right axis
        let forward = (-offset).normalize();
        let right = forward.cross(self.up).normalize();
        let pitch = Quat::from_axis_angle(right, -dy);
        let offset = pitch * offset;

        // Prevent flipping: reject if the new direction is too close to up
        let new_forward = (-offset).normalize();
        if new_forward.dot(Vec3::Y).abs() > 0.99 {
            // Only apply yaw
            let offset = yaw * (self.eye - self.target);
            self.eye = self.target + offset;
        } else {
            self.eye = self.target + offset;
        }

        // Keep distance stable (avoid drift from floating-point)
        let dir = (self.eye - self.target).normalize();
        self.eye = self.target + dir * distance;
    }

    /// Pan the camera (translate both eye and target).
    ///
    /// `dx` moves along the camera-local right axis.
    /// `dy` moves along the camera-local up axis.
    pub fn pan(&mut self, dx: f32, dy: f32) {
        let forward = (self.target - self.eye).normalize();
        let right = forward.cross(self.up).normalize();
        let cam_up = right.cross(forward).normalize();

        let offset = right * dx + cam_up * dy;
        self.eye += offset;
        self.target += offset;
    }

    /// Zoom by moving the eye toward or away from the target.
    ///
    /// Positive `delta` moves closer; negative moves farther.
    /// The eye-to-target distance is clamped to [`MIN_DISTANCE`].
    pub fn zoom(&mut self, delta: f32) {
        let offset = self.eye - self.target;
        let distance = offset.length();
        if distance < f32::EPSILON {
            return;
        }

        let new_distance = (distance - delta).max(MIN_DISTANCE);
        let dir = offset.normalize();
        self.eye = self.target + dir * new_distance;
    }

    /// Snap the camera to a named preset view direction.
    ///
    /// Preserves the current eye-to-target distance and target position.
    pub fn set_view_preset(&mut self, preset: ViewPreset) {
        let distance = (self.eye - self.target).length();
        let dist = if distance < f32::EPSILON {
            PRESET_DISTANCE
        } else {
            distance
        };

        let (dir, up) = match preset {
            ViewPreset::Front => (Vec3::Z, Vec3::Y),
            ViewPreset::Back => (Vec3::NEG_Z, Vec3::Y),
            ViewPreset::Top => (Vec3::Y, Vec3::NEG_Z),
            ViewPreset::Bottom => (Vec3::NEG_Y, Vec3::Z),
            ViewPreset::Left => (Vec3::NEG_X, Vec3::Y),
            ViewPreset::Right => (Vec3::X, Vec3::Y),
            ViewPreset::Iso => (Vec3::new(1.0, 1.0, 1.0).normalize(), Vec3::Y),
        };

        self.eye = self.target + dir * dist;
        self.up = up;
    }

    /// Adjust the camera so the given axis-aligned bounding box fills the
    /// viewport.
    ///
    /// The camera direction is preserved; only the target and eye distance
    /// change.
    pub fn fit_to_bounds(&mut self, min: Vec3, max: Vec3) {
        let center = (min + max) * 0.5;
        let half_extent = (max - min) * 0.5;
        let radius = half_extent.length();

        // Move target to box centre, keep look direction
        let dir = (self.eye - self.target).normalize();
        self.target = center;

        // Distance so the sphere fits in the vertical FOV
        let half_fov = self.fov_y * 0.5;
        let sin_half = half_fov.sin();
        let distance = if sin_half.abs() < f32::EPSILON {
            radius * 10.0
        } else {
            radius / sin_half
        };

        self.eye = center + dir * distance;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_has_expected_values() {
        let cam = Camera::default_iso();
        assert_eq!(cam.eye, Vec3::new(3.0, 3.0, 3.0));
        assert_eq!(cam.target, Vec3::ZERO);
        assert_eq!(cam.up, Vec3::Y);
    }
}
