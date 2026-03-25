//! Render backend selection — GPU or software.

/// Available rendering backends.
pub enum RenderBackend {
    /// GPU-accelerated rendering via wgpu.
    #[cfg(feature = "gpu")]
    Gpu,
    /// CPU software rasterizer.
    Software,
}

/// Select the rendering backend.
///
/// If `force_software` is true, always returns [`RenderBackend::Software`].
/// Otherwise, attempts GPU first and falls back to software with a warning
/// if no GPU adapter is available.
pub fn resolve_backend(force_software: bool) -> RenderBackend {
    if force_software {
        return RenderBackend::Software;
    }

    #[cfg(feature = "gpu")]
    {
        if gpu_available() {
            return RenderBackend::Gpu;
        }
        log::warn!("No GPU adapter found, falling back to software renderer");
    }

    #[cfg(not(feature = "gpu"))]
    {
        log::info!("Built without GPU support, using software renderer");
    }

    RenderBackend::Software
}

/// Probe whether a wgpu GPU adapter is available.
#[cfg(feature = "gpu")]
fn gpu_available() -> bool {
    let instance = wgpu::Instance::default();
    pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        ..Default::default()
    }))
    .is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn software_flag_returns_software() {
        let backend = resolve_backend(true);
        assert!(matches!(backend, RenderBackend::Software));
    }

    #[test]
    fn no_flag_resolves_without_panic() {
        let backend = resolve_backend(false);
        match backend {
            #[cfg(feature = "gpu")]
            RenderBackend::Gpu => {}
            RenderBackend::Software => {}
        }
    }
}
