//! Pixel format conversion for softbuffer presentation.

/// Convert an RGBA pixel `[R, G, B, A]` to softbuffer's native u32 (0x00RRGGBB).
pub fn rgba_to_u32(rgba: [u8; 4]) -> u32 {
    (rgba[0] as u32) << 16 | (rgba[1] as u32) << 8 | rgba[2] as u32
}

/// Convert an RGBA byte buffer to a u32 buffer for softbuffer.
pub fn rgba_buf_to_u32_buf(rgba: &[u8]) -> Vec<u32> {
    rgba.chunks_exact(4)
        .map(|px| rgba_to_u32([px[0], px[1], px[2], px[3]]))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rgba_to_native_pixel() {
        assert_eq!(rgba_to_u32([255, 128, 64, 255]), 0x00FF_8040);
    }

    #[test]
    fn convert_buffer_two_pixels() {
        let rgba = vec![255, 0, 0, 255, 0, 255, 0, 255];
        let native = rgba_buf_to_u32_buf(&rgba);
        assert_eq!(native.len(), 2);
        assert_eq!(native[0], 0x00FF_0000);
        assert_eq!(native[1], 0x0000_FF00);
    }

    #[test]
    fn black_pixel() {
        assert_eq!(rgba_to_u32([0, 0, 0, 255]), 0x0000_0000);
    }

    #[test]
    fn white_pixel() {
        assert_eq!(rgba_to_u32([255, 255, 255, 255]), 0x00FF_FFFF);
    }

    #[test]
    fn alpha_is_ignored() {
        // Two pixels with different alpha but same RGB produce identical output.
        assert_eq!(
            rgba_to_u32([100, 200, 50, 0]),
            rgba_to_u32([100, 200, 50, 255])
        );
    }

    #[test]
    fn empty_buffer_produces_empty_vec() {
        let native = rgba_buf_to_u32_buf(&[]);
        assert!(native.is_empty());
    }

    #[test]
    fn partial_trailing_bytes_ignored() {
        // chunks_exact(4) drops the trailing 2 bytes.
        let rgba = vec![255, 0, 0, 255, 0, 128];
        let native = rgba_buf_to_u32_buf(&rgba);
        assert_eq!(native.len(), 1);
        assert_eq!(native[0], 0x00FF_0000);
    }
}
