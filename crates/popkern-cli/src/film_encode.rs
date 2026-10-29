//! Native film encoders (cli-render-export.md §5): the page's presets
//! (`film-export.js` PRESETS) as PNG (`image`), GIF (`gif`, a NeuQuant
//! palette per frame) and MP4 H.264 (OpenH264 behind the `h264`
//! feature, muxed with `mp4`). WebM stays the browser's (`--browser`).

use std::io::Cursor;

/// An export preset.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Preset {
    /// The preset's name.
    pub name: &'static str,
    /// Output kind.
    pub format: Format,
    /// A fixed frame size, else the film's (or `width` wide).
    pub size: Option<[u32; 2]>,
    /// A fixed width keeping the film's aspect.
    pub width: Option<u32>,
    /// A frame rate other than the film's.
    pub fps: Option<u32>,
    /// Video bit rate, bits per second.
    pub bitrate: u32,
}

/// What a preset writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// MP4 / H.264.
    Mp4,
    /// GIF89a.
    Gif,
    /// One PNG frame.
    Png,
    /// WebM / VP9 (the browser path only).
    Webm,
}

/// The page's preset table.
pub const PRESETS: &[Preset] = &[
    Preset {
        name: "mp4-1080p",
        format: Format::Mp4,
        size: Some([1920, 1080]),
        width: None,
        fps: None,
        bitrate: 8_000_000,
    },
    Preset {
        name: "mp4-vertical",
        format: Format::Mp4,
        size: Some([1080, 1920]),
        width: None,
        fps: None,
        bitrate: 8_000_000,
    },
    Preset {
        name: "mp4-720p",
        format: Format::Mp4,
        size: Some([1280, 720]),
        width: None,
        fps: None,
        bitrate: 4_000_000,
    },
    Preset {
        name: "webm",
        format: Format::Webm,
        size: Some([1920, 1080]),
        width: None,
        fps: None,
        bitrate: 6_000_000,
    },
    Preset {
        name: "gif",
        format: Format::Gif,
        size: None,
        width: Some(640),
        fps: Some(15),
        bitrate: 0,
    },
    Preset {
        name: "png",
        format: Format::Png,
        size: None,
        width: None,
        fps: None,
        bitrate: 0,
    },
];

/// The preset named `name`.
pub fn preset(name: &str) -> Option<Preset> {
    PRESETS.iter().copied().find(|p| p.name == name)
}

impl Preset {
    /// The frame size for a film of `film_size` (`presetSize`): a GIF
    /// keeps the film's aspect at its width, its height even.
    pub fn frame_size(&self, film_size: [u32; 2]) -> [u32; 2] {
        if let Some(size) = self.size {
            return size;
        }
        if let Some(w) = self.width {
            let h = ((w as f64 * film_size[1] as f64 / film_size[0].max(1) as f64 / 2.0).round()
                * 2.0) as u32;
            return [w, h.max(2)];
        }
        film_size
    }
}

/// PNG bytes of an RGBA frame.
pub fn png(rgba: &[u8], width: u32, height: u32) -> Result<Vec<u8>, String> {
    let image = image::RgbaImage::from_raw(width, height, rgba.to_vec())
        .ok_or("the frame does not match its size")?;
    let mut out = Cursor::new(Vec::new());
    image
        .write_to(&mut out, image::ImageFormat::Png)
        .map_err(|e| format!("png: {e}"))?;
    Ok(out.into_inner())
}

/// A GIF89a being written frame by frame, looping forever.
pub struct GifWriter {
    encoder: gif::Encoder<Vec<u8>>,
    width: u16,
    height: u16,
    delay: u16,
}

impl GifWriter {
    /// A GIF of `width × height` at `fps`.
    pub fn new(width: u32, height: u32, fps: u32) -> Result<GifWriter, String> {
        let (w, h) = (
            u16::try_from(width).map_err(|_| "gif too wide")?,
            u16::try_from(height).map_err(|_| "gif too tall")?,
        );
        let mut encoder =
            gif::Encoder::new(Vec::new(), w, h, &[]).map_err(|e| format!("gif: {e}"))?;
        encoder
            .set_repeat(gif::Repeat::Infinite)
            .map_err(|e| format!("gif: {e}"))?;
        Ok(GifWriter {
            encoder,
            width: w,
            height: h,
            // Hundredths of a second, as gifenc rounds 1000/fps ms.
            delay: ((1000.0 / fps.max(1) as f64).round() / 10.0).round() as u16,
        })
    }

    /// Append an RGBA frame.
    pub fn frame(&mut self, rgba: &[u8]) -> Result<(), String> {
        let mut pixels = rgba.to_vec();
        let mut frame = gif::Frame::from_rgba_speed(self.width, self.height, &mut pixels, 10);
        frame.delay = self.delay;
        self.encoder
            .write_frame(&frame)
            .map_err(|e| format!("gif: {e}"))
    }

    /// The finished file.
    pub fn finish(self) -> Result<Vec<u8>, String> {
        self.encoder.into_inner().map_err(|e| format!("gif: {e}"))
    }
}

/// The NAL units of an Annex-B stream, start codes removed.
#[cfg(any(feature = "h264", test))]
fn nal_units(stream: &[u8]) -> Vec<&[u8]> {
    let mut starts = Vec::new();
    let mut i = 0;
    while i + 3 <= stream.len() {
        if stream[i] == 0 && stream[i + 1] == 0 && stream[i + 2] == 1 {
            starts.push(i + 3);
            i += 3;
        } else {
            i += 1;
        }
    }
    let mut out = Vec::new();
    for (k, &s) in starts.iter().enumerate() {
        let mut end = starts.get(k + 1).map_or(stream.len(), |&n| n - 3);
        // A four-byte start code's leading zero belongs to it.
        while end > s && stream[end - 1] == 0 && starts.get(k + 1).is_some() {
            end -= 1;
        }
        if end > s {
            out.push(&stream[s..end]);
        }
    }
    out
}

/// One encoded access unit: AVCC bytes and whether it is a key frame.
struct Sample {
    bytes: Vec<u8>,
    sync: bool,
}

/// An MP4 / H.264 film being written frame by frame.
pub struct Mp4Writer {
    #[cfg(feature = "h264")]
    encoder: openh264::encoder::Encoder,
    width: u32,
    height: u32,
    fps: u32,
    sps: Option<Vec<u8>>,
    pps: Option<Vec<u8>>,
    samples: Vec<Sample>,
}

impl Mp4Writer {
    /// An H.264 stream of `width × height` (both even) at `fps` and
    /// `bitrate`, a key frame every two seconds.
    #[cfg(feature = "h264")]
    pub fn new(width: u32, height: u32, fps: u32, bitrate: u32) -> Result<Mp4Writer, String> {
        use openh264::encoder::{
            BitRate, Encoder, EncoderConfig, FrameRate, IntraFramePeriod, RateControlMode,
        };
        if !width.is_multiple_of(2) || !height.is_multiple_of(2) {
            return Err(format!(
                "an MP4 frame must be even-sized, not {width}x{height}"
            ));
        }
        let config = EncoderConfig::new()
            .bitrate(BitRate::from_bps(bitrate))
            .max_frame_rate(FrameRate::from_hz(fps as f32))
            .rate_control_mode(RateControlMode::Bitrate)
            .skip_frames(false)
            .scene_change_detect(false)
            .intra_frame_period(IntraFramePeriod::from_num_frames(2 * fps));
        let encoder = Encoder::with_api_config(openh264::OpenH264API::from_source(), config)
            .map_err(|e| format!("h264 encoder: {e}"))?;
        Ok(Mp4Writer {
            encoder,
            width,
            height,
            fps,
            sps: None,
            pps: None,
            samples: Vec::new(),
        })
    }

    /// Built without H.264: MP4 is unavailable.
    #[cfg(not(feature = "h264"))]
    pub fn new(_width: u32, _height: u32, _fps: u32, _bitrate: u32) -> Result<Mp4Writer, String> {
        Err(
            "this popkern was built without the h264 feature: MP4 is unavailable \
             (use --preset gif or png, or --browser)"
                .into(),
        )
    }

    /// Append an RGBA frame.
    #[cfg(feature = "h264")]
    pub fn frame(&mut self, rgba: &[u8]) -> Result<(), String> {
        use openh264::encoder::FrameType;
        use openh264::formats::{RgbaSliceU8, YUVBuffer};
        let (w, h) = (self.width as usize, self.height as usize);
        let yuv = YUVBuffer::from_rgba8_source(RgbaSliceU8::new(rgba, (w, h)));
        let stream = self
            .encoder
            .encode(&yuv)
            .map_err(|e| format!("h264: {e}"))?;
        let sync = matches!(stream.frame_type(), FrameType::IDR | FrameType::I);
        let annex_b = stream.to_vec();
        let mut bytes = Vec::with_capacity(annex_b.len());
        for nal in nal_units(&annex_b) {
            match nal.first().map(|b| b & 0x1f) {
                Some(7) => {
                    self.sps.get_or_insert_with(|| nal.to_vec());
                }
                Some(8) => {
                    self.pps.get_or_insert_with(|| nal.to_vec());
                }
                Some(9) | None => {}
                Some(_) => {
                    bytes.extend_from_slice(&(nal.len() as u32).to_be_bytes());
                    bytes.extend_from_slice(nal);
                }
            }
        }
        if bytes.is_empty() {
            return Err(format!(
                "h264: frame {} produced no picture",
                self.samples.len()
            ));
        }
        self.samples.push(Sample { bytes, sync });
        Ok(())
    }

    /// Built without H.264: unreachable (construction fails).
    #[cfg(not(feature = "h264"))]
    pub fn frame(&mut self, _rgba: &[u8]) -> Result<(), String> {
        Err("h264 unavailable".into())
    }

    /// The finished file: one video track, a sample per frame.
    pub fn finish(self) -> Result<Vec<u8>, String> {
        use mp4::{AvcConfig, FourCC, MediaConfig, Mp4Config, Mp4Sample, TrackConfig, TrackType};
        let (Some(sps), Some(pps)) = (self.sps, self.pps) else {
            return Err("h264: the stream carried no parameter sets".into());
        };
        let config = Mp4Config {
            major_brand: FourCC::from(*b"isom"),
            minor_version: 512,
            compatible_brands: [*b"isom", *b"iso2", *b"avc1", *b"mp41"]
                .map(FourCC::from)
                .to_vec(),
            timescale: 1000,
        };
        let err = |e: mp4::Error| format!("mp4: {e}");
        let mut writer =
            mp4::Mp4Writer::write_start(Cursor::new(Vec::new()), &config).map_err(err)?;
        let timescale = self.fps * 1000;
        writer
            .add_track(&TrackConfig {
                track_type: TrackType::Video,
                timescale,
                language: "und".into(),
                media_conf: MediaConfig::AvcConfig(AvcConfig {
                    width: u16::try_from(self.width).map_err(|_| "mp4 too wide")?,
                    height: u16::try_from(self.height).map_err(|_| "mp4 too tall")?,
                    seq_param_set: sps,
                    pic_param_set: pps,
                }),
            })
            .map_err(err)?;
        for (k, sample) in self.samples.into_iter().enumerate() {
            writer
                .write_sample(
                    1,
                    &Mp4Sample {
                        start_time: k as u64 * 1000,
                        duration: 1000,
                        rendering_offset: 0,
                        is_sync: sample.sync,
                        bytes: bytes::Bytes::from(sample.bytes),
                    },
                )
                .map_err(err)?;
        }
        writer.write_end().map_err(err)?;
        Ok(writer.into_writer().into_inner())
    }
}

#[cfg(test)]
#[path = "film_encode_tests.rs"]
mod tests;
