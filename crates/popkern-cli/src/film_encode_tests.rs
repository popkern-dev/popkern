//! The encoders on small synthetic frames.

use super::*;

fn frame(w: u32, h: u32, k: u32) -> Vec<u8> {
    let mut out = Vec::with_capacity((w * h * 4) as usize);
    for y in 0..h {
        for x in 0..w {
            out.extend_from_slice(&[(x * 4 + k * 40) as u8, (y * 5) as u8, 128, 255]);
        }
    }
    out
}

#[test]
fn nal_units_split_both_start_code_lengths() {
    let stream = [
        0, 0, 0, 1, 0x67, 1, 2, 0, 0, 1, 0x68, 3, 0, 0, 0, 1, 0x65, 9, 9,
    ];
    let units = nal_units(&stream);
    assert_eq!(
        units,
        vec![&[0x67, 1, 2][..], &[0x68, 3][..], &[0x65, 9, 9][..]]
    );
}

#[test]
fn presets_size_frames_as_the_page() {
    assert_eq!(
        preset("mp4-vertical").map(|p| p.frame_size([1920, 1080])),
        Some([1080, 1920])
    );
    assert_eq!(
        preset("gif").map(|p| p.frame_size([1920, 1080])),
        Some([640, 360])
    );
    assert_eq!(
        preset("png").map(|p| p.frame_size([800, 600])),
        Some([800, 600])
    );
    assert!(preset("mov").is_none());
}

#[test]
fn png_and_gif_write_their_headers_and_frames() {
    let bytes = png(&frame(32, 24, 0), 32, 24).expect("png");
    let decoded = image::load_from_memory(&bytes).expect("decodes");
    assert_eq!((decoded.width(), decoded.height()), (32, 24));
    let mut gif = GifWriter::new(32, 24, 15).expect("gif");
    for k in 0..4 {
        gif.frame(&frame(32, 24, k)).expect("frame");
    }
    let bytes = gif.finish().expect("finish");
    assert!(bytes.starts_with(b"GIF89a"));
    let mut options = gif::DecodeOptions::new();
    options.set_color_output(gif::ColorOutput::RGBA);
    let mut decoder = options.read_info(Cursor::new(&bytes)).expect("reads");
    let mut frames = 0;
    while let Some(f) = decoder.read_next_frame().expect("frame") {
        assert_eq!(f.delay, 7, "15 fps: 67 ms → 7 cs");
        frames += 1;
    }
    assert_eq!(frames, 4);
}

#[cfg(feature = "h264")]
#[test]
fn mp4_holds_one_sample_per_frame_and_the_film_duration() {
    let (w, h, fps) = (64, 48, 10);
    let mut mp4 = Mp4Writer::new(w, h, fps, 500_000).expect("encoder");
    for k in 0..25 {
        mp4.frame(&frame(w, h, k)).expect("frame");
    }
    let bytes = mp4.finish().expect("finish");
    assert_eq!(&bytes[4..8], b"ftyp");
    let size = bytes.len() as u64;
    let reader = mp4::Mp4Reader::read_header(Cursor::new(bytes), size).expect("parses");
    let track = reader.tracks().get(&1).expect("video track");
    assert_eq!(track.sample_count(), 25);
    assert_eq!(track.width(), 64);
    assert!(
        (track.duration().as_secs_f64() - 2.5).abs() < 1e-6,
        "{:?}",
        track.duration()
    );
    // A key frame every two seconds: frames 0 and 20.
    let mut reader = reader;
    let syncs: Vec<u32> = (1..=25)
        .filter(|&id| {
            reader
                .read_sample(1, id)
                .ok()
                .flatten()
                .is_some_and(|s| s.is_sync)
        })
        .collect();
    assert_eq!(syncs, vec![1, 21]);
}

#[test]
fn an_odd_mp4_size_is_refused() {
    if cfg!(feature = "h264") {
        assert!(Mp4Writer::new(63, 48, 10, 500_000).is_err());
    }
}
