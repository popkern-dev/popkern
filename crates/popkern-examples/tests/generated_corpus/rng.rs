//! A small deterministic generator (SplitMix64): every scene is a pure
//! function of its cell, seed and index, on every platform.

/// SplitMix64 state.
#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    /// A generator for `(seed, stream, index)`: independent streams per
    /// cell and scene, so adding a cell never shifts another's scenes.
    pub fn new(seed: u64, stream: u64, index: u64) -> Rng {
        let mut r = Rng(seed ^ stream.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        for _ in 0..(index % 7) {
            r.next_u64();
        }
        r.0 ^= index.wrapping_mul(0xD1B5_4A32_D192_ED03);
        r.next_u64();
        r
    }

    /// The next 64 random bits.
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[lo, hi)`.
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// Uniform index below `n` (`n > 0`).
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n.max(1) as u64) as usize
    }

    /// A uniformly distributed unit vector.
    pub fn direction(&mut self) -> [f64; 3] {
        loop {
            let v = [
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
            ];
            let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if (0.1..=1.0).contains(&n) {
                return [v[0] / n, v[1] / n, v[2] / n];
            }
        }
    }
}

/// A stable 64-bit hash of a string (FNV-1a), for stream ids.
pub fn stream_id(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
        (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3)
    })
}
