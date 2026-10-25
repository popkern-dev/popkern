use super::*;

/// Every worker count gives the sequential map's results in order,
/// including an uneven last block, more workers than items and an
/// empty or offset range.
#[test]
fn every_worker_count_gives_the_sequential_results_in_order() {
    let f = |i: usize| (i as f64 * 0.1).sqrt().to_bits() ^ (i as u64 * 7919);
    for range in [0..1001, 17..40, 5..5, 3..4] {
        let sequential: Vec<u64> = range.clone().map(f).collect();
        for workers in [1, 2, 3, 8, 64] {
            assert_eq!(
                map_range_on(range.clone(), workers, f),
                sequential,
                "{range:?} on {workers}"
            );
        }
        assert_eq!(map_range(range.clone(), 1, f), sequential);
    }
}
