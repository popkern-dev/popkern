use super::{ChainTrim, discover_sub_chains};

type Trim = ChainTrim<u32, u32>;

/// Trim convention used throughout: fraction `f` along the trim
/// pairs side A's `a0 + (a1-a0)·f` with side B's
/// `b0 + (b1-b0)·(1-f)`; `node_start` sits at `f = 0`,
/// `node_end` at `f = 1`. Rows glued in the SAME chain direction
/// therefore encode side B's interval reversed.
fn trim(
    row_a: u32,
    interval_a: (f64, f64),
    row_b: u32,
    interval_b: (f64, f64),
    node_start: u32,
    node_end: u32,
) -> Trim {
    ChainTrim {
        row_a,
        interval_a,
        row_b,
        interval_b,
        node_start,
        node_end,
    }
}

/// Two full rows glued by one full-span trim: one sub-chain with
/// both rows embedded over the same interval (opposite
/// directions) and both endpoint nodes present.
#[test]
fn full_span_pair_forms_single_sub_chain() {
    let trims = vec![trim(0, (0.0, 1.0), 1, (0.0, 1.0), 10, 11)];
    let discovery = discover_sub_chains(&trims);
    assert_eq!(discovery.declined_chain_count, 0);
    assert_eq!(discovery.sub_chains.len(), 1);
    let sub = &discovery.sub_chains[0];
    assert_eq!(sub.rows.len(), 2);
    let (lo_a, hi_a) = sub.rows[0].1.chain_interval();
    let (lo_b, hi_b) = sub.rows[1].1.chain_interval();
    assert!((lo_a - lo_b).abs() < 1e-9 && (hi_a - hi_b).abs() < 1e-9);
    assert!(
        sub.rows[0].1.scale * sub.rows[1].1.scale < 0.0,
        "opposite traversal"
    );
    assert_eq!(sub.nodes.len(), 2);
}

/// A T-junction: one big row [0,1] faces two half rows. The
/// interior node (the T-stem at 0.5) is NOT clean (the big row
/// spans it), so the chain must stay one sub-chain containing
/// all three rows and three nodes.
#[test]
fn t_junction_keeps_one_sub_chain() {
    let trims = vec![
        trim(0, (0.0, 0.5), 1, (0.0, 1.0), 20, 21),
        trim(0, (0.5, 1.0), 2, (0.0, 1.0), 21, 22),
    ];
    let discovery = discover_sub_chains(&trims);
    assert_eq!(discovery.declined_chain_count, 0);
    assert_eq!(discovery.sub_chains.len(), 1, "T-stem node is not clean");
    let sub = &discovery.sub_chains[0];
    assert_eq!(sub.rows.len(), 3);
    assert_eq!(sub.nodes.len(), 3);
}

/// Offset both-partial overlap: row 0 covers chain [0,1], row 1
/// covers [0.5,1.5] in the same direction (interval_b reversed
/// in the encoding); their outer halves glue to small full rows.
/// No interior station is clean, so the chain stays whole.
#[test]
fn interleaved_offsets_stay_unsplit() {
    let trims = vec![
        // overlap: A[0.5,1] glued to B's [0,0.5] (same direction)
        trim(0, (0.5, 1.0), 1, (0.5, 0.0), 30, 31),
        // A's lower half vs small row 2 (full, same direction)
        trim(0, (0.0, 0.5), 2, (1.0, 0.0), 32, 30),
        // B's upper half vs small row 3 (full, same direction)
        trim(1, (0.5, 1.0), 3, (1.0, 0.0), 31, 33),
    ];
    let discovery = discover_sub_chains(&trims);
    assert_eq!(
        discovery.declined_chain_count, 0,
        "reason={}",
        discovery.decline_reason
    );
    assert_eq!(discovery.sub_chains.len(), 1);
    assert_eq!(discovery.sub_chains[0].rows.len(), 4);
}

/// A fully interleaved closed ring — outer rows break at 1.0 and
/// 2.0, inner rows at 0.5 and 1.5, period 2 — has no clean cut
/// and must decline; callers fall back to the per-trim
/// operations.
#[test]
fn interleaved_closed_ring_declines() {
    // Chain images (all same direction): A=[0,1], B=[1,2],
    // C=[0.5,1.5], D=[1.5,2.5]≡[1.5,2)∪[0,0.5).
    let trims = vec![
        trim(0, (0.5, 1.0), 2, (0.5, 0.0), 71, 72),
        trim(2, (0.5, 1.0), 1, (0.5, 0.0), 72, 73),
        trim(1, (0.5, 1.0), 3, (0.5, 0.0), 73, 74),
        trim(3, (0.5, 1.0), 0, (0.5, 0.0), 74, 71),
    ];
    let discovery = discover_sub_chains(&trims);
    assert_eq!(discovery.sub_chains.len(), 0);
    assert_eq!(discovery.declined_chain_count, 1);
    assert_eq!(discovery.decline_reason, "closed_chain_without_clean_cut");
}

/// A partially interleaved closed ring: stations 0.5/1.0/1.5 are
/// spanned by some row, but station 0 ≡ 2 is a breakpoint of
/// every row containing it. The period is detected through the
/// shared node id, the rebase cuts at the clean station, and one
/// sub-chain covers all five rows.
#[test]
fn partially_interleaved_ring_rebases_at_clean_station() {
    // Chain images: A=[0,1], B=[1,2], C=[0,0.5], D=[0.5,1.5],
    // E=[1.5,2]; period 2 (node 90 registered at 0 and at 2).
    let trims = vec![
        trim(0, (0.0, 0.5), 2, (1.0, 0.0), 90, 91),
        trim(0, (0.5, 1.0), 3, (0.5, 0.0), 91, 92),
        trim(1, (0.0, 0.5), 3, (1.0, 0.5), 92, 93),
        trim(1, (0.5, 1.0), 4, (1.0, 0.0), 93, 90),
    ];
    let discovery = discover_sub_chains(&trims);
    assert_eq!(
        discovery.declined_chain_count, 0,
        "reason={}",
        discovery.decline_reason
    );
    assert_eq!(discovery.sub_chains.len(), 1);
    assert_eq!(discovery.sub_chains[0].rows.len(), 5);
    // All five stations present as nodes (cut node duplicated at
    // the far end of the rebased window).
    assert!(discovery.sub_chains[0].nodes.len() >= 5);
}
