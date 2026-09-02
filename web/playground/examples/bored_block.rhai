// A block with a through bore: the boolean names the bore's rim
// edges, a label pins them, and the fillet selects by that label so a
// later edit of the bore keeps the rounds on it.
let m = model();
let block = m.profile(xy()).rect(1.2, 0.8).sketch("block_sketch").extrude(blind(0.5));
let bore = m.tool(|m| m.profile(xy_at(-0.1)).circle(0.18).sketch("bore_sketch").extrude(blind(0.7)));
block = block.subtract(bore);
m.label_edge("bore_rim", edges().circular().at_height(0.5));
block = block.fillet(edges().labeled("bore_rim"), 0.05, #{ label: "rim_round" });
block = block.fillet(edges().straight().direction("z"), 0.1, #{ label: "corner_rounds" });
m.assert_watertight();
m
