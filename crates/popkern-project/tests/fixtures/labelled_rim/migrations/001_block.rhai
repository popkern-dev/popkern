let m = model();
let block = m.profile(xy()).rect(2.0, 2.0).sketch("block_sketch").extrude(blind(1.0));
let bore = m.tool(|m| m.profile(xy()).circle(0.4).sketch("bore_sketch").extrude(blind(1.0)));
block.subtract(bore);
m.label_edge("bore_rim", edges().circular().at_height(1.0));
m
