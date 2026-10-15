// A tabletop dice set: the classic seven — d4, d6, d8, d10, d%, d12
// and d20 — with engraved numerals, plus a pipped d6. Each die is
// carved from a block by one slab per pair of opposite faces (the d4's
// faces have no opposites, so its slabs reach from far below), its
// edges bevelled or rounded, and its numbers cut in.
//
//   EDGES         how the die edges are broken: "chamfer" (the
//                 default) bevels every edge by edge_chamfer, a
//                 fraction of a second per die; "round" rounds them
//                 with an edge_round fillet, which costs ~9 minutes
//                 (each corner of a d12 or d20 is a fitted patch, ~5 s
//                 apiece, where an exact rolling-ball sphere would do);
//                 "sharp" leaves them. The pipped d6 is always rounded:
//                 a cube's corners take a second.
//   depth         how deep the numerals are cut
//
// Numbering follows real dice: opposite faces sum to n + 1 (the d10 to
// 9, the d% to 90), a 6 and a 9 are told apart by an underline on the
// dice that carry both, the d4 reads at its corners (the number at the
// top vertex is the roll) and the d% reads in tens. The numerals are
// Inter Bold outlines (SIL Open Font License 1.1, rsms.me/inter),
// flattened to polygons by tools/glyphs/digit-polygons.py; each face's
// numerals are one tool body whose counters (the holes in 0, 4, 6, 8,
// 9) are cut out of it, so the counters stay standing in the face.
let m = model(#{ unit: "mm" });
m.param("edge_chamfer", 0.5, "mm");
m.param("edge_round", 0.7, "mm");
m.param("depth", 0.5, "mm");
let EDGES = "chamfer";
let pi = 3.141592653589793;
let phi = (1.0 + sqrt(5.0)) / 2.0;

// ---- Vectors ----
fn dot(a, b) { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn cross(a, b) { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
fn scale(a, s) { [a[0] * s, a[1] * s, a[2] * s] }
fn add(a, b) { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
fn sub(a, b) { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
fn norm(v) { scale(v, 1.0 / sqrt(dot(v, v))) }

// The in-plane axes `plane_at(o, n)` sketches in (the kernel's
// SketchPlane::from_normal): u = n x z (n x x when n is along z),
// v = n x u, turned so u x v = n.
fn plane_axes(n) {
    let c = if dot(cross(n, [0.0, 0.0, 1.0]), cross(n, [0.0, 0.0, 1.0])) > 1e-12 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] };
    let u = norm(cross(n, c));
    let v = norm(cross(n, u));
    if dot(cross(u, v), n) < 0.0 { v = scale(v, -1.0); }
    [u, v]
}

// ---- A die's polyhedron from its face normals ----
// Every face is the plane n . x = r. The vertices are where three
// faces meet inside every other face's half-space.
fn vertices_of(normals, r) {
    let out = [];
    let k = normals.len();
    for i in 0..k {
        for j in (i + 1)..k {
            for l in (j + 1)..k {
                let a = normals[i];
                let b = normals[j];
                let c = normals[l];
                let det = dot(a, cross(b, c));
                if det.abs() < 1e-9 { continue; }
                let x = scale(add(add(cross(b, c), cross(c, a)), cross(a, b)), r / det);
                let inside = true;
                for n in normals {
                    if dot(n, x) > r * (1.0 + 1e-9) { inside = false; break; }
                }
                if !inside { continue; }
                let fresh = true;
                for y in out {
                    let d = sub(x, y);
                    if dot(d, d) < r * r * 1e-12 { fresh = false; break; }
                }
                if fresh { out.push(x); }
            }
        }
    }
    out
}
// The vertices on face n (n . x = r).
fn face_vertices(verts, n, r) {
    let out = [];
    for x in verts {
        if (dot(n, x) - r).abs() < r * 1e-7 { out.push(x); }
    }
    out
}
// Both normals of each opposite pair, in pair order: [n0, -n0, n1, -n1, ...].
fn with_opposites(normals) {
    let out = [];
    for n in normals { out.push(norm(n)); out.push(scale(norm(n), -1.0)); }
    out
}

// Carve a die centred on the origin: material within r of each face.
fn carve(m, ctx, tag, normals, r, paired, big) {
    let die = m.profile(xy_at(-big)).rect(2.0 * big, 2.0 * big)
        .sketch(tag + "_block").extrude(blind(2.0 * big));
    let i = 0;
    for raw in normals {
        let n = norm(raw);
        let slab = if paired {
            m.tool(|t| t.profile(plane_at(scale(n, -r), n))
                .rect(4.0 * big, 4.0 * big).sketch(tag + "_slab_" + i)
                .extrude(blind(2.0 * r)))
        } else {
            let deep = 3.0 * big;
            m.tool(|t| t.profile(plane_at(scale(n, -deep), n))
                .rect(6.0 * big, 6.0 * big).sketch(tag + "_slab_" + i)
                .extrude(blind(deep + r)))
        };
        die = die.intersect(slab, #{ label: tag + "_cut_" + i });
        i += 1;
    }
    if ctx.edges == "chamfer" {
        die = die.chamfer(edges(), p("edge_chamfer"), #{ label: tag + "_bevels" });
    } else if ctx.edges == "round" {
        die = die.fillet(edges(), p("edge_round"), #{ label: tag + "_rounds" });
    }
    die
}

// ---- Numerals ----
let GLYPHS = #{
    "0": #{ adv: 0.9268, contours: [
        #{ hole: false, pts: [[0.5397, -0.0074], [0.5813, 0.0019], [0.6204, 0.0153], [0.6568, 0.0330], [0.6906, 0.0549], [0.7217, 0.0810], [0.7502, 0.1113], [0.7743, 0.1433], [0.7957, 0.1785], [0.8141, 0.2168], [0.8285, 0.2544], [0.8404, 0.2946], [0.8499, 0.3373], [0.8564, 0.3783], [0.8608, 0.4213], [0.8631, 0.4619], [0.8637, 0.5041], [0.8627, 0.5458], [0.8601, 0.5859], [0.8552, 0.6285], [0.8483, 0.6691], [0.8384, 0.7113], [0.8260, 0.7511], [0.8111, 0.7884], [0.7921, 0.8262], [0.7701, 0.8612], [0.7456, 0.8927], [0.7166, 0.9226], [0.6851, 0.9483], [0.6509, 0.9697], [0.6142, 0.9868], [0.5748, 0.9997], [0.5328, 1.0084], [0.4917, 1.0126], [0.4487, 1.0132], [0.4070, 1.0102], [0.3675, 1.0036], [0.3273, 0.9924], [0.2896, 0.9769], [0.2545, 0.9572], [0.2220, 0.9332], [0.1921, 0.9050], [0.1647, 0.8725], [0.1418, 0.8386], [0.1218, 0.8018], [0.1047, 0.7619], [0.0916, 0.7228], [0.0810, 0.6812], [0.0735, 0.6412], [0.0680, 0.5991], [0.0646, 0.5550], [0.0632, 0.5135], [0.0635, 0.4711], [0.0654, 0.4302], [0.0694, 0.3867], [0.0755, 0.3453], [0.0845, 0.3022], [0.0960, 0.2615], [0.1099, 0.2234], [0.1279, 0.1847], [0.1488, 0.1490], [0.1725, 0.1163], [0.2006, 0.0854], [0.2313, 0.0586], [0.2646, 0.0361], [0.3005, 0.0178], [0.3390, 0.0037], [0.3802, -0.0062], [0.4205, -0.0116], [0.4631, -0.0134]] },
        #{ hole: true, pts: [[0.4081, 0.1656], [0.3719, 0.1852], [0.3423, 0.2139], [0.3199, 0.2479], [0.3030, 0.2854], [0.2905, 0.3254], [0.2817, 0.3669], [0.2761, 0.4089], [0.2730, 0.4505], [0.2718, 0.4909], [0.2723, 0.5323], [0.2748, 0.5750], [0.2797, 0.6183], [0.2869, 0.6581], [0.2974, 0.6972], [0.3130, 0.7372], [0.3330, 0.7730], [0.3593, 0.8042], [0.3933, 0.8284], [0.4326, 0.8416], [0.4744, 0.8440], [0.5142, 0.8364], [0.5510, 0.8178], [0.5812, 0.7899], [0.6041, 0.7566], [0.6216, 0.7197], [0.6346, 0.6801], [0.6439, 0.6390], [0.6499, 0.5974], [0.6534, 0.5562], [0.6549, 0.5161], [0.6547, 0.4745], [0.6527, 0.4312], [0.6486, 0.3910], [0.6420, 0.3505], [0.6321, 0.3106], [0.6184, 0.2723], [0.5991, 0.2350], [0.5737, 0.2021], [0.5423, 0.1771], [0.5042, 0.1613], [0.4631, 0.1564]] }] },
    "1": #{ adv: 0.5926, contours: [
        #{ hole: false, pts: [[0.2832, 1.0000], [0.0624, 0.8450], [0.0624, 0.6577], [0.2980, 0.8195], [0.2980, 0.0000], [0.5027, 0.0000], [0.5027, 1.0000]] }] },
    "2": #{ adv: 0.8658, contours: [
        #{ hole: false, pts: [[0.7906, 0.0000], [0.7906, 0.1685], [0.3678, 0.1758], [0.5421, 0.3464], [0.5732, 0.3759], [0.6045, 0.4068], [0.6326, 0.4358], [0.6617, 0.4678], [0.6888, 0.5002], [0.7147, 0.5350], [0.7370, 0.5704], [0.7554, 0.6081], [0.7683, 0.6478], [0.7756, 0.6878], [0.7770, 0.7304], [0.7726, 0.7717], [0.7628, 0.8105], [0.7462, 0.8493], [0.7234, 0.8851], [0.6953, 0.9170], [0.6645, 0.9431], [0.6293, 0.9656], [0.5928, 0.9832], [0.5532, 0.9969], [0.5142, 1.0059], [0.4726, 1.0114], [0.4324, 1.0134], [0.3887, 1.0122], [0.3474, 1.0074], [0.3051, 0.9981], [0.2654, 0.9846], [0.2287, 0.9670], [0.1933, 0.9443], [0.1622, 0.9177], [0.1352, 0.8871], [0.1125, 0.8527], [0.0951, 0.8155], [0.0828, 0.7755], [0.0761, 0.7356], [0.0738, 0.6933], [0.2707, 0.7145], [0.2789, 0.7547], [0.2978, 0.7917], [0.3270, 0.8212], [0.3626, 0.8400], [0.4025, 0.8487], [0.4425, 0.8488], [0.4818, 0.8410], [0.5196, 0.8227], [0.5504, 0.7945], [0.5703, 0.7582], [0.5781, 0.7185], [0.5766, 0.6783], [0.5661, 0.6387], [0.5464, 0.6018], [0.5217, 0.5699], [0.4954, 0.5396], [0.4670, 0.5094], [0.4376, 0.4799], [0.0832, 0.1477], [0.0832, 0.0000]] }] },
    "3": #{ adv: 0.8872, contours: [
        #{ hole: false, pts: [[0.4935, -0.0113], [0.5347, -0.0062], [0.5773, 0.0028], [0.6176, 0.0155], [0.6548, 0.0314], [0.6907, 0.0513], [0.7250, 0.0760], [0.7549, 0.1041], [0.7801, 0.1354], [0.8009, 0.1714], [0.8151, 0.2099], [0.8227, 0.2510], [0.8239, 0.2912], [0.8193, 0.3318], [0.8076, 0.3715], [0.7890, 0.4075], [0.7638, 0.4395], [0.7323, 0.4664], [0.6970, 0.4865], [0.6589, 0.5007], [0.6280, 0.5267], [0.6656, 0.5410], [0.7013, 0.5618], [0.7323, 0.5885], [0.7577, 0.6206], [0.7753, 0.6566], [0.7854, 0.6967], [0.7879, 0.7387], [0.7835, 0.7814], [0.7724, 0.8218], [0.7544, 0.8599], [0.7301, 0.8950], [0.7006, 0.9260], [0.6682, 0.9513], [0.6312, 0.9730], [0.5936, 0.9894], [0.5533, 1.0015], [0.5135, 1.0089], [0.4714, 1.0128], [0.4305, 1.0132], [0.3899, 1.0105], [0.3477, 1.0038], [0.3076, 0.9934], [0.2696, 0.9792], [0.2324, 0.9606], [0.1968, 0.9370], [0.1654, 0.9097], [0.1384, 0.8786], [0.1156, 0.8428], [0.0993, 0.8045], [0.0894, 0.7635], [0.0859, 0.7228], [0.2864, 0.7492], [0.3030, 0.7866], [0.3319, 0.8154], [0.3683, 0.8356], [0.4068, 0.8466], [0.4490, 0.8496], [0.4903, 0.8444], [0.5290, 0.8287], [0.5601, 0.8029], [0.5810, 0.7681], [0.5891, 0.7284], [0.5853, 0.6867], [0.5685, 0.6495], [0.5400, 0.6202], [0.5050, 0.6003], [0.4661, 0.5893], [0.4242, 0.5859], [0.3315, 0.5859], [0.3315, 0.4336], [0.4242, 0.4336], [0.4648, 0.4312], [0.5054, 0.4228], [0.5438, 0.4061], [0.5765, 0.3803], [0.5996, 0.3468], [0.6100, 0.3076], [0.6085, 0.2675], [0.5942, 0.2290], [0.5684, 0.1982], [0.5347, 0.1757], [0.4968, 0.1622], [0.4559, 0.1567], [0.4144, 0.1578], [0.3742, 0.1658], [0.3367, 0.1819], [0.3041, 0.2066], [0.2811, 0.2410], [0.2725, 0.2805], [0.0687, 0.2479], [0.0767, 0.2085], [0.0916, 0.1689], [0.1131, 0.1322], [0.1389, 0.1004], [0.1693, 0.0721], [0.2041, 0.0475], [0.2404, 0.0279], [0.2776, 0.0126], [0.3172, 0.0008], [0.3592, -0.0075], [0.3998, -0.0119], [0.4423, -0.0134]] }] },
    "4": #{ adv: 0.9295, contours: [
        #{ hole: false, pts: [[0.5463, 0.1779], [0.5463, 0.0000], [0.7430, 0.0000], [0.7430, 0.1779], [0.8678, 0.1779], [0.8678, 0.3430], [0.7430, 0.3430], [0.7430, 1.0000], [0.4859, 1.0000], [0.0678, 0.3396], [0.0678, 0.1779]] },
        #{ hole: true, pts: [[0.2745, 0.3510], [0.5503, 0.7745], [0.5503, 0.3430]] }] },
    "5": #{ adv: 0.8785, contours: [
        #{ hole: false, pts: [[0.5034, -0.0101], [0.5441, -0.0034], [0.5856, 0.0081], [0.6246, 0.0240], [0.6612, 0.0442], [0.6944, 0.0679], [0.7239, 0.0950], [0.7516, 0.1278], [0.7735, 0.1618], [0.7910, 0.1984], [0.8038, 0.2372], [0.8119, 0.2784], [0.8152, 0.3184], [0.8146, 0.3597], [0.8095, 0.4023], [0.7999, 0.4425], [0.7858, 0.4803], [0.7656, 0.5182], [0.7413, 0.5521], [0.7130, 0.5819], [0.6807, 0.6075], [0.6447, 0.6289], [0.6060, 0.6449], [0.5649, 0.6555], [0.5246, 0.6604], [0.4837, 0.6606], [0.4428, 0.6557], [0.4021, 0.6448], [0.3633, 0.6272], [0.3300, 0.6045], [0.2993, 0.5779], [0.3255, 0.8309], [0.7604, 0.8309], [0.7604, 1.0000], [0.1557, 1.0000], [0.1074, 0.4752], [0.2963, 0.4490], [0.3275, 0.4757], [0.3631, 0.4938], [0.4030, 0.5053], [0.4441, 0.5087], [0.4840, 0.5045], [0.5237, 0.4907], [0.5580, 0.4674], [0.5849, 0.4373], [0.6045, 0.4001], [0.6146, 0.3597], [0.6166, 0.3174], [0.6113, 0.2773], [0.5968, 0.2376], [0.5733, 0.2035], [0.5417, 0.1760], [0.5048, 0.1580], [0.4649, 0.1498], [0.4239, 0.1501], [0.3836, 0.1589], [0.3458, 0.1776], [0.3148, 0.2049], [0.2938, 0.2405], [0.2859, 0.2805], [0.0875, 0.2625], [0.0937, 0.2221], [0.1065, 0.1813], [0.1258, 0.1433], [0.1495, 0.1101], [0.1780, 0.0803], [0.2108, 0.0541], [0.2452, 0.0331], [0.2837, 0.0153], [0.3218, 0.0025], [0.3622, -0.0065], [0.4047, -0.0118], [0.4456, -0.0134]] }] },
    "6": #{ adv: 0.8926, contours: [
        #{ hole: false, pts: [[0.5375, -0.0072], [0.5771, 0.0014], [0.6176, 0.0149], [0.6555, 0.0327], [0.6903, 0.0543], [0.7238, 0.0812], [0.7531, 0.1119], [0.7783, 0.1465], [0.7980, 0.1821], [0.8131, 0.2202], [0.8234, 0.2606], [0.8291, 0.3033], [0.8301, 0.3444], [0.8268, 0.3871], [0.8190, 0.4274], [0.8056, 0.4682], [0.7869, 0.5064], [0.7635, 0.5414], [0.7361, 0.5723], [0.7047, 0.5990], [0.6694, 0.6214], [0.6313, 0.6386], [0.5908, 0.6501], [0.5510, 0.6559], [0.5091, 0.6568], [0.4681, 0.6527], [0.4290, 0.6436], [0.3890, 0.6282], [0.3535, 0.6084], [0.3218, 0.5838], [0.2936, 0.5530], [0.2645, 0.5253], [0.2660, 0.5668], [0.2702, 0.6090], [0.2777, 0.6512], [0.2893, 0.6925], [0.3057, 0.7320], [0.3268, 0.7671], [0.3541, 0.7973], [0.3880, 0.8199], [0.4271, 0.8331], [0.4688, 0.8369], [0.5094, 0.8323], [0.5483, 0.8176], [0.5806, 0.7930], [0.6041, 0.7606], [0.6181, 0.7228], [0.8190, 0.7261], [0.8110, 0.7674], [0.7980, 0.8066], [0.7802, 0.8437], [0.7578, 0.8782], [0.7303, 0.9103], [0.6988, 0.9386], [0.6658, 0.9614], [0.6298, 0.9805], [0.5919, 0.9952], [0.5514, 1.0055], [0.5118, 1.0113], [0.4700, 1.0134], [0.4281, 1.0118], [0.3881, 1.0063], [0.3471, 0.9960], [0.3083, 0.9811], [0.2717, 0.9617], [0.2376, 0.9378], [0.2062, 0.9098], [0.1778, 0.8778], [0.1541, 0.8447], [0.1329, 0.8082], [0.1157, 0.7717], [0.1006, 0.7325], [0.0881, 0.6910], [0.0790, 0.6511], [0.0719, 0.6092], [0.0669, 0.5653], [0.0641, 0.5241], [0.0631, 0.4814], [0.0638, 0.4369], [0.0663, 0.3950], [0.0710, 0.3517], [0.0779, 0.3113], [0.0879, 0.2704], [0.1010, 0.2313], [0.1175, 0.1926], [0.1380, 0.1548], [0.1611, 0.1213], [0.1897, 0.0895], [0.2222, 0.0613], [0.2569, 0.0378], [0.2939, 0.0191], [0.3337, 0.0048], [0.3746, -0.0053], [0.4164, -0.0114], [0.4591, -0.0134]] },
        #{ hole: true, pts: [[0.3788, 0.1676], [0.3456, 0.1906], [0.3187, 0.2207], [0.2992, 0.2558], [0.2879, 0.2956], [0.2854, 0.3357], [0.2912, 0.3762], [0.3058, 0.4136], [0.3290, 0.4477], [0.3590, 0.4752], [0.3945, 0.4951], [0.4334, 0.5056], [0.4756, 0.5067], [0.5170, 0.4983], [0.5536, 0.4802], [0.5839, 0.4539], [0.6085, 0.4200], [0.6241, 0.3817], [0.6305, 0.3414], [0.6289, 0.2992], [0.6185, 0.2586], [0.5991, 0.2223], [0.5712, 0.1910], [0.5384, 0.1680], [0.4992, 0.1535], [0.4577, 0.1490]] }] },
    "7": #{ adv: 0.7993, contours: [
        #{ hole: false, pts: [[0.3315, 0.0000], [0.7490, 0.8282], [0.7490, 1.0000], [0.0503, 1.0000], [0.0503, 0.8309], [0.5349, 0.8242], [0.1188, 0.0000]] }] },
    "8": #{ adv: 0.8946, contours: [
        #{ hole: false, pts: [[0.5024, -0.0111], [0.5437, -0.0059], [0.5829, 0.0022], [0.6234, 0.0144], [0.6612, 0.0299], [0.6980, 0.0495], [0.7331, 0.0738], [0.7636, 0.1016], [0.7893, 0.1325], [0.8101, 0.1679], [0.8241, 0.2058], [0.8312, 0.2460], [0.8314, 0.2885], [0.8246, 0.3298], [0.8106, 0.3687], [0.7896, 0.4053], [0.7633, 0.4381], [0.7320, 0.4664], [0.6965, 0.4896], [0.6580, 0.5064], [0.6219, 0.5240], [0.6610, 0.5369], [0.6962, 0.5563], [0.7274, 0.5820], [0.7547, 0.6141], [0.7755, 0.6502], [0.7889, 0.6892], [0.7949, 0.7312], [0.7938, 0.7725], [0.7859, 0.8117], [0.7701, 0.8509], [0.7483, 0.8851], [0.7202, 0.9166], [0.6893, 0.9424], [0.6539, 0.9647], [0.6167, 0.9824], [0.5786, 0.9958], [0.5382, 1.0054], [0.4957, 1.0112], [0.4546, 1.0134], [0.4134, 1.0124], [0.3707, 1.0079], [0.3302, 0.9999], [0.2889, 0.9871], [0.2507, 0.9704], [0.2142, 0.9491], [0.1822, 0.9243], [0.1529, 0.8939], [0.1294, 0.8603], [0.1118, 0.8218], [0.1018, 0.7805], [0.0994, 0.7393], [0.1036, 0.6987], [0.1154, 0.6583], [0.1348, 0.6210], [0.1608, 0.5878], [0.1927, 0.5597], [0.2294, 0.5387], [0.2617, 0.5135], [0.2218, 0.5009], [0.1847, 0.4823], [0.1509, 0.4578], [0.1220, 0.4288], [0.0981, 0.3956], [0.0797, 0.3585], [0.0680, 0.3190], [0.0632, 0.2771], [0.0653, 0.2350], [0.0742, 0.1954], [0.0899, 0.1583], [0.1124, 0.1235], [0.1394, 0.0933], [0.1711, 0.0665], [0.2046, 0.0447], [0.2422, 0.0258], [0.2806, 0.0111], [0.3217, -0.0003], [0.3614, -0.0076], [0.4032, -0.0120], [0.4470, -0.0134]] },
        #{ hole: true, pts: [[0.3805, 0.1475], [0.3432, 0.1642], [0.3109, 0.1907], [0.2880, 0.2251], [0.2769, 0.2655], [0.2772, 0.3064], [0.2883, 0.3460], [0.3105, 0.3801], [0.3416, 0.4071], [0.3790, 0.4256], [0.4197, 0.4347], [0.4620, 0.4358], [0.5016, 0.4295], [0.5407, 0.4139], [0.5742, 0.3893], [0.5996, 0.3574], [0.6147, 0.3183], [0.6180, 0.2780], [0.6109, 0.2375], [0.5918, 0.2011], [0.5624, 0.1722], [0.5266, 0.1522], [0.4870, 0.1413], [0.4470, 0.1383]] },
        #{ hole: true, pts: [[0.3724, 0.6030], [0.3394, 0.6262], [0.3143, 0.6589], [0.3013, 0.6978], [0.3002, 0.7396], [0.3104, 0.7791], [0.3336, 0.8134], [0.3654, 0.8378], [0.4041, 0.8525], [0.4446, 0.8570], [0.4849, 0.8535], [0.5238, 0.8402], [0.5577, 0.8157], [0.5820, 0.7821], [0.5933, 0.7434], [0.5933, 0.7022], [0.5818, 0.6633], [0.5585, 0.6304], [0.5259, 0.6055], [0.4871, 0.5908], [0.4470, 0.5866]] }] },
    "9": #{ adv: 0.8926, contours: [
        #{ hole: false, pts: [[0.5019, -0.0082], [0.5432, 0.0018], [0.5821, 0.0164], [0.6188, 0.0356], [0.6532, 0.0593], [0.6848, 0.0871], [0.7134, 0.1189], [0.7373, 0.1519], [0.7587, 0.1882], [0.7761, 0.2245], [0.7914, 0.2636], [0.8041, 0.3050], [0.8135, 0.3449], [0.8208, 0.3867], [0.8260, 0.4305], [0.8289, 0.4717], [0.8302, 0.5144], [0.8298, 0.5552], [0.8276, 0.5974], [0.8233, 0.6412], [0.8168, 0.6821], [0.8072, 0.7235], [0.7948, 0.7623], [0.7787, 0.8015], [0.7586, 0.8398], [0.7359, 0.8740], [0.7084, 0.9060], [0.6763, 0.9352], [0.6419, 0.9595], [0.6052, 0.9790], [0.5666, 0.9937], [0.5250, 1.0049], [0.4826, 1.0119], [0.4395, 1.0147], [0.3957, 1.0135], [0.3542, 1.0086], [0.3151, 1.0000], [0.2751, 0.9866], [0.2378, 0.9689], [0.2030, 0.9472], [0.1695, 0.9201], [0.1402, 0.8891], [0.1149, 0.8543], [0.0953, 0.8184], [0.0802, 0.7801], [0.0698, 0.7397], [0.0642, 0.6972], [0.0632, 0.6563], [0.0665, 0.6135], [0.0743, 0.5731], [0.0866, 0.5351], [0.1049, 0.4966], [0.1280, 0.4611], [0.1550, 0.4299], [0.1860, 0.4028], [0.2207, 0.3801], [0.2579, 0.3628], [0.2977, 0.3509], [0.3400, 0.3443], [0.3816, 0.3431], [0.4227, 0.3468], [0.4619, 0.3556], [0.5020, 0.3708], [0.5381, 0.3907], [0.5706, 0.4157], [0.5980, 0.4450], [0.6294, 0.4715], [0.6277, 0.4300], [0.6233, 0.3878], [0.6155, 0.3457], [0.6036, 0.3044], [0.5868, 0.2650], [0.5654, 0.2302], [0.5377, 0.2002], [0.5034, 0.1779], [0.4640, 0.1651], [0.4221, 0.1618], [0.3823, 0.1668], [0.3442, 0.1819], [0.3125, 0.2071], [0.2893, 0.2402], [0.2758, 0.2785], [0.0740, 0.2719], [0.0825, 0.2304], [0.0958, 0.1911], [0.1140, 0.1540], [0.1367, 0.1197], [0.1644, 0.0876], [0.1963, 0.0593], [0.2297, 0.0365], [0.2660, 0.0174], [0.3043, 0.0029], [0.3450, -0.0074], [0.3849, -0.0129], [0.4268, -0.0148]] },
        #{ hole: true, pts: [[0.3572, 0.5102], [0.3226, 0.5335], [0.2949, 0.5643], [0.2754, 0.6000], [0.2647, 0.6401], [0.2626, 0.6818], [0.2687, 0.7220], [0.2841, 0.7607], [0.3076, 0.7942], [0.3378, 0.8210], [0.3739, 0.8402], [0.4140, 0.8498], [0.4543, 0.8501], [0.4956, 0.8409], [0.5325, 0.8215], [0.5627, 0.7948], [0.5862, 0.7617], [0.6019, 0.7241], [0.6085, 0.6832], [0.6063, 0.6418], [0.5948, 0.6016], [0.5743, 0.5654], [0.5468, 0.5354], [0.5130, 0.5122], [0.4754, 0.4980], [0.4342, 0.4933]] }] }
}
;

// Glyph contours of `text` laid out on the pen line, cap height 1:
// [[hole, points], ...], centred on the origin, with an underline when
// `mark` (the 6 and 9 of a die that has both).
fn text_contours(glyphs, text, mark) {
    let out = [];
    let pen = 0.0;
    for ch in text.chars() {
        let g = glyphs[ch.to_string()];
        // A glyph's outline first, then the counters inside it (every
        // digit has one outline).
        for hole in [false, true] {
            for c in g.contours {
                if c.hole != hole { continue; }
                let pts = [];
                for q in c.pts { pts.push([q[0] + pen, q[1]]); }
                out.push([c.hole, pts]);
            }
        }
        pen += g.adv;
    }
    // Centre the ink: x over its extent, y on the cap height.
    let lo = 1e9;
    let hi = -1e9;
    for c in out { for q in c[1] { if q[0] < lo { lo = q[0]; } if q[0] > hi { hi = q[0]; } } }
    let dx = -0.5 * (lo + hi);
    let dy = if mark { -0.38 } else { -0.5 };
    let centred = [];
    for c in out {
        let pts = [];
        for q in c[1] { pts.push([q[0] + dx, q[1] + dy]); }
        centred.push([c[0], pts]);
    }
    if mark {
        let w = 0.42 * (hi - lo);
        let y0 = -0.62;
        centred.push([false, [[-w, y0 - 0.1], [w, y0 - 0.1], [w, y0], [-w, y0]]]);
    }
    #{ contours: centred, width: hi - lo }
}

// Cut `text` into the face of `die` with outward normal n, through
// point `at` on the face, reading upright toward `up` (a point on the
// face), cap height h. Each glyph is one cutter: its outline prism with
// its counters (the holes in 0, 4, 6, 8, 9) cut out of it, so the
// counters stay standing in the pocket.
fn engrave(m, ctx, die, tag, n, at, up, text, h) {
    let mark = ctx.mark && (text == "6" || text == "9");
    let laid = text_contours(ctx.glyphs, text, mark);
    // Two-digit numbers shrink to the face's width budget.
    let s = h;
    if laid.width * s > ctx.width * h { s = ctx.width * h / laid.width; }
    let w = sub(up, at);
    w = norm(sub(w, scale(n, dot(w, n))));
    // Reading direction across the face, seen from outside (n toward
    // the viewer, w up: right = w x n).
    let right = cross(w, n);
    // Sketch planes face into the die so the cut runs inward.
    let inward = scale(n, -1.0);
    let axes = plane_axes(inward);
    let a = [dot(w, axes[0]), dot(w, axes[1])];
    let b = [dot(right, axes[0]), dot(right, axes[1])];
    let to_plane = |q| [s * (q[0] * b[0] + q[1] * a[0]), s * (q[0] * b[1] + q[1] * a[1])];
    // The ink starts 0.2 above the face so its walls cross it cleanly;
    // counters overshoot the ink at both ends.
    let lift = 0.2;
    let top = add(at, scale(n, lift));
    let high = add(at, scale(n, 2.0 * lift));
    // Group each outline with the counters inside it (text_contours
    // lists a glyph's outline before its counters).
    let glyphs = [];
    for c in laid.contours {
        let pts = [];
        for q in c[1] { pts.push(to_plane.call(q)); }
        if c[0] && glyphs.len() > 0 { glyphs[glyphs.len() - 1].holes.push(pts); }
        else { glyphs.push(#{ outer: pts, holes: [] }); }
    }
    let k = 0;
    for g in glyphs {
        let outer = g.outer;
        let holes = g.holes;
        let cutter = m.tool(|t| {
            let ink = t.profile(plane_at(top, inward)).polygon(outer)
                .sketch(tag + "_ink_" + k).extrude(blind(expr("depth + 0.2mm")));
            let j = 0;
            for pts in holes {
                ink = ink.subtract(t.profile(plane_at(high, inward)).polygon(pts)
                    .sketch(tag + "_counter_" + k + "_" + j).extrude(blind(expr("depth + 0.6mm"))),
                    #{ label: tag + "_counter_cut_" + k + "_" + j });
                j += 1;
            }
            ink
        });
        die = die.subtract(cutter, #{ label: tag + "_cut_" + k });
        k += 1;
    }
    die
}

// Rest a die on the table: turn face normal `down` to -z, then set it
// on z = 0 at (x, y).
fn rest(m, die, tag, down, r, x, y) {
    let d = norm(down);
    let c = -d[2];
    if c >= 1.0 { c = 1.0 }
    if c <= -1.0 { c = -1.0 }
    let ang = acos(c) * 180.0 / 3.141592653589793;
    let axis = [-d[1], d[0], 0.0];
    if sqrt(axis[0] * axis[0] + axis[1] * axis[1]) < 1e-12 { axis = [1.0, 0.0, 0.0]; }
    if ang > 1e-9 {
        die = die.rotate(axis_at([0.0, 0.0, 0.0], axis), ang, #{ label: tag + "_rest" });
    }
    die.translate(x, y, r, #{ label: tag + "_place" })
}

// A die with one numeral per face: values[i] on normals[i], read
// toward the face vertex farthest from its centre (the kite's apex on
// a d10, any corner elsewhere — the first found).
fn numbered(m, ctx, tag, normals, paired, r, values, h) {
    let all = if paired { with_opposites(normals) } else { normals.map(|n| norm(n)) };
    let die = carve(m, ctx, tag, normals, r, paired, 2.0 * r);
    let verts = vertices_of(all, r);
    for i in 0..all.len() {
        let n = all[i];
        let c = scale(n, r);
        let fv = face_vertices(verts, n, r);
        let up = fv[0];
        let far = 0.0;
        for v in fv { let d = sub(v, c); if dot(d, d) > far + 1e-9 * r * r { far = dot(d, d); up = v; } }
        die = engrave(m, ctx, die, tag + "_" + values[i], n, c, up, values[i], h);
    }
    die
}

// A d4 reads at its corners: every face carries the numbers of its
// three vertices, each near its vertex and upright toward it, so the
// number at the top vertex reads the same on all three upper faces.
fn corner_numbered(m, ctx, tag, normals, r, h) {
    let all = normals.map(|n| norm(n));
    let die = carve(m, ctx, tag, normals, r, false, 3.0 * r);
    let verts = vertices_of(all, r);
    for i in 0..all.len() {
        let n = all[i];
        let c = scale(n, r);
        for j in 0..verts.len() {
            let v = verts[j];
            if (dot(n, v) - r).abs() > r * 1e-7 { continue; }
            let at = add(c, scale(sub(v, c), 0.5));
            die = engrave(m, ctx, die, tag + "_f" + i + "_" + (j + 1), n, at, v, (j + 1).to_string(), h);
        }
    }
    die
}

// The pentagonal trapezohedron's face normals: five upper faces at
// elevation `tilt`, the lower five turned 36 degrees (each lower face
// is the antipode of an upper one, so one slab cuts a pair).
fn d10_normals(tilt) {
    let out = [];
    let ce = cos(tilt * 3.141592653589793 / 180.0);
    let se = sin(tilt * 3.141592653589793 / 180.0);
    for k in 0..5 {
        let a = 2.0 * 3.141592653589793 * k / 5.0;
        out.push([ce * cos(a), ce * sin(a), se]);
    }
    out
}

// ---- The set ----
let ctx = #{ glyphs: GLYPHS, mark: false, width: 1.3, edges: EDGES };
let marked = ctx;
marked.mark = true;

// Values for paired normals [n0, -n0, n1, -n1, ...]: opposite faces sum
// to n + 1.
fn opposite_pairs(count, total) {
    let out = [];
    for k in 1..=count { out.push(k.to_string()); out.push((total - k).to_string()); }
    out
}

// d4: tetrahedron.
let r4 = 4.2;
let d4 = corner_numbered(m, ctx, "d4",
    [[-1.0, -1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]], r4, 2.2);
d4 = rest(m, d4, "d4", [-1.0, -1.0, -1.0], r4, 0.0, 0.0);

// d6: cube, numerals.
let r6 = 8.0;
let d6 = numbered(m, ctx, "d6", [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], true, r6,
    ["1", "6", "2", "5", "3", "4"], 7.0);
d6 = rest(m, d6, "d6", [0.0, 0.0, -1.0], r6, 26.0, 0.0);

// d8: octahedron, faces toward the cube corners.
let r8 = 6.0;
let d8 = numbered(m, ctx, "d8",
    [[1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [-1.0, 1.0, 1.0], [-1.0, -1.0, 1.0]], true, r8,
    opposite_pairs(4, 9), 4.6);
d8 = rest(m, d8, "d8", [-1.0, -1.0, -1.0], r8, 52.0, 0.0);

// d10 and d%: pentagonal trapezohedra, 0..9 and 00..90 (opposite
// faces sum to 9 and 90).
let r10 = 7.6;
let tilt = 42.0;
let d10_values = [];
let pct_values = [];
for k in 0..5 {
    d10_values.push(k.to_string()); d10_values.push((9 - k).to_string());
    pct_values.push(if k == 0 { "00" } else { (k * 10).to_string() });
    pct_values.push(((9 - k) * 10).to_string());
}
let d10 = numbered(m, marked, "d10", d10_normals(tilt), true, r10, d10_values, 4.4);
d10 = rest(m, d10, "d10", scale(d10_normals(tilt)[0], -1.0), r10, 0.0, 28.0);
let pct = numbered(m, ctx, "dpct", d10_normals(tilt), true, r10, pct_values, 3.6);
pct = rest(m, pct, "dpct", scale(d10_normals(tilt)[0], -1.0), r10, 26.0, 28.0);

// d12: dodecahedron, faces toward the icosahedron's vertices.
let r12 = 8.6;
let d12 = numbered(m, marked, "d12",
    [[0.0, 1.0, phi], [0.0, -1.0, phi], [1.0, phi, 0.0],
     [1.0, -phi, 0.0], [phi, 0.0, 1.0], [phi, 0.0, -1.0]], true, r12,
    opposite_pairs(6, 13), 4.8);
d12 = rest(m, d12, "d12", [0.0, -1.0, -phi], r12, 52.0, 28.0);

// d20: icosahedron, faces toward the dodecahedron's vertices.
let r20 = 8.8;
let d20 = numbered(m, marked, "d20",
    [[1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [1.0, -1.0, 1.0], [1.0, -1.0, -1.0],
     [0.0, 1.0 / phi, phi], [0.0, -1.0 / phi, phi], [1.0 / phi, phi, 0.0],
     [-1.0 / phi, phi, 0.0], [phi, 0.0, 1.0 / phi], [phi, 0.0, -1.0 / phi]], true, r20,
    opposite_pairs(10, 21), 3.6);
d20 = rest(m, d20, "d20", [-1.0, -1.0, -1.0], r20, 78.0, 28.0);

// ---- The pipped d6 ----
// One pip: a sphere (a semicircle revolved about z) set just outside
// the face so its cap dips in as a dimple.
fn pip(m, tag, center, rp) {
    m.tool(|t| {
        let ball = t.profile(xz())
            .path(0.0, -rp).arc_to([rp, 0.0], [0.0, rp]).line_to(0.0, -rp)
            .close().sketch(tag)
            .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
        ball.translate(center[0], center[1], center[2])
    })
}
// Pip layouts on a face's (u, v) grid.
fn pips_for(value, g) {
    switch value {
        1 => [[0.0, 0.0]],
        2 => [[-1.0 * g, -1.0 * g], [g, g]],
        3 => [[-1.0 * g, -1.0 * g], [0.0, 0.0], [g, g]],
        4 => [[-1.0 * g, -1.0 * g], [-1.0 * g, g], [g, -1.0 * g], [g, g]],
        5 => [[-1.0 * g, -1.0 * g], [-1.0 * g, g], [0.0, 0.0], [g, -1.0 * g], [g, g]],
        _ => [[-1.0 * g, -1.0 * g], [-1.0 * g, 0.0], [-1.0 * g, g],
              [g, -1.0 * g], [g, 0.0], [g, g]],
    }
}
let half = 8.0;
let pipped = m.profile(xy_at(-half)).rect(2.0 * half, 2.0 * half).sketch("pip_d6_block")
    .extrude(blind(2.0 * half));
// A cube's corners are exact spheres: always rounded.
pipped = pipped.fillet(edges(), p("edge_round"), #{ label: "pip_d6_rounds" });
let rp = 2.1;            // pip sphere radius
let sink = 1.5;          // sphere centre this far outside the face
let g = 4.4;             // pip grid pitch
// Each face: [normal, u axis, v axis, value]; opposite faces sum to 7.
let faces6 = [
    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1],
    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0], 6],
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 3],
    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], 4],
    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2],
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 5],
];
for face in faces6 {
    let n = face[0];
    let u = face[1];
    let v = face[2];
    let value = face[3];
    let k = 0;
    for uv in pips_for(value, g) {
        let c = add(add(scale(n, half + sink), scale(u, uv[0])), scale(v, uv[1]));
        let tag = "pip_" + value + "_" + k;
        pipped = pipped.subtract(pip(m, tag, c, rp), #{ label: tag + "_cut" });
        k += 1;
    }
}
pipped = pipped.translate(78.0, 0.0, half, #{ label: "pip_d6_place" });

m.assert_watertight();
m
