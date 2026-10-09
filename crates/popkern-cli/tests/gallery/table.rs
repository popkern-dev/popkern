//! The gallery table (split from `gallery.rs` for the file-size limit).

/// The gallery: file stem, title, blurb, expected feature count, and
/// the surface areas the script demonstrates (the guide cites them).
pub(crate) const GALLERY: &[(&str, &str, &str, usize, &str)] = &[
    (
        "boss_plate",
        "Boss plate",
        "A plate, a boss, a union and G1 rounds selected by geometry — the reference part.",
        6,
        "extrude, tool bodies, union, fillet with selectors, assert_watertight",
    ),
    (
        "revolved_knob",
        "Revolved knob",
        "A stepped section revolved about Z; exact planes, cylinders and rounded rims.",
        3,
        "path profiles, revolve, circular edge selection, fillet",
    ),
    (
        "shelled_enclosure",
        "Shelled enclosure",
        "A rounded block hollowed from the top to a uniform wall.",
        4,
        "fillet on vertical edges, shell with an open face",
    ),
    (
        "chamfered_bracket",
        "Chamfered bracket",
        "An L-bracket on the YZ plane with equal, two-distance and distance-angle chamfers.",
        5,
        "sketch planes other than XY, chamfer kinds, allow_empty",
    ),
    (
        "bored_block",
        "Bored block",
        "A through bore whose exact rims are labelled, then rounded by that label.",
        7,
        "subtract, exact boolean rims, label_edge, labeled() selection",
    ),
    (
        "tapered_fillet",
        "Tapered fillet",
        "A variable-radius round growing along one edge under a linear law.",
        3,
        "fillet radius laws",
    ),
    (
        "pyramid_apex",
        "Pyramid apex",
        "A pentagonal pyramid carved by five exact planes; its slant edges \
         rounded with curvature continuity (g2), and the five rounds closed \
         by the setback Gregory vertex cap with dual cross-jet corners.",
        18,
        "boolean carve with exact planes, chord-length edge selection, \
         vertex blend cap",
    ),
    (
        "mixed_radius_corner",
        "Mixed-radius corner",
        "A block corner where rounds of three different radii meet, closed \
         by the variational corner fit: three NURBS sectors certified G2.",
        3,
        "fillet_each radius groups, at_vertex edge selection, \
         variational corner fit",
    ),
    (
        "parametric_plate",
        "Parametric plate",
        "Parameters, an expression and a named configuration drive a plate and its rib.",
        6,
        "param, param_expr, p(), expr(), configurations",
    ),
    (
        "project_box",
        "Project box",
        "A printable two-body box: shelled base with screw bosses from a loop, a bar-lip lid with loop-cut vents, exact-cone countersinks and a vents-off configuration.",
        71,
        "loops computing placements, shell, multi-body documents, configurations suppressing loop features, chord_length_between",
    ),
    (
        "ball_bearing",
        "Ball bearing",
        "A deep-groove bearing assembly: races revolved with exact torus grooves and chamfered lead-ins, and one ball built once and placed eight times by trigonometry.",
        8,
        "assemblies (part, place, translate), revolve with arc paths, exact torus and sphere surfaces, trigonometric placement loops, product-structure STEP export",
    ),
    (
        "involute_gear",
        "Involute gear",
        "A twenty-tooth spur gear whose flanks are true involutes computed in the script, bored and keyed in one combined cut.",
        8,
        "computed polygon profiles, script functions, boolean tool unions",
    ),
    (
        "centrifugal_impeller",
        "Centrifugal impeller",
        "Seven backswept blades carved by one star-shaped cut whose spiral walls and tapered trailing edges this script computes.",
        11,
        "script functions, computed star polygons, logarithmic spirals, single-cut cavities",
    ),
    (
        "polyhedral_dice",
        "Polyhedral dice",
        "The five Platonic solids as a dice set — face normals computed from the \
         golden ratio, each die carved from a block by one exact plane cut per \
         face pair, and a filleted d6 dimpled with 21 spherical pips.",
        176,
        "plane_at custom sketch planes, intersect, computed face normals, \
         spherical dimple cuts, multi-body sets",
    ),
    (
        "coil_spring",
        "Coil spring",
        "Three springs from one verb: a ground-end compression spring, a tapered valve \
         spring, and a left-hand spring phased a quarter turn — every path an exact \
         rational helix.",
        12,
        "helical_sweep, exact helix paths, ground flat ends, taper_to, left_hand winding, \
         start_angle phasing",
    ),
    (
        "propeller",
        "Marine propeller",
        "Three lofted blades — stadium sections whose chord and pitch the script \
         computes per station — twisted, phased 120° apart and unioned onto a \
         revolved hub with a tapered nose.",
        26,
        "loft through computed sections, rational stadium profiles, rigid rotations, \
         chart-wrapping boolean unions on closed lofted walls",
    ),
    (
        "keyed_shaft",
        "Keyed shaft",
        "A shaft with a keyway cut by a moved and turned tool, and a chamfered lead-in.",
        8,
        "translate, rotate, subtract, chamfer on a circular edge",
    ),
    (
        "planetary_gearbox",
        "Planetary gearbox",
        "A 4:1 planetary set as an assembly — sun, one planet placed three times at its \
         mesh phase, an internal ring gear cut by a toothed disc, and a spider carrier — \
         whose involute teeth and planet phasing are all computed in the script.",
        35,
        "computed involute teeth, internal ring gears, meshing-assembly phasing, \
         the (Zs+Zr)%N constraint, single-polygon spiders, assemblies (one part placed \
         many times with rotate then translate)",
    ),
    (
        "threaded_stud",
        "Threaded stud",
        "The smallest real-thread part: an M10×1.5 stud from one verb — a coarse, \
         short coil that serves in seconds, with the tip left square for speed.",
        1,
        "threaded_rod, ISO thread designations, tip_chamfer option",
    ),
    (
        "hex_bolt",
        "Hex bolt",
        "A single short M6 hex bolt — head, plain shank and a real modeled ISO \
         thread from one verb, sized to build fast.",
        1,
        "bolt generator, thread_length, tip_chamfer option",
    ),
    (
        "countersunk_joint",
        "Countersunk joint",
        "A plate with a flush exact-cone clearance seat and the flat-head M5 screw \
         that fits it — one quick thread coil in the whole build.",
        5,
        "clearance_hole with countersink, countersunk head bolts, axis_at placement",
    ),
    (
        "fastener_kit",
        "Fastener kit",
        "An M6 hex bolt with its nut, a socket-head cap screw, and a countersunk \
         screw as an assembly — the three screws are one parametric part family — \
         with real modeled ISO threads, heads, tip lead-ins, and runout from the \
         nominal tables.",
        4,
        "bolt and nut generators, ISO thread designations, head styles, \
         printed-fit clearance, part families (one parametric part, several sizes)",
    ),
    (
        "angled_fitting",
        "Angled fitting",
        "A wedge bracket tapped normal to its inclined face, a countersunk \
         clearance hole through the base, and the flat-head screw that fits it.",
        6,
        "tapped_hole and clearance_hole on arbitrary axes, axis_at placement, \
         countersinks, thread fit classes",
    ),
    (
        "door_hinge",
        "Door hinge",
        "A 76 butt hinge: two leaves with five knuckles, a separate pin, and \
         three countersunk M3 screw holes a leaf, opened to 150°.",
        3,
        "hinge generators (butt), screw holes from a designation, sizes that \
         follow the length, the open angle",
    ),
    (
        "piano_hinge",
        "Piano hinge",
        "A continuous hinge: 25 knuckles along 200, one comb cut a leaf \
         whatever the count.",
        3,
        "hinge generators at high knuckle counts",
    ),
    (
        "living_hinges",
        "Living hinges",
        "A straight flexure with rounded roots for flexible plastics, and a \
         slotted serpentine that lets a rigid plate bend.",
        3,
        "living_hinge styles, single-profile construction",
    ),
    (
        "hinged_box",
        "Hinged box",
        "A shelled box and its lid joined by a print-in-place hinge grown \
         out of their back walls along the seam.",
        9,
        "hinge_bodies on existing bodies, print-in-place pins, shells",
    ),
    (
        "hinge_assembly",
        "Hinge assembly",
        "A print-in-place hinge as an assembly: the leaves mated about the \
         pin axis, the open angle a slider, and a swing check for the first \
         collision.",
        2,
        "assembly hinges, mates driven by a parameter, limits, swing",
    ),
    (
        "hinged_case",
        "Hinged case",
        "The hinged box as an assembly: box and lid are separate parts \
         hinged along their seam, and the lid opens on a slider without \
         rebuilding either part.",
        9,
        "assembly seam hinges (a.hinge_bodies), mates on a literal axis, \
         a driven open angle, swing",
    ),
    (
        "cycloidal_joint",
        "Cycloidal joint modules",
        "The three joint sizes of a six-axis arm (L, M, S), each a motor driving a \
         cycloidal reducer: finned housing, slotted stator and windings, a rotor \
         with 180°-opposed cams, two biarc-fitted cycloid discs riding the cams \
         round the ring pins, output pins, and a cross-roller bearing with its \
         rollers at ±45°. One slider per module turns the output; rotor, discs \
         and roller cage follow it.",
        426,
        "computed cycloid flanks fitted with G1 biarcs (deviation checked against \
         the backlash), revolute joints with followers at gear ratios, \
         sub-assemblies, part appearance, finned castings with rounded roots",
    ),
];
