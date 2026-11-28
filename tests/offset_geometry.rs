use open_evdop::{Arithmetic, OffsetCurve, OffsetError, OffsetParameters, Offsetter, Orientation};
use std::f64::consts::{PI, TAU};

const LINE_TOLERANCE: f64 = 1e-9;
const SAMPLES_PER_CURVE: usize = 16;

#[derive(Debug, Clone, Copy)]
enum Input {
    Line {
        start: (f64, f64),
        end: (f64, f64),
    },
    Arc {
        start: (f64, f64),
        end: (f64, f64),
        center: (f64, f64),
        orientation: Orientation,
    },
}

fn dist(a: (f64, f64), b: (f64, f64)) -> f64 {
    ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt()
}

fn angle_of(p: (f64, f64), c: (f64, f64)) -> f64 {
    (p.1 - c.1).atan2(p.0 - c.0)
}

fn ccw_sweep(from: f64, to: f64) -> f64 {
    (to - from).rem_euclid(TAU)
}

fn input_distance(input: &Input, p: (f64, f64)) -> f64 {
    match *input {
        Input::Line { start, end } => {
            let (dx, dy) = (end.0 - start.0, end.1 - start.1);
            let t = (((p.0 - start.0) * dx + (p.1 - start.1) * dy) / (dx * dx + dy * dy))
                .clamp(0.0, 1.0);
            dist(p, (start.0 + t * dx, start.1 + t * dy))
        }
        Input::Arc {
            start,
            end,
            center,
            orientation,
        } => {
            let (from, to) = if orientation == Orientation::CounterClockwise {
                (start, end)
            } else {
                (end, start)
            };
            let (a0, a1, a) = (
                angle_of(from, center),
                angle_of(to, center),
                angle_of(p, center),
            );
            if ccw_sweep(a0, a) <= ccw_sweep(a0, a1) {
                (dist(p, center) - dist(start, center)).abs()
            } else {
                dist(p, start).min(dist(p, end))
            }
        }
    }
}

fn nearest_input_distance(inputs: &[Input], p: (f64, f64)) -> f64 {
    inputs
        .iter()
        .map(|i| input_distance(i, p))
        .fold(f64::INFINITY, f64::min)
}

/// Points along an output curve; arcs are swept in their reported orientation, except a vertex arc
/// at a tangent joint, whose ends lie within rounding of each other and which has no sweep.
fn curve_samples(curve: &OffsetCurve) -> Vec<(f64, f64)> {
    (0..=SAMPLES_PER_CURVE)
        .map(|k| {
            let t = k as f64 / SAMPLES_PER_CURVE as f64;
            match *curve {
                OffsetCurve::Line { start, end } => (
                    start.0 + t * (end.0 - start.0),
                    start.1 + t * (end.1 - start.1),
                ),
                OffsetCurve::Arc {
                    start,
                    end,
                    center,
                    orientation,
                } => {
                    let (a0, a1) = (angle_of(start, center), angle_of(end, center));
                    let ccw = ccw_sweep(a0, a1);
                    let sweep = match orientation {
                        _ if dist(start, end) <= LINE_TOLERANCE => 0.0,
                        Orientation::CounterClockwise => ccw,
                        Orientation::Clockwise => ccw - TAU,
                    };
                    let r = dist(start, center);
                    let a = a0 + t * sweep;
                    (center.0 + r * a.cos(), center.1 + r * a.sin())
                }
            }
        })
        .collect()
}

fn same_point(a: (f64, f64), b: (f64, f64)) -> bool {
    (a.0.to_bits(), a.1.to_bits()) == (b.0.to_bits(), b.1.to_bits())
}

fn curve_endpoints(curve: &OffsetCurve) -> ((f64, f64), (f64, f64)) {
    match *curve {
        OffsetCurve::Line { start, end } | OffsetCurve::Arc { start, end, .. } => (start, end),
    }
}

fn run(inputs: &[Input], params: OffsetParameters) -> Offsetter {
    let mut v = Offsetter::new();
    v.parameters = params;
    for input in inputs {
        match *input {
            Input::Line { start, end } => {
                v.add_line(start, end);
            }
            Input::Arc {
                start,
                end,
                center,
                orientation,
            } => {
                v.add_arc(start, end, center, orientation);
            }
        }
    }
    v.compute_offset().expect("compute_offset");
    v
}

const ARITHMETICS: [Arithmetic; 2] = [Arithmetic::Double, Arithmetic::Extended];

/// Checks the offsets in every arithmetic, which must give the same number of loops; returns the
/// offset distance of each loop in `Double`.
fn check_offsets(
    name: &str,
    inputs: &[Input],
    params: OffsetParameters,
    tolerance: f64,
) -> Vec<f64> {
    let mut per_arithmetic: Vec<Vec<f64>> = ARITHMETICS
        .iter()
        .map(|&arithmetic| {
            let label = format!("{name} [{arithmetic:?}]");
            check_offsets_in(
                &label,
                inputs,
                OffsetParameters {
                    arithmetic,
                    ..params
                },
                tolerance,
            )
        })
        .collect();
    assert_eq!(
        per_arithmetic[0].len(),
        per_arithmetic[1].len(),
        "{name}: the arithmetics give different numbers of loops"
    );
    per_arithmetic.swap_remove(0)
}

/// Checks every loop is closed and every sampled point lies at the expected offset distance.
/// Returns the offset distance of each loop.
fn check_offsets_in(
    name: &str,
    inputs: &[Input],
    params: OffsetParameters,
    tolerance: f64,
) -> Vec<f64> {
    let v = run(inputs, params);
    let mut errors = Vec::new();
    let mut loop_dists = Vec::new();
    println!("{name}: {} loops", v.num_offset_loops());
    for l in 0..v.num_offset_loops() {
        let n = v.num_offset_curves(l).expect("loop exists");
        let curves: Vec<OffsetCurve> = (0..n)
            .map(|c| v.curve(l, c).expect("curve exists"))
            .collect();
        let first = curve_samples(&curves[0])[0];
        let expected = nearest_input_distance(inputs, first);
        loop_dists.push(expected);
        println!("  loop {l}: {n} curves at distance {expected:.9}");
        for (c, curve) in curves.iter().enumerate() {
            let (_, end) = curve_endpoints(curve);
            let (next_start, _) = curve_endpoints(&curves[(c + 1) % n]);
            if !same_point(end, next_start) {
                errors.push(format!(
                    "loop {l} curve {c}: end {end:?} != next start {next_start:?}, {:e} away",
                    dist(end, next_start)
                ));
            }
            if let OffsetCurve::Arc {
                start, end, center, ..
            } = *curve
            {
                let (rs, re) = (dist(start, center), dist(end, center));
                if (rs - re).abs() > tolerance {
                    errors.push(format!("loop {l} curve {c}: arc radii differ {rs} vs {re}"));
                }
            }
            for p in curve_samples(curve) {
                let d = nearest_input_distance(inputs, p);
                if (d - expected).abs() > tolerance {
                    errors.push(format!("loop {l} curve {c} {curve:?}: sample {p:?} at distance {d:.9}, expected {expected:.9}"));
                }
            }
        }
    }
    for e in &errors {
        println!("  {e}");
    }
    assert!(
        errors.is_empty(),
        "{name}: {} geometry errors",
        errors.len()
    );
    loop_dists
}

fn polygon(points: &[(f64, f64)]) -> Vec<Input> {
    (0..points.len())
        .map(|i| Input::Line {
            start: points[i],
            end: points[(i + 1) % points.len()],
        })
        .collect()
}

fn params(radius: f64, left: bool, right: bool) -> OffsetParameters {
    OffsetParameters {
        offset_radius: radius,
        left_offset: left,
        right_offset: right,
        ..OffsetParameters::default()
    }
}

const SQUARE: [(f64, f64); 4] = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
const L_SHAPE: [(f64, f64); 6] = [
    (0.0, 0.0),
    (10.0, 0.0),
    (10.0, 4.0),
    (4.0, 4.0),
    (4.0, 10.0),
    (0.0, 10.0),
];

fn star(spikes: usize, outer: f64, inner: f64) -> Vec<(f64, f64)> {
    (0..2 * spikes)
        .map(|k| {
            let r = if k % 2 == 0 { outer } else { inner };
            let a = PI * k as f64 / spikes as f64 + 0.1;
            (r * a.cos(), r * a.sin())
        })
        .collect()
}

fn slot() -> Vec<Input> {
    vec![
        Input::Line {
            start: (0.0, 0.0),
            end: (10.0, 0.0),
        },
        Input::Arc {
            start: (10.0, 0.0),
            end: (10.0, 4.0),
            center: (10.0, 2.0),
            orientation: Orientation::CounterClockwise,
        },
        Input::Line {
            start: (10.0, 4.0),
            end: (0.0, 4.0),
        },
        Input::Arc {
            start: (0.0, 4.0),
            end: (0.0, 0.0),
            center: (0.0, 2.0),
            orientation: Orientation::CounterClockwise,
        },
    ]
}

fn circle(r: f64, orientation: Orientation) -> Vec<Input> {
    let (a, b) = ((r, 0.0), (-r, 0.0));
    vec![
        Input::Arc {
            start: a,
            end: b,
            center: (0.0, 0.0),
            orientation,
        },
        Input::Arc {
            start: b,
            end: a,
            center: (0.0, 0.0),
            orientation,
        },
    ]
}

const CCW: Orientation = Orientation::CounterClockwise;
const CW: Orientation = Orientation::Clockwise;

#[test]
fn square_outer_offset_is_rounded_square() {
    let dists = check_offsets(
        "square outer",
        &polygon(&SQUARE),
        params(1.0, false, true),
        LINE_TOLERANCE,
    );
    assert_eq!(dists.len(), 1);
    let v = run(&polygon(&SQUARE), params(1.0, false, true));
    let curves: Vec<_> = (0..v.num_offset_curves(0).unwrap())
        .map(|c| v.curve(0, c).unwrap())
        .collect();
    let arcs = curves
        .iter()
        .filter(|c| matches!(c, OffsetCurve::Arc { .. }))
        .count();
    let lines = curves.len() - arcs;
    assert_eq!(
        (lines, arcs),
        (4, 4),
        "expected 4 straight sides and 4 corner arcs, got {curves:#?}"
    );
}

#[test]
fn square_inner_offset_is_smaller_square() {
    let dists = check_offsets(
        "square inner",
        &polygon(&SQUARE),
        params(1.0, true, false),
        LINE_TOLERANCE,
    );
    assert_eq!(dists.len(), 1);
    let v = run(&polygon(&SQUARE), params(1.0, true, false));
    let corners: Vec<_> = (0..4)
        .map(|c| curve_endpoints(&v.curve(0, c).unwrap()).0)
        .collect();
    for expected in [(1.0, 1.0), (9.0, 1.0), (9.0, 9.0), (1.0, 9.0)] {
        assert!(
            corners.iter().any(|&p| dist(p, expected) < LINE_TOLERANCE),
            "corner {expected:?} missing from {corners:?}"
        );
    }
}

#[test]
fn both_sides_selected_offsets_both_sides() {
    let dists = check_offsets(
        "square both",
        &polygon(&SQUARE),
        params(1.0, true, true),
        LINE_TOLERANCE,
    );
    assert_eq!(dists.len(), 2, "inner and outer loop expected");
}

#[test]
fn l_shape_offsets_round_only_the_reflex_corner_inside() {
    check_offsets(
        "L inner",
        &polygon(&L_SHAPE),
        params(1.0, true, false),
        LINE_TOLERANCE,
    );
    check_offsets(
        "L outer",
        &polygon(&L_SHAPE),
        params(1.5, false, true),
        LINE_TOLERANCE,
    );
}

#[test]
fn star_offsets_both_sides_at_several_radii() {
    let s = polygon(&star(7, 10.0, 4.0));
    for r in [0.25, 0.8, 2.0] {
        check_offsets(
            &format!("star r={r}"),
            &s,
            params(r, true, true),
            LINE_TOLERANCE,
        );
    }
}

#[test]
fn inner_offset_vanishes_beyond_inradius() {
    let v = run(&polygon(&SQUARE), params(6.0, true, false));
    assert_eq!(
        v.num_offset_loops(),
        0,
        "no inner offset exists at distance 6 inside a 10x10 square"
    );
}

#[test]
fn depletion_produces_successive_offsets() {
    let p = OffsetParameters {
        depletion: 1.0,
        ..params(1.0, true, false)
    };
    let dists = check_offsets("square inner stepped", &polygon(&SQUARE), p, LINE_TOLERANCE);
    let expected = [1.0, 2.0, 3.0, 4.0];
    assert_eq!(dists.len(), expected.len(), "distances {dists:?}");
    for (d, e) in dists.iter().zip(expected) {
        assert!(
            (d - e).abs() < LINE_TOLERANCE,
            "distances {dists:?}, expected {expected:?}"
        );
    }
}

#[test]
fn open_chain_offsets_around_both_sides() {
    let chain = vec![
        Input::Line {
            start: (0.0, 0.0),
            end: (5.0, 0.0),
        },
        Input::Line {
            start: (5.0, 0.0),
            end: (5.0, 5.0),
        },
    ];
    let dists = check_offsets(
        "open chain",
        &chain,
        params(1.0, true, true),
        LINE_TOLERANCE,
    );
    assert_eq!(
        dists.len(),
        1,
        "an open chain has a single offset loop around it"
    );
}

#[test]
fn circle_offsets_are_concentric_circles() {
    let apx = OffsetParameters::default().apx_arcs;
    for (orientation, left, expected_r) in [
        (CCW, false, 6.0),
        (CCW, true, 4.0),
        (CW, true, 6.0),
        (CW, false, 4.0),
    ] {
        let p = params(1.0, left, !left);
        check_offsets(
            &format!("circle {orientation:?} left={left}"),
            &circle(5.0, orientation),
            p,
            2.0 * apx,
        );
        let v = run(&circle(5.0, orientation), p);
        assert_eq!(v.num_offset_loops(), 1);
        for c in 0..v.num_offset_curves(0).unwrap() {
            if let OffsetCurve::Arc {
                start, end, center, ..
            } = v.curve(0, c).unwrap()
                && dist(center, (0.0, 0.0)) < LINE_TOLERANCE
            {
                assert!(
                    (dist(start, center) - expected_r).abs() < 2.0 * apx,
                    "radius {}",
                    dist(start, center)
                );
                assert!((dist(end, center) - expected_r).abs() < 2.0 * apx);
            }
        }
    }
}

#[test]
fn slot_offsets_follow_lines_and_arcs() {
    let apx = OffsetParameters::default().apx_arcs;
    check_offsets("slot outer", &slot(), params(1.0, false, true), 2.0 * apx);
    check_offsets("slot inner", &slot(), params(1.0, true, false), 2.0 * apx);
}

#[test]
fn unrecovered_arcs_are_offset_as_their_polygonal_approximation() {
    let p = OffsetParameters {
        recover_arcs: false,
        ..params(1.0, false, true)
    };
    let apx = p.apx_arcs;
    check_offsets("circle no recovery", &circle(5.0, CCW), p, 2.0 * apx);
    let v = run(&circle(5.0, CCW), p);
    let n = v.num_offset_curves(0).unwrap();
    assert!(
        n > 8,
        "without recovery the offset follows the approximating polygon, got {n} curves"
    );
}

#[test]
fn arcs_without_approximation_are_rejected_and_input_is_kept() {
    let mut v = Offsetter::new();
    v.parameters.apx_arcs = 0.0;
    for input in circle(5.0, CCW) {
        if let Input::Arc {
            start,
            end,
            center,
            orientation,
        } = input
        {
            v.add_arc(start, end, center, orientation);
        }
    }
    assert_eq!(v.compute_offset(), Err(OffsetError::ArcsNotApproximated));
    assert_eq!(
        v.num_offset_loops(),
        0,
        "no offsets are reported after an error"
    );
    v.parameters.apx_arcs = OffsetParameters::default().apx_arcs;
    v.compute_offset()
        .expect("the same input computes once arcs may be approximated");
    assert_eq!(
        v.num_offset_loops(),
        1,
        "circle input was kept across the error"
    );
}

#[test]
fn reset_input_forgets_the_snapping_end_point() {
    let mut v = Offsetter::new();
    v.parameters.tolerance = 0.01;
    v.add_line((0.0, 0.0), (10.0, 0.0));
    v.reset_input();
    assert_eq!(
        v.add_line((10.005, 0.0), (10.0, 10.0)),
        None,
        "discarded geometry must not attract the next start point"
    );
}

#[test]
fn tangent_joints_keep_their_vertex_arcs() {
    let joints = [(5.0, 0.0), (-5.0, 0.0)];
    for arithmetic in ARITHMETICS {
        let v = run(
            &circle(5.0, CCW),
            OffsetParameters {
                arithmetic,
                ..params(1.0, false, true)
            },
        );
        let curves: Vec<_> = (0..v.num_offset_curves(0).unwrap())
            .map(|c| v.curve(0, c).unwrap())
            .collect();
        println!("[{arithmetic:?}] {curves:#?}");
        assert_eq!(
            curves.len(),
            4,
            "[{arithmetic:?}] the two recovered arcs and a vertex arc at each joint"
        );
        let vertex_arcs: Vec<_> = curves
            .iter()
            .filter_map(|curve| match *curve {
                OffsetCurve::Arc {
                    start,
                    end,
                    center,
                    orientation: CW,
                } if joints.contains(&center) => Some(dist(start, end)),
                _ => None,
            })
            .collect();
        assert_eq!(
            vertex_arcs.len(),
            joints.len(),
            "[{arithmetic:?}] one clockwise arc about each joint"
        );
        assert!(
            vertex_arcs.iter().all(|&length| length <= LINE_TOLERANCE),
            "[{arithmetic:?}] the vertex arcs have no length: {vertex_arcs:?}"
        );
    }
}

/// A stored FomCam contour (`BO4465`, machining "Arslan Samlple"): four arcs, counter-clockwise,
/// meeting tangentially; each row is the bits of its start, end and centre.
const TANGENT_ARCS_BITS: [[u64; 6]; 4] = [
    [
        0x405226666666665f,
        0x400bfffffffffc38,
        0xc05226666666663c,
        0x400c000000000612,
        0xbd7af0f885e98a7e,
        0xc06c6905bd6ac337,
    ],
    [
        0xc05226666666663c,
        0x400c000000000612,
        0xc0522666666666a6,
        0xc00bfffffffffda0,
        0xc051dfeef8d02013,
        0x3d511bbe12287f6c,
    ],
    [
        0xc0522666666666a6,
        0xc00bfffffffffda0,
        0x40522666666665fd,
        0xc00c000000000439,
        0xbd037dae9dda1edd,
        0x406c6905bd6ac330,
    ],
    [
        0x40522666666665fd,
        0xc00c000000000439,
        0x405226666666665f,
        0x400bfffffffffc38,
        0x4051dfeef8d01f64,
        0xbd5099e0072083e7,
    ],
];
const TANGENT_ARCS_RADIUS: f64 = 2.5;
const TANGENT_ARCS_APX: f64 = 0.1;
/// How far a vertex arc's centre may lie from its input vertex after scaling and unscaling.
const VERTEX_CENTRE_TOLERANCE: f64 = 1e-12;

fn tangent_arcs() -> Vec<Input> {
    let point = |x: u64, y: u64| (f64::from_bits(x), f64::from_bits(y));
    TANGENT_ARCS_BITS
        .iter()
        .map(|&[xs, ys, xe, ye, xc, yc]| Input::Arc {
            start: point(xs, ys),
            end: point(xe, ye),
            center: point(xc, yc),
            orientation: CCW,
        })
        .collect()
}

#[test]
fn curves_shorter_than_zero_keep_the_loop_closed_bit_for_bit() {
    let inputs = tangent_arcs();
    let vertices: Vec<(f64, f64)> = inputs
        .iter()
        .map(|input| match *input {
            Input::Arc { start, .. } | Input::Line { start, .. } => start,
        })
        .collect();
    for arithmetic in ARITHMETICS {
        let v = run(
            &inputs,
            OffsetParameters {
                left_offset: true,
                right_offset: false,
                offset_radius: TANGENT_ARCS_RADIUS,
                recover_arcs: true,
                apx_arcs: TANGENT_ARCS_APX,
                depletion: 0.0,
                tolerance: 0.0,
                arithmetic,
            },
        );
        assert_eq!(v.num_offset_loops(), 1, "[{arithmetic:?}] one loop");
        let n = v.num_offset_curves(0).expect("loop exists");
        let curves: Vec<OffsetCurve> = (0..n)
            .map(|c| v.curve(0, c).expect("curve exists"))
            .collect();
        let mut gaps = Vec::new();
        for (c, curve) in curves.iter().enumerate() {
            let (start, end) = curve_endpoints(curve);
            let (next_start, _) = curve_endpoints(&curves[(c + 1) % n]);
            println!(
                "[{arithmetic:?}] curve {c}: {curve:?}, length {:e}",
                dist(start, end)
            );
            if !same_point(end, next_start) {
                gaps.push(format!(
                    "curve {c} ends at {end:?}, curve {} starts at {next_start:?}, {:e} away",
                    (c + 1) % n,
                    dist(end, next_start)
                ));
            }
        }
        for gap in &gaps {
            println!("[{arithmetic:?}] gap: {gap}");
        }
        assert!(
            gaps.is_empty(),
            "[{arithmetic:?}] {} gaps in the loop",
            gaps.len()
        );
        let vertex_arcs: Vec<(f64, f64)> = curves
            .iter()
            .filter_map(|curve| match *curve {
                OffsetCurve::Arc {
                    center,
                    orientation: CW,
                    ..
                } => Some(center),
                _ => None,
            })
            .collect();
        assert_eq!(
            vertex_arcs.len(),
            vertices.len(),
            "[{arithmetic:?}] one clockwise vertex arc at each tangent joint: {vertex_arcs:?}"
        );
        for centre in &vertex_arcs {
            assert!(
                vertices
                    .iter()
                    .any(|&vertex| dist(*centre, vertex) <= VERTEX_CENTRE_TOLERANCE),
                "[{arithmetic:?}] vertex arc about {centre:?} lies at no input vertex {vertices:?}"
            );
        }
    }
}

#[test]
fn start_points_within_tolerance_snap_to_previous_end() {
    let mut v = Offsetter::new();
    v.parameters.tolerance = 0.01;
    assert_eq!(
        v.add_line((0.0, 0.0), (10.0, 0.0)),
        None,
        "first line has no previous end"
    );
    assert_eq!(
        v.add_line((10.0, 0.0), (10.0, 10.0)),
        None,
        "exact continuation is not a snap"
    );
    let gap = v
        .add_line((10.005, 10.0), (0.0, 10.0))
        .expect("gap within tolerance is snapped");
    assert!((gap - 0.005).abs() < 1e-12, "gap {gap}");
    v.add_line((0.0, 10.0), (0.0, 0.0));
    v.parameters.offset_radius = 1.0;
    v.compute_offset().expect("compute");
    assert_eq!(
        v.num_offset_loops(),
        1,
        "snapped input forms a closed square"
    );
    let corners: Vec<_> = (0..v.num_offset_curves(0).unwrap())
        .filter_map(|c| match v.curve(0, c).unwrap() {
            OffsetCurve::Arc { center, .. } => Some(center),
            OffsetCurve::Line { .. } => None,
        })
        .collect();
    assert!(
        corners
            .iter()
            .any(|&c| dist(c, (10.0, 10.0)) < LINE_TOLERANCE),
        "snapped corner (10,10) expected among {corners:?}"
    );

    let mut v = Offsetter::new();
    v.parameters.tolerance = 0.01;
    v.add_line((0.0, 0.0), (10.0, 0.0));
    assert_eq!(
        v.add_line((10.5, 0.0), (10.0, 10.0)),
        None,
        "gap beyond tolerance is kept"
    );
}

#[test]
fn out_of_range_queries_return_none() {
    let v = run(&polygon(&SQUARE), params(1.0, false, true));
    assert_eq!(v.num_offset_curves(1), None);
    let n = v.num_offset_curves(0).unwrap();
    assert!(v.curve(0, n).is_none());
    assert!(v.curve(1, 0).is_none());
}

#[test]
fn compute_offset_can_be_repeated_with_new_parameters() {
    let mut v = run(&polygon(&SQUARE), params(1.0, false, true));
    assert_eq!(v.num_offset_loops(), 1);
    v.parameters.offset_radius = 0.0;
    v.compute_offset()
        .expect("non-positive radius computes no offset");
    assert_eq!(
        v.num_offset_loops(),
        0,
        "stale loops from the previous run must not remain"
    );
    v.parameters = params(2.0, true, false);
    v.compute_offset().expect("recompute");
    let n = v.num_offset_curves(0).unwrap();
    let corners: Vec<_> = (0..n)
        .map(|c| curve_endpoints(&v.curve(0, c).unwrap()).0)
        .collect();
    assert!(
        corners
            .iter()
            .any(|&p| dist(p, (2.0, 2.0)) < LINE_TOLERANCE),
        "inner offset at 2 expected, got {corners:?}"
    );
}

#[test]
fn neither_side_selected_offsets_both_sides() {
    let dists = check_offsets(
        "square neither",
        &polygon(&SQUARE),
        params(1.0, false, false),
        LINE_TOLERANCE,
    );
    assert_eq!(dists.len(), 2, "inner and outer loop expected");
}

/// Every curve of every loop as the bit patterns of its coordinates.
fn curve_bits(v: &Offsetter) -> Vec<Vec<Vec<u64>>> {
    (0..v.num_offset_loops())
        .map(|l| {
            (0..v.num_offset_curves(l).expect("loop exists"))
                .map(|c| match v.curve(l, c).expect("curve exists") {
                    OffsetCurve::Line { start, end } => [start.0, start.1, end.0, end.1]
                        .iter()
                        .map(|x| x.to_bits())
                        .collect(),
                    OffsetCurve::Arc {
                        start,
                        end,
                        center,
                        orientation,
                    } => [start.0, start.1, end.0, end.1, center.0, center.1]
                        .iter()
                        .map(|x| x.to_bits())
                        .chain([u64::from(orientation == CW)])
                        .collect(),
                })
                .collect()
        })
        .collect()
}

/// A star from the stress generator whose offset at 0.5 has different bits in the two
/// arithmetics.
const ARITHMETIC_SENSITIVE_STAR: [(f64, f64); 5] = [
    (4.426874060425207, 1.7480605463141181),
    (0.3477302992469246, 3.5066862468025337),
    (-7.153602676193704, 2.5183478991067965),
    (-3.806764447049871, -4.095410996898055),
    (3.8196208126703555, -8.637793073017678),
];
const ARITHMETIC_SENSITIVE_RADIUS: f64 = 0.5;

#[test]
fn extended_differs_from_double() {
    let star = polygon(&ARITHMETIC_SENSITIVE_STAR);
    let base = params(ARITHMETIC_SENSITIVE_RADIUS, false, true);
    let [double, extended] = ARITHMETICS
        .map(|arithmetic| curve_bits(&run(&star, OffsetParameters { arithmetic, ..base })));
    let differing: usize = double
        .iter()
        .zip(&extended)
        .map(|(d, e)| d.iter().zip(e).filter(|(x, y)| x != y).count())
        .sum();
    println!(
        "star at {ARITHMETIC_SENSITIVE_RADIUS}: {differing} curves differ between the arithmetics"
    );
    assert!(differing > 0, "Extended gave the same bits as Double");
}

#[test]
fn switching_arithmetic_between_calls() {
    let star = polygon(&ARITHMETIC_SENSITIVE_STAR);
    let base = params(ARITHMETIC_SENSITIVE_RADIUS, false, true);
    let mut reused = run(&star, base);
    for arithmetic in [Arithmetic::Double, Arithmetic::Extended, Arithmetic::Double] {
        let p = OffsetParameters { arithmetic, ..base };
        reused.parameters = p;
        reused.compute_offset().expect("compute_offset");
        let fresh = run(&star, p);
        let (again, alone) = (curve_bits(&reused), curve_bits(&fresh));
        for (l, (a, b)) in again.iter().zip(&alone).enumerate() {
            for (c, (x, y)) in a.iter().zip(b).enumerate() {
                if x != y {
                    println!("[{arithmetic:?}] loop {l} curve {c}: reused {x:x?}, fresh {y:x?}");
                }
            }
        }
        println!("[{arithmetic:?}]: {} loops", again.len());
        assert_eq!(
            again, alone,
            "[{arithmetic:?}] reused Offsetter differs from a fresh one"
        );
    }
}
