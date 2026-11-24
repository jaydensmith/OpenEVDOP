use openevdop::{OffsetCurve, Offsetter, Orientation};
use std::f64::consts::TAU;

const TOLERANCE: f64 = 1e-7;

struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A simple polygon: vertices at strictly increasing angles with random radii.
fn random_star(rng: &mut Lcg, n: usize, r_min: f64, r_max: f64) -> Vec<(f64, f64)> {
    (0..n)
        .map(|k| {
            let a = TAU * (k as f64 + 0.5 * rng.next_f64()) / n as f64;
            let r = r_min + (r_max - r_min) * rng.next_f64();
            (r * a.cos(), r * a.sin())
        })
        .collect()
}

fn seg_dist(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
    ((a.0 + t * dx - p.0).powi(2) + (a.1 + t * dy - p.1).powi(2)).sqrt()
}

fn sample(curve: OffsetCurve, t: f64) -> (f64, f64) {
    match curve {
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
            let a0 = (start.1 - center.1).atan2(start.0 - center.0);
            let a1 = (end.1 - center.1).atan2(end.0 - center.0);
            let ccw = (a1 - a0).rem_euclid(TAU);
            let sweep = match orientation {
                Orientation::CounterClockwise => ccw,
                Orientation::Clockwise => ccw - TAU,
            };
            let r = ((start.0 - center.0).powi(2) + (start.1 - center.1).powi(2)).sqrt();
            let a = a0 + t * sweep;
            (center.0 + r * a.cos(), center.1 + r * a.sin())
        }
    }
}

fn check_polygon(label: &str, poly: &[(f64, f64)], radius: f64) -> usize {
    let mut v = Offsetter::new();
    v.parameters.offset_radius = radius;
    v.parameters.left_offset = true;
    v.parameters.right_offset = true;
    for i in 0..poly.len() {
        v.add_line(poly[i], poly[(i + 1) % poly.len()]);
    }
    v.compute_offset()
        .unwrap_or_else(|e| panic!("{label}: {e}"));
    let mut worst = 0.0_f64;
    let mut curves = 0;
    for l in 0..v.num_offset_loops() {
        for c in 0..v.num_offset_curves(l).unwrap() {
            curves += 1;
            let curve = v.curve(l, c).unwrap();
            for k in 0..=4 {
                let p = sample(curve, k as f64 / 4.0);
                let d = (0..poly.len())
                    .map(|i| seg_dist(p, poly[i], poly[(i + 1) % poly.len()]))
                    .fold(f64::INFINITY, f64::min);
                worst = worst.max((d - radius).abs());
            }
        }
    }
    println!(
        "{label}: {} loops, {curves} curves, worst distance error {worst:.3e}",
        v.num_offset_loops()
    );
    assert!(
        worst <= TOLERANCE * radius.max(1.0),
        "{label}: worst distance error {worst:e}"
    );
    assert!(
        v.num_offset_loops() >= 1,
        "{label}: outer offset always exists"
    );
    curves
}

#[test]
fn random_star_polygons_offset_at_exact_distance() {
    let mut rng = Lcg(0x5eed);
    for (n, trials) in [(5, 40), (50, 20), (500, 4)] {
        for trial in 0..trials {
            let poly = random_star(&mut rng, n, 3.0, 10.0);
            for radius in [0.05, 0.5, 2.5] {
                check_polygon(&format!("n={n} trial={trial} r={radius}"), &poly, radius);
            }
        }
    }
}

#[test]
fn large_polygon_completes() {
    let mut rng = Lcg(42);
    let poly = random_star(&mut rng, 5000, 50.0, 100.0);
    let curves = check_polygon("n=5000", &poly, 1.0);
    assert!(curves > 0);
}

#[test]
fn regular_polygon_with_collinear_and_cocircular_vertices() {
    let n = 64;
    let poly: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            (
                10.0 * (TAU * k as f64 / n as f64).cos(),
                10.0 * (TAU * k as f64 / n as f64).sin(),
            )
        })
        .collect();
    check_polygon("regular 64-gon", &poly, 1.0);
    let mut grid_square = Vec::new();
    for k in 0..10 {
        grid_square.push((k as f64, 0.0));
    }
    for k in 0..10 {
        grid_square.push((10.0, k as f64));
    }
    for k in 0..10 {
        grid_square.push((10.0 - k as f64, 10.0));
    }
    for k in 0..10 {
        grid_square.push((0.0, 10.0 - k as f64));
    }
    check_polygon("square with collinear vertices", &grid_square, 1.0);
}
