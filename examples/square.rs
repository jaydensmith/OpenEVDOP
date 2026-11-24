use open_evdop::{OffsetCurve, Offsetter};

fn main() {
    let mut v = Offsetter::new();
    v.parameters.offset_radius = 1.0;
    let sq = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
    for i in 0..4 {
        v.add_line(sq[i], sq[(i + 1) % 4]);
    }
    let r = v.compute_offset();
    println!("result {r:?} loops {}", v.num_offset_loops());
    for l in 0..v.num_offset_loops() {
        let n = v.num_offset_curves(l).unwrap();
        println!("loop {l}: {n} curves");
        for c in 0..n {
            match v.curve(l, c).unwrap() {
                OffsetCurve::Line { start, end } => println!("  line {start:?} -> {end:?}"),
                OffsetCurve::Arc {
                    start,
                    end,
                    center,
                    orientation,
                } => println!("  arc  {start:?} -> {end:?} c={center:?} {orientation:?}"),
            }
        }
    }
}
