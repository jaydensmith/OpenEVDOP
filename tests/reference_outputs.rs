//! Replays tests/reference/*.script and compares the offsets with the recorded reference output
//! (*.expected).
//!
//! Script lines: `P left right radius recover_arcs apx_arcs depletion tolerance`,
//! `L xs ys xe ye`, `A xs ys xe ye xc yc clockwise`, `R` (reset input), `C` (compute and dump).
use openevdop::{OffsetCurve, OffsetParameters, Offsetter, Orientation};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

const REFERENCE_DIR: &str = "tests/reference";
/// Values agree bit for bit except vertices of approximated arcs, which come from sin/cos and may
/// differ in the last bit between math libraries.
const RELATIVE_TOLERANCE: f64 = 1e-12;

fn parse_bool(token: &str) -> bool {
    token != "0"
}

fn numbers(fields: &[&str]) -> Vec<f64> {
    fields
        .iter()
        .map(|f| f.parse().unwrap_or_else(|_| panic!("bad number {f}")))
        .collect()
}

fn dump(v: &Offsetter, out: &mut String) {
    writeln!(out, "loops {}", v.num_offset_loops()).unwrap();
    for l in 0..v.num_offset_loops() {
        let n = v.num_offset_curves(l).unwrap();
        writeln!(out, "loop {l} curves {n}").unwrap();
        for c in 0..n {
            match v.curve(l, c).unwrap() {
                OffsetCurve::Line { start, end } => {
                    writeln!(out, "line {} {} {} {}", start.0, start.1, end.0, end.1).unwrap()
                }
                OffsetCurve::Arc {
                    start,
                    end,
                    center,
                    orientation,
                } => writeln!(
                    out,
                    "arc {} {} {} {} {} {} {}",
                    start.0,
                    start.1,
                    end.0,
                    end.1,
                    center.0,
                    center.1,
                    u8::from(orientation == Orientation::Clockwise)
                )
                .unwrap(),
            }
        }
    }
}

fn replay(script: &str) -> String {
    let mut v = Offsetter::new();
    let mut out = String::new();
    for line in script.lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        match f.first().copied() {
            Some("P") => {
                let n = numbers(&f[1..]);
                v.parameters = OffsetParameters {
                    left_offset: parse_bool(f[1]),
                    right_offset: parse_bool(f[2]),
                    offset_radius: n[2],
                    recover_arcs: parse_bool(f[4]),
                    apx_arcs: n[4],
                    depletion: n[5],
                    tolerance: n[6],
                };
            }
            Some("L") => {
                let n = numbers(&f[1..]);
                if let Some(gap) = v.add_line((n[0], n[1]), (n[2], n[3])) {
                    writeln!(out, "snap line {gap:.6}").unwrap();
                }
            }
            Some("A") => {
                let n = numbers(&f[1..7]);
                let orientation = if parse_bool(f[7]) {
                    Orientation::Clockwise
                } else {
                    Orientation::CounterClockwise
                };
                if let Some(gap) = v.add_arc((n[0], n[1]), (n[2], n[3]), (n[4], n[5]), orientation)
                {
                    writeln!(out, "snap arc {gap:.6}").unwrap();
                }
            }
            Some("R") => v.reset_input(),
            Some("C") => {
                if v.compute_offset().is_err() {
                    writeln!(out, "error").unwrap();
                }
                writeln!(out, "compute").unwrap();
                dump(&v, &mut out);
            }
            _ => {}
        }
    }
    out
}

fn tokens_match(expected: &str, actual: &str) -> bool {
    match (expected.parse::<f64>(), actual.parse::<f64>()) {
        (Ok(e), Ok(a)) => (e - a).abs() <= RELATIVE_TOLERANCE * e.abs().max(1.0),
        _ => expected == actual,
    }
}

/// Line-by-line comparison; returns human-readable mismatches.
fn compare(expected: &str, actual: &str) -> Vec<String> {
    let (e_lines, a_lines): (Vec<&str>, Vec<&str>) =
        (expected.lines().collect(), actual.lines().collect());
    let mut diffs = Vec::new();
    for i in 0..e_lines.len().max(a_lines.len()) {
        let (e, a) = (
            e_lines.get(i).copied().unwrap_or("<missing>"),
            a_lines.get(i).copied().unwrap_or("<missing>"),
        );
        let (et, at): (Vec<&str>, Vec<&str>) = (
            e.split_whitespace().collect(),
            a.split_whitespace().collect(),
        );
        if et.len() != at.len() || !et.iter().zip(&at).all(|(x, y)| tokens_match(x, y)) {
            diffs.push(format!(
                "line {}:\n    expected: {e}\n    actual:   {a}",
                i + 1
            ));
        }
    }
    diffs
}

#[test]
fn matches_reference_outputs() {
    let mut scripts: Vec<_> = fs::read_dir(REFERENCE_DIR)
        .expect("tests/reference exists")
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "script"))
        .collect();
    scripts.sort();
    assert!(!scripts.is_empty(), "no reference scripts found");
    let mut failed = Vec::new();
    for script in &scripts {
        let name = script.file_stem().unwrap().to_string_lossy().into_owned();
        let reference = Path::new(REFERENCE_DIR).join(format!("{name}.expected"));
        let expected = fs::read_to_string(&reference)
            .unwrap_or_else(|_| panic!("missing {}", reference.display()));
        let actual = replay(&fs::read_to_string(script).unwrap());
        let diffs = compare(&expected, &actual);
        println!(
            "{name}: {} lines, {} mismatches",
            expected.lines().count(),
            diffs.len()
        );
        for d in diffs.iter().take(10) {
            println!("  {d}");
        }
        if !diffs.is_empty() {
            failed.push(name);
        }
    }
    assert!(
        failed.is_empty(),
        "scripts differing from the reference: {failed:?}"
    );
}
