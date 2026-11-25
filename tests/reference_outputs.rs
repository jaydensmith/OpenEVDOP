//! Replays tests/reference/*.script and compares the offsets with the recorded reference output:
//! *.expected for `Arithmetic::Double`, *.extended.expected for `Arithmetic::Extended`.
//!
//! Script lines: `P left right radius recover_arcs apx_arcs depletion tolerance`,
//! `L xs ys xe ye`, `A xs ys xe ye xc yc clockwise`, `R` (reset input), `C` (compute and dump).
use open_evdop::{Arithmetic, OffsetCurve, OffsetParameters, Offsetter, Orientation};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

const REFERENCE_DIR: &str = "tests/reference";
const ARITHMETICS: [Arithmetic; 2] = [Arithmetic::Double, Arithmetic::Extended];

fn reference_path(name: &str, arithmetic: Arithmetic) -> std::path::PathBuf {
    let suffix = match arithmetic {
        Arithmetic::Double => "expected",
        Arithmetic::Extended => "extended.expected",
    };
    Path::new(REFERENCE_DIR).join(format!("{name}.{suffix}"))
}

fn script_names() -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(REFERENCE_DIR)
        .expect("tests/reference exists")
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "script"))
        .map(|p| p.file_stem().unwrap().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert!(!names.is_empty(), "no reference scripts found");
    names
}

fn read_reference(name: &str, arithmetic: Arithmetic) -> String {
    let path = reference_path(name, arithmetic);
    fs::read_to_string(&path).unwrap_or_else(|_| panic!("missing {}", path.display()))
}

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

fn replay(script: &str, arithmetic: Arithmetic) -> String {
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
                    arithmetic,
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

/// Numbers agree when they are the same double; the files print some with 17 significant digits.
fn tokens_match(expected: &str, actual: &str) -> bool {
    match (expected.parse::<f64>(), actual.parse::<f64>()) {
        (Ok(e), Ok(a)) => e.to_bits() == a.to_bits(),
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
    let mut failed = Vec::new();
    for arithmetic in ARITHMETICS {
        for name in script_names() {
            let expected = read_reference(&name, arithmetic);
            let script = Path::new(REFERENCE_DIR).join(format!("{name}.script"));
            let actual = replay(&fs::read_to_string(script).unwrap(), arithmetic);
            let diffs = compare(&expected, &actual);
            println!(
                "{name} ({arithmetic:?}): {} lines, {} mismatches",
                expected.lines().count(),
                diffs.len()
            );
            for d in diffs.iter().take(10) {
                println!("  {d}");
            }
            if !diffs.is_empty() {
                failed.push(format!("{name} ({arithmetic:?})"));
            }
        }
    }
    assert!(
        failed.is_empty(),
        "scripts differing from the reference: {failed:?}"
    );
}

#[test]
fn extended_differs_from_double() {
    let mut differing = Vec::new();
    for name in script_names() {
        let diffs = compare(
            &read_reference(&name, Arithmetic::Double),
            &read_reference(&name, Arithmetic::Extended),
        );
        println!("{name}: {} lines differ between the modes", diffs.len());
        for d in diffs.iter().take(3) {
            println!("  {d}");
        }
        if !diffs.is_empty() {
            differing.push(name);
        }
    }
    assert!(
        !differing.is_empty(),
        "no reference script tells the arithmetic modes apart"
    );
}
