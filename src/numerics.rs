use crate::consts::*;
use crate::real::Real;
use crate::state::*;

const NUM_SAMPLES: usize = 11;
const MAX_EXPANSIONS: usize = 50;

/// Gaussian elimination with partial pivoting.
/// Returns the solution and the eliminated pivot; every singular or
/// underdetermined system yields `None`.
pub(crate) fn solve_2x2<R: Real>(a: [[f64; 2]; 2], b: [f64; 2]) -> Option<([f64; 2], f64)> {
    let nonzero = |v: f64| v > PIVOT_TINY || v < -PIVOT_TINY;
    let (col, other) = if nonzero(a[0][0]) || nonzero(a[1][0]) {
        (0, 1)
    } else if nonzero(a[0][1]) || nonzero(a[1][1]) {
        (1, 0)
    } else {
        return None;
    };
    let (p, q) = if a[0][col].abs() <= a[1][col].abs() {
        (1, 0)
    } else {
        (0, 1)
    };
    let ratio = R::of(a[q][col]) / a[p][col];
    let elim = R::of(a[q][other]) - R::of(a[p][other]) * ratio;
    if !(elim > PIVOT_TINY || elim < -PIVOT_TINY) {
        return None;
    }
    let other_value = (R::of(b[q]) - R::of(b[p]) * ratio) / elim;
    let mut x = [0.0; 2];
    x[other] = other_value.stored();
    x[col] = ((R::of(b[p]) - other_value * a[p][other]) / a[p][col]).stored();
    Some((x, elim.stored()))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Qualification {
    First,
    Second,
    Both,
    Neither,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum QuadRoots {
    NoRoot,
    One(f64),
    Two(f64, f64),
}

/// Roots of `a t^2 + b t + c = 0` with coefficient rescaling and
/// numerically stable root selection.
pub(crate) fn solve_quadratic<R: Real>(mut a: f64, mut b: f64, mut c: f64) -> QuadRoots {
    if a <= PARALLEL_EPS {
        loop {
            if a < -PARALLEL_EPS
                || b > PARALLEL_EPS
                || b < -PARALLEL_EPS
                || c > PARALLEL_EPS
                || c < -PARALLEL_EPS
                || (a == 0.0 && b == 0.0 && c == 0.0)
            {
                break;
            }
            a = (R::of(a) * RESCALE).stored();
            b = (R::of(b) * RESCALE).stored();
            c = (R::of(c) * RESCALE).stored();
            if a > PARALLEL_EPS {
                break;
            }
        }
    }
    if a > PIVOT_TINY || a < -PIVOT_TINY {
        let disc = (R::of(b) * b - R::of(a) * 4.0 * c).stored();
        if disc <= 0.0 {
            if disc < -DISC_ZERO {
                return QuadRoots::NoRoot;
            }
            return QuadRoots::One((-(R::of(b) / (R::of(a) * 2.0))).stored());
        }
        let root = R::of(disc).sqrt().stored();
        let qq = (if b <= 0.0 {
            R::of(b) - root
        } else {
            R::of(root) + b
        } * -0.5)
            .stored();
        if qq > PIVOT_TINY || qq < -PIVOT_TINY {
            QuadRoots::Two((R::of(qq) / a).stored(), (R::of(c) / qq).stored())
        } else {
            QuadRoots::One((R::of(qq) / a).stored())
        }
    } else if b > PIVOT_TINY || b < -PIVOT_TINY {
        QuadRoots::One((-(R::of(c) / b)).stored())
    } else {
        QuadRoots::NoRoot
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RootTarget {
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct RootCurve {
    pub base: Point,
    pub dir: Point,
}

impl RootCurve {
    fn at<R: Real>(self, t: f64) -> Point {
        Point {
            x: (R::of(self.dir.x) * t + self.base.x).stored(),
            y: (R::of(self.dir.y) * t + self.base.y).stored(),
        }
    }
}

impl RootTarget {
    fn value<R: Real>(self, p: Point) -> R {
        R::of(p.x) * self.a + R::of(p.y) * self.b + self.c
    }
}

fn root_sign(f: f64) -> i32 {
    if f > SMALL {
        1
    } else if f >= -SMALL {
        0
    } else {
        -1
    }
}

pub(crate) fn exact_sign(f: f64) -> i32 {
    if f > 0.0 {
        1
    } else if f >= 0.0 {
        0
    } else {
        -1
    }
}

/// The squared distance from the differences as they are in the registers.
pub(crate) fn dist2<R: Real>(a: Point, b: Point) -> R {
    let dx = R::of(a.x) - b.x;
    let dy = R::of(a.y) - b.y;
    dy * dy + dx * dx
}

impl<R: Real> Engine<R> {
    pub(crate) fn seg_start(&self, s: usize) -> Point {
        self.pnt_point(self.segs[s].i1)
    }

    pub(crate) fn seg_end(&self, s: usize) -> Point {
        self.pnt_point(self.segs[s].i2)
    }

    pub(crate) fn seg_has_endpoint(&self, s: usize, p: usize) -> bool {
        self.segs[s].i1 == p || self.segs[s].i2 == p
    }

    /// Parameter of the projection of `p` onto the supporting line of segment `s`.
    pub(crate) fn seg_param(&self, s: usize, p: Point) -> R {
        let seg = self.segs[s];
        let start = self.seg_start(s);
        (R::of(p.x) - start.x) * seg.b - (R::of(p.y) - start.y) * seg.a
    }

    /// Classifies which of two candidate centres project onto segment `s`,
    /// relaxing the tolerance from ZERO up to ZERO_MAX.
    pub(crate) fn qualify_on_seg(&self, s: usize, c1: Point, c2: Point) -> Qualification {
        let lgth = self.segs[s].lgth;
        let p1 = self.seg_param(s, c1).stored();
        let mut eps = ZERO;
        while eps <= ZERO_MAX {
            let hi = (R::of(eps) + lgth).stored();
            let p2 = self.seg_param(s, c2).stored();
            if p1 <= hi && -eps <= p1 {
                return if p2 > hi || p2 < -eps {
                    Qualification::First
                } else {
                    Qualification::Both
                };
            }
            if p2 <= hi && -eps <= p2 {
                return Qualification::Second;
            }
            eps = (R::of(eps) * 10.0).stored();
        }
        Qualification::Neither
    }

    pub(crate) fn qualify_between_nodes(&self, e: usize, c1: Point, c2: Point) -> Qualification {
        let mut eps = ZERO;
        while eps <= ZERO_MAX {
            if self.is_pnt_between_nodes(e, c1.x, c1.y, false, eps) {
                return if self.is_pnt_between_nodes(e, c2.x, c2.y, false, eps) {
                    Qualification::Both
                } else {
                    Qualification::First
                };
            }
            if self.is_pnt_between_nodes(e, c2.x, c2.y, false, eps) {
                return Qualification::Second;
            }
            eps = (R::of(eps) * 10.0).stored();
        }
        Qualification::Neither
    }

    /// Final choice between two candidates after the between-nodes test.
    /// Returns `(choose_second, ambiguous)`.
    pub(crate) fn resolve_candidates(
        &self,
        e: usize,
        c1: Point,
        c2: Point,
        q: Qualification,
        first_was_both: bool,
    ) -> (bool, bool) {
        match q {
            Qualification::First => (false, false),
            Qualification::Second => (true, false),
            Qualification::Both => (self.second_nearer_to_deleted_node(e, c1, c2), false),
            Qualification::Neither => (
                self.second_nearer_to_deleted_node(e, c1, c2),
                !first_was_both,
            ),
        }
    }

    fn second_nearer_to_deleted_node(&self, e: usize, c1: Point, c2: Point) -> bool {
        let mut n = self.edges[e].n1;
        if self.nodes[n].status != NodeStatus::Deleted {
            n = self.edges[e].n2;
        }
        let node = Point {
            x: self.nodes[n].x,
            y: self.nodes[n].y,
        };
        dist2::<R>(node, c2).stored() <= dist2::<R>(node, c1).stored()
    }

    /// Point `k/10` of the way along segment `s`, for k = 0..=10.
    pub(crate) fn seg_samples(&self, s: usize) -> [Point; NUM_SAMPLES] {
        let start = self.seg_start(s);
        let end = self.seg_end(s);
        let dx = (R::of(end.x) - start.x).stored();
        let dy = (R::of(end.y) - start.y).stored();
        std::array::from_fn(|k| {
            let t = (R::of(k as f64) / 10.0).stored();
            Point {
                x: (R::of(t) * dx + start.x).stored(),
                y: (R::of(t) * dy + start.y).stored(),
            }
        })
    }

    /// Brute-force search used by the desperate node routines: circumcircles of
    /// all sample triples, keeping the one that penetrates the probe samples least.
    pub(crate) fn desperate_circle_search(
        &self,
        first: &[Point],
        second: &[Point],
        third: &[Point],
        probes: &[&[Point]],
    ) -> (Point, f64) {
        let mut best = f64::MAX;
        let mut best_center = Point { x: 0.0, y: 0.0 };
        let mut best_r = 0.0;
        let mut last_center = best_center;
        let mut last_r2 = 0.0;
        for &a in first {
            for &b in second {
                for &c in third {
                    let (ok, center, r2) = crate::pnt_pnt_pnt::circumcircle_3pts::<R>(a, b, c);
                    last_center = center;
                    last_r2 = r2;
                    if !ok {
                        continue;
                    }
                    let r = R::of(r2).sqrt().stored();
                    let mut pen = 0.0;
                    for &s in probes.iter().flat_map(|p| p.iter()) {
                        let d = dist2::<R>(center, s).sqrt().stored();
                        if d < r && R::of(r) - d >= pen {
                            pen = (R::of(r) - d).stored();
                        }
                    }
                    if pen < best {
                        best_center = center;
                        best_r = r;
                        best = pen;
                        if (-ZERO..=ZERO).contains(&pen) {
                            return (best_center, best_r);
                        }
                    }
                }
            }
        }
        if best == f64::MAX {
            return (last_center, R::of(last_r2).sqrt().stored());
        }
        (best_center, best_r)
    }

    pub(crate) fn find_root_on_curve(
        &self,
        target: RootTarget,
        curve: RootCurve,
        mut t_s: f64,
        mut t_e: f64,
        t: &mut f64,
    ) -> Option<Point> {
        let probe = |tt: f64| {
            let p = curve.at::<R>(tt);
            (p, target.value::<R>(p).stored())
        };
        let mut t_cur = *t;
        debug_assert!(t_s <= t_e, "t_s <= t_e");
        debug_assert!(
            R::of(t_s) - ZERO <= t_cur && R::of(t_e) + ZERO >= t_cur,
            "(t_s - ZERO <= t) && (t <= t_e + ZERO)"
        );

        let (p, f) = probe(t_cur);
        let mut sign_cur = root_sign(f);
        if sign_cur == 0 {
            *t = t_cur;
            return Some(p);
        }
        let mut best = f;
        let (p_s, f_s) = probe(t_s);
        let mut sign_lo = root_sign(f_s);
        if sign_lo == 0 {
            *t = t_s;
            return Some(p_s);
        }
        let (p_e, f_e) = probe(t_e);
        let sign_e = root_sign(f_e);
        if sign_e == 0 {
            *t = t_e;
            return Some(p_e);
        }
        best = best.abs();
        if f_s.abs() < best {
            best = f_s.abs();
            *t = t_s;
        }
        if f_e.abs() < best {
            best = f_e.abs();
            *t = t_e;
        }

        let mut last_p = p;
        if sign_lo == sign_e {
            if sign_lo == sign_cur {
                let mut step = ((R::of(t_e) - t_s) / ROOT_SAMPLES).stored();
                if step <= ZERO {
                    step = ZERO;
                }
                loop {
                    t_cur = (R::of(step) + t_s).stored();
                    if t_e <= t_cur {
                        let p = curve.at::<R>(*t);
                        if best > ZERO || best < -ZERO {
                            return None;
                        }
                        return Some(p);
                    }
                    let (p, f) = probe(t_cur);
                    last_p = p;
                    sign_cur = root_sign(f);
                    if f.abs() < best {
                        best = f.abs();
                        *t = t_cur;
                    }
                    t_s = t_cur;
                    if sign_cur != sign_lo {
                        break;
                    }
                }
            }
            if sign_cur == 0 {
                *t = t_cur;
                return Some(last_p);
            }
            if *t <= t_cur {
                t_s = *t;
            } else {
                t_s = t_cur;
                t_cur = *t;
                sign_lo = sign_cur;
            }
            let mid = ((R::of(t_s) + t_cur) * 0.5).stored();
            let (_, f) = probe(mid);
            t_e = t_cur;
            t_cur = mid;
            sign_cur = root_sign(f);
        }

        let t_minus = if R::of(t_cur) - ZERO <= t_s {
            t_s
        } else {
            (R::of(t_cur) - ZERO).stored()
        };
        let (p, f) = probe(t_minus);
        let sign_minus = root_sign(f);
        if sign_minus == 0 {
            *t = t_minus;
            return Some(p);
        }
        if f.abs() < best {
            best = f.abs();
            *t = t_minus;
        }
        if sign_minus == sign_cur {
            let t_plus = if R::of(t_cur) + ZERO >= t_e {
                t_e
            } else {
                (R::of(t_cur) + ZERO).stored()
            };
            let (p, f) = probe(t_plus);
            let sign_plus = root_sign(f);
            if sign_plus == 0 {
                *t = t_plus;
                return Some(p);
            }
            if f.abs() < best {
                best = f.abs();
                *t = t_plus;
            }
            if sign_plus != sign_cur {
                t_e = t_plus;
                t_s = t_cur;
                t_cur = ((R::of(t_plus) + t_cur) * 0.5).stored();
                sign_lo = sign_cur;
                sign_cur = root_sign(probe(t_cur).1);
            }
        } else {
            t_s = t_minus;
            t_e = t_cur;
            t_cur = ((R::of(t_minus) + t_cur) * 0.5).stored();
            sign_lo = sign_minus;
            sign_cur = root_sign(probe(t_cur).1);
        }

        let mut iterations = 0;
        loop {
            if sign_cur == sign_lo {
                t_s = t_cur;
            } else {
                t_e = t_cur;
            }
            t_cur = ((R::of(t_s) + t_e) * 0.5).stored();
            let (p, f) = probe(t_cur);
            sign_cur = root_sign(f);
            if sign_cur == 0 {
                *t = t_cur;
                return Some(p);
            }
            if f.abs() < best {
                best = f.abs();
                *t = t_cur;
            }
            iterations += 1;
            if t_cur <= t_s || t_e < t_cur {
                iterations = MAX_ROOT_ITERATIONS;
            }
            if iterations >= MAX_ROOT_ITERATIONS {
                if best > ZERO || best < -ZERO {
                    return None;
                }
                *t = t_cur;
                return Some(p);
            }
        }
    }

    pub(crate) fn seg_seg_bisector(
        &self,
        i: usize,
        j: usize,
        ni: Point,
        nj: Point,
        cj: f64,
    ) -> Option<(Point, Point)> {
        let s = (R::of(ni.x) * nj.x + 1.0 + R::of(ni.y) * nj.y).stored();
        let not_antiparallel = s > PARALLEL_EPS || s < -PARALLEL_EPS;
        let dir = if not_antiparallel {
            Point {
                x: ((R::of(ni.x) + nj.x) * 0.5).stored(),
                y: ((R::of(ni.y) + nj.y) * 0.5).stored(),
            }
        } else {
            Point {
                x: ((R::of(ni.y) - nj.y) * 0.5).stored(),
                y: ((R::of(nj.x) - ni.x) * 0.5).stored(),
            }
        };

        let (si1, si2) = (self.segs[i].i1, self.segs[i].i2);
        if self.seg_has_endpoint(j, si1) {
            return Some((self.pnt_point(si1), dir));
        }
        if self.seg_has_endpoint(j, si2) {
            return Some((self.pnt_point(si2), dir));
        }
        let arc = self.segs[i].apx_arc;
        if arc != NIL && arc == self.segs[j].apx_arc {
            let center = Point {
                x: self.arcs[arc].cx,
                y: self.arcs[arc].cy,
            };
            return Some((center, dir));
        }

        let start = self.seg_start(i);
        let line_j = |p: Point| R::of(p.y) * nj.y + R::of(p.x) * nj.x + cj;
        if !not_antiparallel {
            let t0 = (-((R::of(start.x) * nj.x + R::of(start.y) * nj.y + cj) / (R::of(s) - 1.0)))
                .stored();
            let q = Point {
                x: (R::of(t0) * ni.x + start.x).stored(),
                y: (R::of(t0) * ni.y + start.y).stored(),
            };
            let a = [
                [
                    (R::of(ni.x) + ni.y).stored(),
                    (-(R::of(ni.y) + nj.x)).stored(),
                ],
                [
                    (R::of(ni.y) - ni.x).stored(),
                    (-(R::of(nj.y) - ni.x)).stored(),
                ],
            ];
            let (x, _) = solve_2x2::<R>(
                a,
                [
                    (R::of(q.x) - start.x).stored(),
                    (R::of(q.y) - start.y).stored(),
                ],
            )?;
            let pnt = Point {
                x: (R::of(a[0][0]) * x[0] + start.x).stored(),
                y: (R::of(x[0]) * a[1][0] + start.y).stored(),
            };
            return Some((pnt, dir));
        }

        let cross = (R::of(ni.x) * nj.y - R::of(ni.y) * nj.x).stored();
        if cross > PARALLEL_EPS || cross < -PARALLEL_EPS {
            let t = (line_j(start) / cross).stored();
            let pnt = Point {
                x: (R::of(t) * ni.y + start.x).stored(),
                y: (R::of(t) * -ni.x + start.y).stored(),
            };
            return Some((pnt, dir));
        }

        let seg = self.segs[i];
        let u = Point {
            x: seg.b,
            y: -seg.a,
        };
        let mut t_s = 0.0;
        let mut t_e = seg.lgth;
        let mut p0 = start;
        let mut p1 = self.seg_end(i);
        let mut sign0 = exact_sign(line_j(p0).stored());
        let mut sign1 = exact_sign(line_j(p1).stored());
        if sign0 == 0 {
            return Some((p0, dir));
        }
        if sign1 == 0 {
            return Some((p1, dir));
        }
        let mut expansions = 0;
        while sign0 == sign1 && expansions < MAX_EXPANSIONS {
            let d = (R::of(t_e) - t_s).stored();
            t_s = (R::of(t_s) - d).stored();
            t_e = (R::of(t_e) + d).stored();
            p0 = Point {
                x: (R::of(t_s) * u.x + start.x).stored(),
                y: (R::of(t_s) * u.y + start.y).stored(),
            };
            p1 = Point {
                x: (R::of(t_e) * u.x + start.x).stored(),
                y: (R::of(t_e) * u.y + start.y).stored(),
            };
            sign0 = exact_sign(line_j(p0).stored());
            sign1 = exact_sign(line_j(p1).stored());
            if sign0 == 0 {
                return Some((p0, dir));
            }
            if sign1 == 0 {
                return Some((p1, dir));
            }
            expansions += 1;
        }
        if sign0 != sign1 {
            let mut t = ((R::of(t_s) + t_e) * 0.5).stored();
            let target = RootTarget {
                a: nj.x,
                b: nj.y,
                c: cj,
            };
            let curve = RootCurve {
                base: start,
                dir: u,
            };
            if let Some(pnt) = self.find_root_on_curve(target, curve, t_s, t_e, &mut t) {
                return Some((pnt, dir));
            }
        }

        let a0 = self.seg_start(i);
        let a1 = self.seg_end(i);
        let b0 = self.seg_start(j);
        let b1 = self.seg_end(j);
        let d00 = dist2::<R>(a0, b0).stored();
        let d01 = dist2::<R>(a0, b1).stored();
        let (mut best_b, best) = if d01 <= d00 { (b1, d01) } else { (b0, d00) };
        let d10 = dist2::<R>(a1, b0).stored();
        let d11 = dist2::<R>(a1, b1).stored();
        let best_a = if d11 <= d10 {
            if best <= d11 {
                a0
            } else {
                best_b = b1;
                a1
            }
        } else if best <= d10 {
            a0
        } else {
            best_b = b0;
            a1
        };
        let pnt = Point {
            x: ((R::of(best_a.x) + best_b.x) * 0.5).stored(),
            y: ((R::of(best_b.y) + best_a.y) * 0.5).stored(),
        };
        Some((pnt, dir))
    }

    pub(crate) fn pnt_ray_cntr(&mut self, i1: usize, v: Point, e: usize) -> Option<(Point, f64)> {
        let edge = self.edges[e];
        let lft_is_i1 = edge.lft == i1 && edge.ltype == SiteType::Pnt;
        debug_assert!(
            lft_is_i1 || (edge.rgt == i1 && edge.rtype == SiteType::Pnt),
            "IsLftRgtSite(e, i1, PNT)"
        );
        debug_assert!(
            ((v.x * v.x + v.y * v.y).sqrt() - 1.0).abs() <= ZERO,
            "eq(VecLen(v), 1.0)"
        );
        let (other, other_type) = if lft_is_i1 {
            (edge.rgt, edge.rtype)
        } else {
            (edge.lft, edge.ltype)
        };
        let p = self.pnt_point(i1);
        let on_ray = |t: f64| Point {
            x: (R::of(v.x) * t + p.x).stored(),
            y: (R::of(v.y) * t + p.y).stored(),
        };
        match other_type {
            SiteType::Seg => {
                if self.seg_has_endpoint(other, i1) {
                    return Some((p, 0.0));
                }
                let seg = self.segs[other];
                let d = (R::of(seg.a) * p.x + R::of(p.y) * seg.b + seg.c).stored();
                if d > ZERO || d < -ZERO {
                    let den = (R::of(seg.a) * v.x + R::of(seg.b) * v.y - f64::from(exact_sign(d)))
                        .stored();
                    if den > ZERO || den < -ZERO {
                        let t = (-(R::of(d) / den)).stored();
                        if t >= -ZERO {
                            return Some((on_ray(t), t.abs()));
                        }
                    }
                } else {
                    self.restart = true;
                }
                None
            }
            SiteType::Pnt => {
                let q = self.pnt_point(other);
                let w = Point {
                    x: (R::of(q.x) - p.x).stored(),
                    y: (R::of(q.y) - p.y).stored(),
                };
                let den = (R::of(w.x) * v.x + R::of(w.y) * v.y).stored();
                if den > ZERO || den < -ZERO {
                    let t = ((R::of(w.x) * w.x + R::of(w.y) * w.y) / (R::of(den) + den)).stored();
                    if t >= -ZERO {
                        return Some((on_ray(t), t.abs()));
                    }
                }
                None
            }
            SiteType::Arc | SiteType::NoSite => None,
        }
    }
}
