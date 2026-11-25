use crate::consts::*;
use crate::numerics::{QuadRoots, Qualification, solve_quadratic};
use crate::real::Real;
use crate::state::*;

impl<R: Real> Engine<R> {
    pub(crate) fn seg_pnt_cntr(&self, i: usize, p: Point) -> Option<(Point, f64)> {
        let seg = self.segs[i];
        debug_assert!(seg.lgth > ZERO, "gt(GetSegLgth(i), ZERO)");
        let t = self.seg_param(i, p).stored();
        if R::of(seg.lgth) + ZERO < t || t < -ZERO {
            return None;
        }
        let h = ((R::of(p.x) * seg.a + R::of(p.y) * seg.b + seg.c) * 0.5).stored();
        let center = Point {
            x: (R::of(p.x) - R::of(seg.a) * h).stored(),
            y: (R::of(p.y) - R::of(seg.b) * h).stored(),
        };
        Some((center, h.abs()))
    }

    pub(crate) fn seg_pnt_pnt_cntr(
        &self,
        i: usize,
        j: usize,
        k: usize,
        e: usize,
        out: &mut NodeCandidate,
    ) -> bool {
        debug_assert!(self.segs[i].lgth > ZERO, "gt(GetSegLgth(i), ZERO)");
        if j == k {
            return match self.seg_pnt_cntr(i, self.pnt_point(j)) {
                Some((center, r)) => {
                    out.center = center;
                    out.r = r;
                    true
                }
                None => false,
            };
        }
        let (j, k) = if self.seg_has_endpoint(i, k) {
            (k, j)
        } else {
            (j, k)
        };
        if !self.seg_has_endpoint(i, j) {
            let (lo, hi) = if k < j { (k, j) } else { (j, k) };
            return match self.seg_pnt_pnt_general(i, self.pnt_point(lo), self.pnt_point(hi), e) {
                Some((center, r, ambiguous)) => {
                    out.center = center;
                    out.r = r;
                    if ambiguous {
                        out.ambiguous = true;
                    }
                    true
                }
                None => false,
            };
        }

        let pj = self.pnt_point(j);
        let pk = self.pnt_point(k);
        let seg = self.segs[i];
        let (mut a, mut b) = (seg.a, seg.b);
        let d_k = (R::of(pk.x) * a + R::of(pk.y) * b + seg.c).stored();
        if d_k <= ZERO && -ZERO < d_k {
            let t = self.seg_param(i, pk).stored();
            if R::of(seg.lgth) + ZERO >= t && -ZERO <= t {
                if self.seg_has_endpoint(i, k) {
                    out.center = Point {
                        x: ((R::of(pj.x) + pk.x) * 0.5).stored(),
                        y: ((R::of(pj.y) + pk.y) * 0.5).stored(),
                    };
                    out.r = 0.0;
                    out.ambiguous = true;
                    return true;
                }
                return false;
            }
        }
        let wx = R::of(pj.x) - pk.x;
        let wy = R::of(pj.y) - pk.y;
        let w2 = wx * wx + wy * wy;
        if w2 > ZERO_SQ || w2 < -ZERO_SQ {
            if d_k > 0.0 {
                a = -a;
                b = -b;
            }
            let den = wx * a + wy * b;
            if den > ZERO || den < -ZERO {
                let t = (wx * wx + wy * wy) * -0.5 / den;
                out.center = Point {
                    x: (t * a + pj.x).stored(),
                    y: (t * b + pj.y).stored(),
                };
                out.r = t.abs().stored();
                return true;
            }
        }
        false
    }

    /// Node equidistant to segment `i` and two points that are not its endpoints:
    /// the centre lies on the perpendicular bisector of `pj`,`pk`.
    pub(crate) fn seg_pnt_pnt_general(
        &self,
        i: usize,
        pj: Point,
        pk: Point,
        e: usize,
    ) -> Option<(Point, f64, bool)> {
        let seg = self.segs[i];
        debug_assert!(seg.lgth > ZERO, "gt(GetSegLgth(i), ZERO)");
        let dx = R::of(pk.x) - pj.x;
        let dy = R::of(pk.y) - pj.y;
        let w2 = dx * dx + dy * dy;
        let v = Point {
            x: (-(dy * 0.5)).stored(),
            y: (dx * 0.5).stored(),
        };
        let h2 = (w2 * 0.25).stored();

        let d_j = (R::of(seg.a) * pj.x + R::of(seg.b) * pj.y + seg.c).stored();
        if (-ZERO..=ZERO).contains(&d_j) {
            let t = self.seg_param(i, pj).stored();
            if t < seg.lgth && 0.0 < t {
                return None;
            }
        }
        let d_k = (R::of(seg.a) * pk.x + R::of(seg.b) * pk.y + seg.c).stored();
        if (-ZERO..=ZERO).contains(&d_k) {
            let t = self.seg_param(i, pk).stored();
            if t < seg.lgth && 0.0 < t {
                return None;
            }
        }

        let d_mid = (R::of(d_k) + d_j) * 0.5;
        let q = R::of(seg.b) * v.y + R::of(seg.a) * v.x;
        let coef_b = d_mid * 2.0 * q;
        let coef_a = q * q - h2;
        let coef_c = d_mid * d_mid - h2;
        let m = Point {
            x: ((R::of(pj.x) + pk.x) * 0.5).stored(),
            y: ((R::of(pj.y) + pk.y) * 0.5).stored(),
        };
        let at = |t: R| Point {
            x: (R::of(v.x) * t + m.x).stored(),
            y: (t * v.y + m.y).stored(),
        };

        let (center, ambiguous) = match solve_quadratic(coef_a, coef_b, coef_c) {
            QuadRoots::NoRoot => return None,
            QuadRoots::One(t) => (at(t), false),
            QuadRoots::Two(t1, t2) => {
                let c1 = at(t1);
                let c2 = at(t2);
                let qual = self.qualify_on_seg(i, self.seg_param(i, c1), self.seg_param(i, c2));
                let (second, ambiguous) = match qual {
                    Qualification::First => (false, false),
                    Qualification::Second => (true, false),
                    Qualification::Both | Qualification::Neither => {
                        let q2 = self.qualify_between_nodes(e, c1, c2);
                        self.resolve_candidates(e, c1, c2, q2, qual == Qualification::Both)
                    }
                };
                (if second { c2 } else { c1 }, ambiguous)
            }
        };
        let r = (R::of(center.y) * seg.b + R::of(center.x) * seg.a + seg.c)
            .abs()
            .stored();
        Some((center, r, ambiguous))
    }

    pub(crate) fn seg_pnt_pnt_desperate(&self, i: usize, pj: Point, pk: Point) -> (Point, f64) {
        let samples = self.seg_samples(i);
        self.desperate_circle_search(&[pj], &[pk], &samples, &[&samples])
    }
}
