use crate::consts::*;
use crate::numerics::{QuadRoots, Qualification, solve_quadratic};
use crate::state::*;

impl Engine {
    pub(crate) fn seg_pnt_cntr(&self, i: usize, p: Point) -> Option<(Point, f64)> {
        let seg = self.segs[i];
        debug_assert!(seg.lgth > ZERO, "gt(GetSegLgth(i), ZERO)");
        let t = self.seg_param(i, p);
        if t > seg.lgth + ZERO || t < -ZERO {
            return None;
        }
        let h = (p.x * seg.a + p.y * seg.b + seg.c) * 0.5;
        let center = Point {
            x: p.x - seg.a * h,
            y: p.y - seg.b * h,
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
        let d_k = pk.x * a + pk.y * b + seg.c;
        if d_k <= ZERO && -ZERO < d_k {
            let t = self.seg_param(i, pk);
            if t <= seg.lgth + ZERO && -ZERO <= t {
                if self.seg_has_endpoint(i, k) {
                    out.center = Point {
                        x: (pj.x + pk.x) * 0.5,
                        y: (pj.y + pk.y) * 0.5,
                    };
                    out.r = 0.0;
                    out.ambiguous = true;
                    return true;
                }
                return false;
            }
        }
        let wx = pj.x - pk.x;
        let wy = pj.y - pk.y;
        let w2 = wx * wx + wy * wy;
        if w2 > ZERO_SQ || w2 < -ZERO_SQ {
            if d_k > 0.0 {
                a = -a;
                b = -b;
            }
            let den = wx * a + wy * b;
            if den > ZERO || den < -ZERO {
                let t = ((wx * wx + wy * wy) * -0.5) / den;
                out.center = Point {
                    x: t * a + pj.x,
                    y: t * b + pj.y,
                };
                out.r = t.abs();
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
        let dx = pk.x - pj.x;
        let dy = pk.y - pj.y;
        let w2 = dx * dx + dy * dy;
        let v = Point {
            x: -(dy * 0.5),
            y: dx * 0.5,
        };
        let h2 = w2 * 0.25;

        let d_j = seg.a * pj.x + seg.b * pj.y + seg.c;
        if (-ZERO..=ZERO).contains(&d_j) {
            let t = self.seg_param(i, pj);
            if t < seg.lgth && 0.0 < t {
                return None;
            }
        }
        let d_k = seg.a * pk.x + seg.b * pk.y + seg.c;
        if (-ZERO..=ZERO).contains(&d_k) {
            let t = self.seg_param(i, pk);
            if t < seg.lgth && 0.0 < t {
                return None;
            }
        }

        let d_mid = (d_k + d_j) * 0.5;
        let q = seg.b * v.y + seg.a * v.x;
        let coef_b = d_mid * 2.0 * q;
        let coef_a = q * q - h2;
        let coef_c = d_mid * d_mid - h2;
        let m = Point {
            x: (pj.x + pk.x) * 0.5,
            y: (pj.y + pk.y) * 0.5,
        };
        let at = |t: f64| Point {
            x: v.x * t + m.x,
            y: t * v.y + m.y,
        };

        let (center, ambiguous) = match solve_quadratic(coef_a, coef_b, coef_c) {
            QuadRoots::NoRoot => return None,
            QuadRoots::One(t) => (at(t), false),
            QuadRoots::Two(t1, t2) => {
                let c1 = at(t1);
                let c2 = at(t2);
                let qual = self.qualify_on_seg(i, c1, c2);
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
        let r = (center.y * seg.b + center.x * seg.a + seg.c).abs();
        Some((center, r, ambiguous))
    }

    pub(crate) fn seg_pnt_pnt_desperate(&self, i: usize, pj: Point, pk: Point) -> (Point, f64) {
        let samples = self.seg_samples(i);
        self.desperate_circle_search(&[pj], &[pk], &samples, &[&samples])
    }
}
