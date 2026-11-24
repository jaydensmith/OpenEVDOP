use crate::consts::*;
use crate::numerics::{QuadRoots, Qualification, exact_sign, solve_quadratic};
use crate::state::*;

impl Engine {
    /// Orients the line of segment `m` so that `pk` (or, if `pk` lies on the
    /// line, the reference node or segment `other`) is on its positive side.
    fn orient_seg_line(
        &mut self,
        m: usize,
        pk: Point,
        refp: Point,
        other: usize,
    ) -> (f64, f64, f64) {
        let seg = self.segs[m];
        let (a, b, c) = (seg.a, seg.b, seg.c);
        let value = |p: Point| p.x * a + p.y * b + c;
        let d = value(pk);
        let flip = if d > ZERO || d < -ZERO {
            d < 0.0
        } else {
            self.restart = true;
            let d_ref = value(refp);
            if d_ref > ZERO || d_ref < -ZERO {
                d_ref < 0.0
            } else {
                let e1 = value(self.seg_start(other));
                let e2 = value(self.seg_end(other));
                if e1 < 0.0 && e2 < 0.0 {
                    true
                } else if e1 > 0.0 && e2 > 0.0 {
                    false
                } else {
                    d_ref < 0.0
                }
            }
        };
        if flip { (-a, -b, -c) } else { (a, b, c) }
    }

    pub(crate) fn seg_seg_pnt_cntr(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        e: usize,
        out: &mut NodeCandidate,
    ) -> bool {
        debug_assert!(self.segs[i].lgth > ZERO, "gt(GetSegLgth(i), ZERO)");
        debug_assert!(self.segs[j].lgth > ZERO, "gt(GetSegLgth(j), ZERO)");
        out.coincident = NIL;
        if i == j {
            return false;
        }
        if self.seg_has_endpoint(j, self.segs[i].i1) && self.seg_has_endpoint(j, self.segs[i].i2) {
            return false;
        }
        let pk = self.pnt_point(k);
        let (i, j) = if self.seg_has_endpoint(j, k) {
            (j, i)
        } else {
            (i, j)
        };
        if self.seg_has_endpoint(j, k) {
            out.center = pk;
            out.r = 0.0;
            out.coincident = k;
            return true;
        }
        if self.seg_has_endpoint(i, k) {
            self.seg_seg_pnt_at_endpoint(i, j, pk, e, out)
        } else {
            self.seg_seg_pnt_general(i, j, pk, e, out)
        }
    }

    /// `pk` is an endpoint of segment `i`: the centre lies on the normal of `i` through `pk`.
    fn seg_seg_pnt_at_endpoint(
        &self,
        i: usize,
        j: usize,
        pk: Point,
        e: usize,
        out: &mut NodeCandidate,
    ) -> bool {
        let si = self.segs[i];
        let sj = self.segs[j];
        let d = pk.x * sj.a + pk.y * sj.b + sj.c;
        if (-ZERO..=ZERO).contains(&d) {
            return false;
        }
        let q = sj.b * si.b + sj.a * si.a;
        let mut roots = [0.0; 2];
        let mut count = 0;
        if q > ZERO || q < -ZERO {
            let s = f64::from(exact_sign(d));
            let den1 = s + q;
            if den1 > ZERO || den1 < -ZERO {
                roots[0] = -(d / den1);
                count = 1;
            }
            let den2 = s - q;
            if den2 > ZERO || den2 < -ZERO {
                roots[count] = d / den2;
                count += 1;
            }
        } else {
            roots = [-d, d];
            count = 2;
        }
        let at = |t: f64| Point {
            x: t * si.a + pk.x,
            y: si.b * t + pk.y,
        };
        let t = match count {
            2 => {
                let c1 = at(roots[0]);
                let c2 = at(roots[1]);
                let qual = self.qualify_on_seg(j, c1, c2);
                let (second, ambiguous) = match qual {
                    Qualification::First => (false, false),
                    Qualification::Second => (true, false),
                    Qualification::Both | Qualification::Neither => {
                        let q2 = self.qualify_between_nodes(e, c1, c2);
                        self.resolve_candidates(e, c1, c2, q2, qual == Qualification::Both)
                    }
                };
                if ambiguous {
                    out.ambiguous = true;
                }
                if second { roots[1] } else { roots[0] }
            }
            1 => roots[0],
            _ => return false,
        };
        out.center = at(t);
        out.r = t.abs();
        true
    }

    /// `pk` is not an endpoint of either segment: the centre lies on the
    /// bisector of the two oriented lines and is equidistant to `pk`.
    fn seg_seg_pnt_general(
        &mut self,
        i: usize,
        j: usize,
        pk: Point,
        e: usize,
        out: &mut NodeCandidate,
    ) -> bool {
        let (i, j) = if j < i { (j, i) } else { (i, j) };
        let mut n = self.edges[e].n1;
        if self.nodes[n].status == NodeStatus::Deleted {
            n = self.edges[e].n2;
        }
        let refp = Point {
            x: self.nodes[n].x,
            y: self.nodes[n].y,
        };
        let (ai, bi, ci) = self.orient_seg_line(i, pk, refp, j);
        let (aj, bj, cj) = self.orient_seg_line(j, pk, refp, i);
        if self.restart {
            return false;
        }
        let Some((p, d)) =
            self.seg_seg_bisector(i, j, Point { x: ai, y: bi }, Point { x: aj, y: bj }, cj)
        else {
            return false;
        };
        let wx = p.x - pk.x;
        let wy = p.y - pk.y;
        let dp = (p.y * bi + p.x * ai) + ci;
        let q = ai * d.x + bi * d.y;
        let coef_a = (d.x * d.x + d.y * d.y) - q * q;
        let coef_b = ((wy * d.y + wx * d.x) - q * dp) * 2.0;
        let coef_c = (wy * wy + wx * wx) - dp * dp;
        let at = |t: f64| Point {
            x: d.x * t + p.x,
            y: t * d.y + p.y,
        };
        let center = match solve_quadratic(coef_a, coef_b, coef_c) {
            QuadRoots::NoRoot => return false,
            QuadRoots::One(t) => at(t),
            QuadRoots::Two(t1, t2) => {
                let c1 = at(t1);
                let c2 = at(t2);
                let qi = self.qualify_on_seg(i, c1, c2);
                let qj = self.qualify_on_seg(j, c1, c2);
                let first_was_both = match (qi, qj) {
                    (Qualification::Neither, _) | (_, Qualification::Neither) => Some(false),
                    (Qualification::Both, _)
                    | (_, Qualification::Both)
                    | (Qualification::First, Qualification::Second)
                    | (Qualification::Second, Qualification::First) => Some(true),
                    _ => None,
                };
                let second = match first_was_both {
                    None => qi == Qualification::Second,
                    Some(fwb) => {
                        let q2 = self.qualify_between_nodes(e, c1, c2);
                        let (second, ambiguous) = self.resolve_candidates(e, c1, c2, q2, fwb);
                        if ambiguous {
                            out.ambiguous = true;
                        }
                        second
                    }
                };
                if second { c2 } else { c1 }
            }
        };
        out.center = center;
        out.r = (center.x * ai + center.y * bi + ci).abs();
        true
    }

    pub(crate) fn seg_seg_pnt_desperate(&self, i: usize, j: usize, pk: Point) -> (Point, f64) {
        let si = self.seg_samples(i);
        let sj = self.seg_samples(j);
        self.desperate_circle_search(&[pk], &si, &sj, &[&si, &sj])
    }
}
