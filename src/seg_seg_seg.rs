use crate::consts::*;
use crate::numerics::solve_2x2;
use crate::real::Real;
use crate::state::*;

#[derive(Debug, Clone, Copy)]
struct Line {
    a: f64,
    b: f64,
    c: f64,
}

impl Line {
    fn value<R: Real>(self, p: Point) -> R {
        R::of(p.x) * self.a + R::of(p.y) * self.b + self.c
    }

    fn flipped(self) -> Line {
        Line {
            a: -self.a,
            b: -self.b,
            c: -self.c,
        }
    }

    fn normal(self) -> Point {
        Point {
            x: self.a,
            y: self.b,
        }
    }

    fn cos<R: Real>(self, other: Line) -> R {
        R::of(other.a) * self.a + R::of(other.b) * self.b
    }
}

impl<R: Real> Engine<R> {
    fn seg_line(&self, s: usize) -> Line {
        let seg = self.segs[s];
        Line {
            a: seg.a,
            b: seg.b,
            c: seg.c,
        }
    }

    fn is_seg_site_of_edge(&self, e: usize, s: usize) -> bool {
        let edge = self.edges[e];
        (edge.lft == s && edge.ltype == SiteType::Seg)
            || (edge.rgt == s && edge.rtype == SiteType::Seg)
    }

    fn segs_share_endpoint(&self, s: usize, t: usize) -> bool {
        self.seg_has_endpoint(t, self.segs[s].i1) || self.seg_has_endpoint(t, self.segs[s].i2)
    }

    /// Orients the line of segment `m` so that edge `e` lies on its positive side,
    /// falling back to the endpoints of the other two segments when both edge
    /// nodes lie on the line.
    fn orient_sss_line(
        &mut self,
        e: usize,
        m: usize,
        node_a: Point,
        node_b: Point,
        others: [usize; 2],
    ) -> Line {
        let line = self.seg_line(m);
        let nonzero = |d: R| d > ZERO || d < -ZERO;
        let (first, second) = if self.is_seg_site_of_edge(e, m) {
            (node_b, node_a)
        } else {
            (node_a, node_b)
        };
        let d1 = line.value::<R>(first);
        let flip = if nonzero(d1) {
            d1 < 0.0
        } else {
            let d2 = line.value::<R>(second);
            if nonzero(d2) {
                d2 < 0.0
            } else {
                self.restart = true;
                self.orient_by_endpoints(line, others, d2)
            }
        };
        if flip { line.flipped() } else { line }
    }

    fn orient_by_endpoints(&mut self, line: Line, others: [usize; 2], fallback: R) -> bool {
        for s in others {
            let e1 = line.value::<R>(self.seg_start(s));
            let e2 = line.value::<R>(self.seg_end(s));
            if e1 < 0.0 && e2 < 0.0 {
                return true;
            }
            if e1 > 0.0 && e2 > 0.0 {
                return false;
            }
        }
        self.restart = true;
        fallback < 0.0
    }

    pub(crate) fn seg_seg_seg_cntr(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        e: usize,
        out: &mut NodeCandidate,
    ) -> bool {
        debug_assert!(self.segs[i].lgth > ZERO, "gt(GetSegLgth(i), ZERO)");
        debug_assert!(self.segs[j].lgth > ZERO, "gt(GetSegLgth(j), ZERO)");
        debug_assert!(self.segs[k].lgth > ZERO, "gt(GetSegLgth(k), ZERO)");
        out.coincident = NIL;
        if i == j || j == k {
            out.center = self.seg_start(i);
            out.r = 0.0;
            return false;
        }
        let arc = self.segs[i].apx_arc;
        if arc != NIL && arc == self.segs[j].apx_arc && arc == self.segs[k].apx_arc {
            out.center = Point {
                x: self.arcs[arc].cx,
                y: self.arcs[arc].cy,
            };
            out.r = self.arcs[arc].r;
            return true;
        }

        let (k1, k2) = (self.segs[k].i1, self.segs[k].i2);
        if self.seg_has_endpoint(i, k1) {
            if self.seg_has_endpoint(j, k1) {
                self.set_coincident(out, k1);
                return !(self.seg_has_endpoint(i, k2) || self.seg_has_endpoint(j, k2));
            }
            if self.seg_has_endpoint(i, k2) {
                return false;
            }
        } else if self.seg_has_endpoint(i, k2) {
            if self.seg_has_endpoint(j, k2) {
                self.set_coincident(out, k2);
                return !self.seg_has_endpoint(j, k1);
            }
        } else if self.seg_has_endpoint(j, k2) && self.seg_has_endpoint(j, k1) {
            return false;
        }
        self.seg_seg_seg_general(i, j, k, e, out)
    }

    fn set_coincident(&self, out: &mut NodeCandidate, p: usize) {
        out.center = self.pnt_point(p);
        out.r = 0.0;
        out.coincident = p;
    }

    fn seg_seg_seg_general(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        e: usize,
        out: &mut NodeCandidate,
    ) -> bool {
        let new_site = k;
        let mut idx = [i, j, k];
        idx.sort_unstable();
        let [i, j, k] = idx;

        let edge = self.edges[e];
        let (na, nb) = if self.nodes[edge.n1].status == NodeStatus::Deleted {
            (edge.n2, edge.n1)
        } else {
            (edge.n1, edge.n2)
        };
        let node_a = Point {
            x: self.nodes[na].x,
            y: self.nodes[na].y,
        };
        let node_b = Point {
            x: self.nodes[nb].x,
            y: self.nodes[nb].y,
        };
        let (r_lo, r_hi) = if self.nodes[nb].r <= self.nodes[na].r {
            (self.nodes[nb].r, self.nodes[na].r)
        } else {
            (self.nodes[na].r, self.nodes[nb].r)
        };

        let mut li = self.orient_sss_line(e, i, node_a, node_b, [j, k]);
        let mut lj = self.orient_sss_line(e, j, node_a, node_b, [i, k]);
        let mut lk = self.orient_sss_line(e, k, node_a, node_b, [i, j]);
        if self.restart {
            return false;
        }

        let share_ij = self.segs_share_endpoint(i, j);
        let share_ik = self.segs_share_endpoint(i, k);
        let share_jk = self.segs_share_endpoint(j, k);
        let acceptable = |s: &Self, center: Point, r: R| {
            !(r + ZERO_MAX < r_lo || R::of(r_hi) + ZERO_MAX < r)
                && s.is_pnt_between_nodes(e, center.x, center.y, false, ZERO_MAX)
        };
        let mut previous: Option<(Point, f64)> = None;

        loop {
            let mut near_singular: Option<(Point, f64)> = None;
            let a = [
                [(R::of(li.a) - lj.a).stored(), (R::of(li.b) - lj.b).stored()],
                [(R::of(li.a) - lk.a).stored(), (R::of(li.b) - lk.b).stored()],
            ];
            let rhs = [(R::of(lj.c) - li.c).stored(), (R::of(lk.c) - li.c).stored()];
            if let Some((x, pivot)) = solve_2x2::<R>(a, rhs) {
                let center = Point { x: x[0], y: x[1] };
                let r = li.value::<R>(center).abs();
                out.center = center;
                out.r = r.stored();
                if acceptable(self, center, r) {
                    if pivot > ZERO_MAX || pivot < -ZERO_MAX {
                        return true;
                    }
                    near_singular = Some((center, r.stored()));
                }
            }

            let c_ik = lk.cos::<R>(li).stored();
            let c_ij = lj.cos::<R>(li).stored();
            let c_jk = lk.cos::<R>(lj).stored();
            let (v_ik, v_ij, v_jk, max) = match (share_ij, share_ik, share_jk) {
                (false, false, true) => (c_ik, c_ij, 1.0, 1.0),
                (false, true, _) => (1.0, c_ij, if share_jk { -1.0 } else { c_jk }, 1.0),
                (true, false, _) | (true, true, false) => (
                    if share_ik { -1.0 } else { c_ik },
                    1.0,
                    if share_jk { -1.0 } else { c_jk },
                    1.0,
                ),
                (false, false, false) | (true, true, true) => {
                    let max = if c_ij <= c_ik {
                        if c_jk <= c_ik { c_ik } else { c_jk }
                    } else if c_jk <= c_ij {
                        c_ij
                    } else {
                        c_jk
                    };
                    (c_ik, c_ij, c_jk, max)
                }
            };
            let pair_for = |v: f64| {
                if v == v_ik {
                    (i, li, k, lk)
                } else if v == v_ij {
                    (i, li, j, lj)
                } else {
                    (j, lj, k, lk)
                }
            };
            let min = if v_ij <= v_ik {
                if v_jk <= v_ij { v_jk } else { v_ij }
            } else if v_jk <= v_ik {
                v_jk
            } else {
                v_ik
            };

            let bisector = |s: &Self, (s1, l1, s2, l2): (usize, Line, usize, Line)| {
                s.seg_seg_bisector(s1, s2, l1.normal(), l2.normal(), l2.c)
            };
            let bisectors = bisector(self, pair_for(max))
                .and_then(|b1| bisector(self, pair_for(min)).map(|b2| (b1, b2)));
            let Some(((p1, d1), (p2, d2))) = bisectors else {
                return match near_singular {
                    Some((center, r)) => {
                        out.center = center;
                        out.r = r;
                        true
                    }
                    None => false,
                };
            };
            let a = [[d1.x, -d2.x], [d1.y, -d2.y]];
            let rhs = [(R::of(p2.x) - p1.x).stored(), (R::of(p2.y) - p1.y).stored()];
            if let Some((x, _)) = solve_2x2::<R>(a, rhs) {
                let center = Point {
                    x: (R::of(d1.x) * x[0] + p1.x).stored(),
                    y: (R::of(d1.y) * x[0] + p1.y).stored(),
                };
                let r = li.value::<R>(center).abs();
                out.center = center;
                out.r = r.stored();
                if acceptable(self, center, r) {
                    return true;
                }
            }
            if let Some((center, r)) = near_singular {
                out.center = center;
                out.r = r;
                return true;
            }

            match previous {
                None => {
                    debug_assert!(
                        !self.is_seg_site_of_edge(e, new_site),
                        "!IsLftRgtSite(e, K, SEG)"
                    );
                    previous = Some((out.center, out.r));
                    if new_site == i {
                        li = li.flipped();
                    } else if new_site == j {
                        lj = lj.flipped();
                    } else {
                        lk = lk.flipped();
                    }
                }
                Some((prev_center, prev_r)) => {
                    let deviation = |r: f64| {
                        if r_lo <= r {
                            R::of(r) - r_hi
                        } else {
                            R::of(r_lo) - r
                        }
                    };
                    let keep_current = deviation(out.r) <= deviation(prev_r);
                    if !keep_current {
                        out.center = prev_center;
                        out.r = prev_r;
                    }
                    return true;
                }
            }
        }
    }

    pub(crate) fn seg_seg_seg_desperate(&self, i: usize, j: usize, k: usize) -> (Point, f64) {
        let si = self.seg_samples(i);
        let sj = self.seg_samples(j);
        let sk = self.seg_samples(k);
        self.desperate_circle_search(&si, &sj, &sk, &[&si, &sj, &sk])
    }
}
