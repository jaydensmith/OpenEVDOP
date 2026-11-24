use crate::consts::*;
use crate::data::dist;
use crate::numerics::{RootCurve, RootTarget};
use crate::state::*;

const QSORT_CUTOFF: isize = 8;

/// Unstable quicksort with insertion sort for short ranges; its element order decides which duplicate
/// record survives in clean_data.
pub(crate) fn quicksort<T: Copy>(v: &mut [T], cmp: impl Fn(&T, &T) -> i32) {
    if v.len() < 2 {
        return;
    }
    let c = |v: &[T], a: isize, b: isize| cmp(&v[a as usize], &v[b as usize]);
    let swap = |v: &mut [T], a: isize, b: isize| v.swap(a as usize, b as usize);
    let mut stack: Vec<(isize, isize)> = Vec::new();
    let mut lo: isize = 0;
    let mut hi: isize = v.len() as isize - 1;
    loop {
        let size = hi - lo + 1;
        if size <= QSORT_CUTOFF {
            shortsort(v, lo, hi, &c);
        } else {
            let mut mid = lo + size / 2;
            if c(v, lo, mid) > 0 {
                swap(v, lo, mid);
            }
            if c(v, lo, hi) > 0 {
                swap(v, lo, hi);
            }
            if c(v, mid, hi) > 0 {
                swap(v, mid, hi);
            }
            let mut loguy = lo;
            let mut higuy = hi;
            loop {
                if mid > loguy {
                    loop {
                        loguy += 1;
                        if !(loguy < mid && c(v, loguy, mid) <= 0) {
                            break;
                        }
                    }
                }
                if mid <= loguy {
                    loop {
                        loguy += 1;
                        if !(loguy <= hi && c(v, loguy, mid) <= 0) {
                            break;
                        }
                    }
                }
                loop {
                    higuy -= 1;
                    if !(higuy > mid && c(v, higuy, mid) > 0) {
                        break;
                    }
                }
                if higuy < loguy {
                    break;
                }
                swap(v, loguy, higuy);
                if mid == higuy {
                    mid = loguy;
                }
            }
            higuy += 1;
            if mid < higuy {
                loop {
                    higuy -= 1;
                    if !(higuy > mid && c(v, higuy, mid) == 0) {
                        break;
                    }
                }
            }
            if mid >= higuy {
                loop {
                    higuy -= 1;
                    if !(higuy > lo && c(v, higuy, mid) == 0) {
                        break;
                    }
                }
            }
            if higuy - lo >= hi - loguy {
                if lo < higuy {
                    stack.push((lo, higuy));
                }
                if loguy < hi {
                    lo = loguy;
                    continue;
                }
            } else {
                if loguy < hi {
                    stack.push((loguy, hi));
                }
                if lo < higuy {
                    hi = higuy;
                    continue;
                }
            }
        }
        match stack.pop() {
            Some((l, h)) => {
                lo = l;
                hi = h;
            }
            None => return,
        }
    }
}

fn shortsort<T>(v: &mut [T], lo: isize, mut hi: isize, c: &impl Fn(&[T], isize, isize) -> i32) {
    while hi > lo {
        let mut max = lo;
        for p in lo + 1..=hi {
            if c(v, p, max) > 0 {
                max = p;
            }
        }
        v.swap(max as usize, hi as usize);
        hi -= 1;
    }
}

fn sign_eps(v: f64, eps: f64) -> i32 {
    if eps < v {
        1
    } else if -eps <= v {
        0
    } else {
        -1
    }
}

impl Engine {
    pub(crate) fn insert_pnt_pnt_mid_nodes(&mut self) {
        let num_edges = self.edges.len();
        for e in 0..num_edges {
            let edge = self.edges[e];
            if edge.deleted
                || self.nodes[edge.n1].deg2
                || self.nodes[edge.n2].deg2
                || edge.ltype != SiteType::Pnt
                || edge.rtype != SiteType::Pnt
            {
                continue;
            }
            let pi = self.pnts[edge.lft];
            let pj = self.pnts[edge.rgt];
            let dx = pi.x - pj.x;
            let n1 = self.nodes[edge.n1];
            let det1 = (pj.x - n1.x) * (pj.y - pi.y) + (pj.y - n1.y) * dx;
            if !(det1 < 0.0) {
                continue;
            }
            let n2 = self.nodes[edge.n2];
            let det2 = (pj.x - n2.x) * (pj.y - pi.y) + (pj.y - n2.y) * dx;
            if 0.0 < det2 {
                let mx = (pj.x + pi.x) * 0.5;
                let my = (pi.y + pj.y) * 0.5;
                let r = dist(dx, pi.y - pj.y) * 0.5;
                let node = self.store_node(mx, my, r);
                self.split_edge(e, node, false);
            }
        }
        self.deg2_nodes_added = true;
    }

    /// Tests whether the projection of (x, y) onto the bisector direction of edge `e` lies between
    /// the projections of its two nodes.
    pub(crate) fn is_pnt_between_nodes(
        &self,
        e: usize,
        x: f64,
        y: f64,
        strict: bool,
        eps: f64,
    ) -> bool {
        let edge = self.edges[e];
        let n1 = self.nodes[edge.n1];
        let n2 = self.nodes[edge.n2];
        let (v1, v2) = match (edge.ltype, edge.rtype) {
            (SiteType::Pnt, SiteType::Pnt) => {
                let (lo, hi) = if edge.rgt < edge.lft {
                    (edge.rgt, edge.lft)
                } else {
                    (edge.lft, edge.rgt)
                };
                let q = self.pnts[hi];
                let dv = q.y - self.pnts[lo].y;
                let du = self.pnts[lo].x - q.x;
                let sp = (q.x - x) * dv + (q.y - y) * du;
                let v1 = ((q.x - n1.x) * dv + (q.y - n1.y) * du) - sp;
                let v2 = ((q.y - n2.y) * du + (q.x - n2.x) * dv) - sp;
                (v1, v2)
            }
            (SiteType::Pnt, SiteType::Seg) => {
                self.seg_pnt_between(edge.rgt, edge.lft, &n1, &n2, x, y)
            }
            (SiteType::Seg, SiteType::Pnt) => {
                self.seg_pnt_between(edge.lft, edge.rgt, &n1, &n2, x, y)
            }
            (SiteType::Seg, SiteType::Seg) => {
                let s = if edge.rgt < edge.lft {
                    edge.rgt
                } else {
                    edge.lft
                };
                self.seg_dir_between(s, &n1, &n2, x, y)
            }
            _ => unreachable!("arcs are approximated by segments before Voronoi insertion"),
        };
        let s1 = if eps < v1 {
            1
        } else if -eps < v1 {
            0
        } else {
            -1
        };
        let s2 = if v2 > eps || v2.is_nan() || eps.is_nan() {
            1
        } else if -eps <= v2 {
            0
        } else {
            -1
        };
        let prod = s2 * s1;
        if strict { prod < 0 } else { prod < 1 }
    }

    fn seg_dir_between(&self, s: usize, n1: &Node, n2: &Node, x: f64, y: f64) -> (f64, f64) {
        let seg = self.segs[s];
        let b = seg.b;
        let na = -seg.a;
        let p = self.pnts[seg.i1];
        let sp = (x - p.x) * b + (y - p.y) * na;
        let v1 = ((n1.x - p.x) * b + (n1.y - p.y) * na) - sp;
        let v2 = ((n2.x - p.x) * b + (n2.y - p.y) * na) - sp;
        (v1, v2)
    }

    fn seg_pnt_between(
        &self,
        s: usize,
        pnt: usize,
        n1: &Node,
        n2: &Node,
        x: f64,
        y: f64,
    ) -> (f64, f64) {
        let seg = self.segs[s];
        if seg.i1 != pnt && seg.i2 != pnt {
            return self.seg_dir_between(s, n1, n2, x, y);
        }
        let a = self.pnts[seg.i1];
        let b = self.pnts[seg.i2];
        let dy = b.y - a.y;
        let dx = a.x - b.x;
        let sp = (b.x - x) * dy + (b.y - y) * dx;
        let v1 = if n1.r == 0.0 {
            -sp
        } else {
            ((b.x - n1.x) * dy + (b.y - n1.y) * dx) - sp
        };
        let v2 = if n2.r == 0.0 {
            -sp
        } else {
            ((b.y - n2.y) * dx + (b.x - n2.x) * dy) - sp
        };
        (v1, v2)
    }

    pub(crate) fn check_pnt_on_seg(&mut self, seg: usize, pnt: usize) -> bool {
        let s = self.segs[seg];
        if s.i1 == pnt || s.i2 == pnt {
            return false;
        }
        let p = self.pnts[pnt];
        let d = p.x * s.a + p.y * s.b + s.c;
        if TINY < d || d < -TINY {
            return false;
        }
        let start = self.pnts[s.i1];
        let t = (p.x - start.x) * s.b - (p.y - start.y) * s.a;
        if s.lgth + ZERO < t || t < -ZERO {
            return false;
        }
        let d_start = dist(start.x - p.x, start.y - p.y);
        if !(TINY < d_start || d_start < -TINY) {
            self.merge_pnts(s.i1, pnt);
            return true;
        }
        let end = self.pnts[s.i2];
        let d_end = dist(end.x - p.x, end.y - p.y);
        if !(TINY < d_end || d_end < -TINY) {
            self.merge_pnts(s.i2, pnt);
            return true;
        }
        self.split_seg(seg, pnt);
        true
    }

    pub(crate) fn check_seg_seg_intersection(&mut self, i1: usize, i2: usize) -> bool {
        if i1 == i2 {
            return false;
        }
        let (lo, mut hi) = if i2 < i1 { (i2, i1) } else { (i1, i2) };
        let s1 = self.segs[lo];
        let s2 = self.segs[hi];
        let (st1, en1) = (self.pnts[s1.i1], self.pnts[s1.i2]);
        let (st2, en2) = (self.pnts[s2.i1], self.pnts[s2.i2]);

        let mut shared = 0;
        let sign_s = if s1.i1 == s2.i1 || s1.i2 == s2.i1 {
            shared += 1;
            0
        } else {
            sign_eps(st2.x * s1.a + st2.y * s1.b + s1.c, SMALL)
        };
        let sign_e = if s1.i1 == s2.i2 || s1.i2 == s2.i2 {
            shared += 1;
            0
        } else {
            sign_eps(en2.x * s1.a + en2.y * s1.b + s1.c, SMALL)
        };

        if shared == 2 {
            if s2.apx_arc != NIL {
                if s1.apx_arc == NIL {
                    hi = lo;
                } else {
                    self.pnts[s2.i1].apx_arc = NIL;
                    self.pnts[s2.i2].apx_arc = NIL;
                    self.segs[lo].apx_arc = NIL;
                }
            }
            self.delete_seg(hi);
            return true;
        }
        if sign_e * sign_s != -1 {
            return false;
        }

        let d1 = st1.x * s2.a + st1.y * s2.b + s2.c;
        let side1 = if SMALL < d1 {
            1
        } else if -SMALL < d1 {
            0
        } else {
            -1
        };
        let d2 = en1.x * s2.a + en1.y * s2.b + s2.c;
        let side2 = if d2 > SMALL || d2.is_nan() {
            1
        } else if -SMALL <= d2 {
            0
        } else {
            -1
        };
        if side2 * side1 != -1 {
            return false;
        }

        let dir = Point { x: s1.b, y: -s1.a };
        let t = (d1.abs() * s1.lgth) / (d2.abs() + d1.abs());
        let mut q = Point {
            x: t * dir.x + st1.x,
            y: t * dir.y + st1.y,
        };
        let check = q.x * s2.a + q.y * s2.b + s2.c;
        if SMALL < check || check < -SMALL {
            let mut t = s1.lgth * 0.5;
            if let Some(root) = self.find_root_on_curve(
                RootTarget {
                    a: s2.a,
                    b: s2.b,
                    c: s2.c,
                },
                RootCurve {
                    base: Point { x: st1.x, y: st1.y },
                    dir,
                },
                0.0,
                s1.lgth,
                &mut t,
            ) {
                q = root;
            }
        }
        let new_pnt = self.store_pnt(q.x, q.y);
        self.split_seg(lo, new_pnt);
        self.split_seg(hi, new_pnt);
        true
    }

    pub(crate) fn check_pnt_pnt_duplicate(&mut self, i1: usize, i2: usize) -> bool {
        if i1 == i2 {
            return false;
        }
        let (lo, hi) = if i2 < i1 { (i2, i1) } else { (i1, i2) };
        let d = dist(
            self.pnts[lo].x - self.pnts[hi].x,
            self.pnts[lo].y - self.pnts[hi].y,
        );
        if (-TINY..=TINY).contains(&d) {
            self.merge_pnts(lo, hi);
            return true;
        }
        false
    }

    /// Looks for an input defect involving site `i` (duplicate point, point on segment, crossing
    /// segments) and repairs the first one found.
    pub(crate) fn check_site_intersections(&mut self, i: usize, t: SiteType) -> bool {
        match t {
            SiteType::Pnt => (0..self.pnts.len())
                .any(|j| i != j && !self.pnts[j].deleted && self.check_pnt_pnt_duplicate(i, j)),
            SiteType::Seg => {
                (0..self.pnts.len()).any(|j| !self.pnts[j].deleted && self.check_pnt_on_seg(i, j))
                    || (0..self.segs.len()).any(|j| i != j && self.check_seg_seg_intersection(i, j))
            }
            SiteType::Arc | SiteType::NoSite => {
                unreachable!("arcs are approximated by segments before Voronoi insertion")
            }
        }
    }

    pub(crate) fn insert_seg_endpoint_nodes(&mut self) {
        let mut i = 0;
        while i < self.segs.len() {
            let seg = self.segs[i];
            let p_start = self.pnts[seg.i1];
            let p_end = self.pnts[seg.i2];
            let e0 = seg.vd_edge;
            debug_assert!(
                self.is_lft_site(e0, i, SiteType::Seg) || self.is_rgt_site(e0, i, SiteType::Seg),
                "IsLftSite(e, i, SEG) || IsRgtSite(e, i, SEG)"
            );
            let mut k = if self.is_rgt_site(e0, i, SiteType::Seg) {
                self.edges[e0].n2
            } else {
                self.edges[e0].n1
            };
            let mut e = e0;
            loop {
                let edge = self.edges[e];
                let (other, other_type) = if edge.lft == i && edge.ltype == SiteType::Seg {
                    (edge.rgt, edge.rtype)
                } else {
                    debug_assert!(
                        edge.rgt == i && edge.rtype == SiteType::Seg,
                        "(edges[e].rgt == i) && (edges[e].rtype == SEG)"
                    );
                    (edge.lft, edge.ltype)
                };
                if other_type == SiteType::Pnt && (seg.i1 == other || seg.i2 == other) {
                    let k1 = self.other_node(e, k);
                    let nk = self.nodes[k];
                    let nk1 = self.nodes[k1];
                    if 0.0 < nk.r
                        && 0.0 < nk1.r
                        && (nk.x * seg.a + nk.y * seg.b + seg.c)
                            * (nk1.x * seg.a + nk1.y * seg.b + seg.c)
                            < 0.0
                    {
                        let p = if seg.i1 == other { p_start } else { p_end };
                        let node = self.store_node(p.x, p.y, 0.0);
                        self.split_edge(e, node, false);
                        self.nodes[node].status = NodeStatus::Unchecked;
                    }
                }
                debug_assert!(
                    self.edges[e].n1 == k || self.edges[e].n2 == k,
                    "(edges[e].n1 == k) || (edges[e].n2 == k)"
                );
                e = self.ccw_edge(e, k);
                k = self.other_node(e, k);
                if e == e0 {
                    break;
                }
            }
            i += 1;
        }
    }
}
