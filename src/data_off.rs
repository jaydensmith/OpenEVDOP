use crate::consts::*;
use crate::real::Real;
use crate::state::*;

const DEGENERATE_SUBDIVISIONS: f64 = 10.0;
const DEGENERATE_STEP: f64 = 0.1;
const DEGENERATE_MAX_BUCKET: i32 = 10;
const FIRST_REGULAR_EDGE: usize = 4;

impl<R: Real> Engine<R> {
    pub(crate) fn get_num_edges(&self) -> usize {
        self.edges.len()
    }

    fn is_node_on_seg_side(&self, seg: usize, node: usize, left: bool) -> bool {
        let n = &self.nodes[node];
        let v = self.segs[seg].signed_dist::<R>(n.x, n.y).stored();
        if left { v > 0.0 || v.is_nan() } else { v < 0.0 }
    }

    fn is_offset_edge(&self, e: usize) -> bool {
        let edge = &self.edges[e];
        let inside = |site: usize| self.min_bbox_pnt < site && site < self.max_bbox_pnt;
        !edge.deleted
            && edge.lft != NIL
            && edge.rgt != NIL
            && (edge.ltype != SiteType::Pnt || inside(edge.lft))
            && (edge.rtype != SiteType::Pnt || inside(edge.rgt))
    }

    fn scan_offset_edges(&mut self, start_e: usize, start_n: usize, on_side: bool) {
        let mut stack: Vec<(usize, usize)> = Vec::new();
        let mut e = start_e;
        let mut n = start_n;
        loop {
            if !self.edge_flags[e].available {
                return;
            }
            self.edge_flags[e].available = false;
            n = self.other_node(e, n);
            let r = self.nodes[n].r;
            if on_side && self.is_offset_edge(e) {
                self.active_edges.push(e);
            }
            if r <= ZERO {
                while !self.edge_flags[e].available {
                    let Some((pe, pn)) = stack.pop() else { break };
                    e = pe;
                    n = pn;
                }
            } else {
                let e1 = self.ccw_edge(e, n);
                if !self.edge_flags[e1].available {
                    e = self.cw_edge(e, n);
                    while !self.edge_flags[e].available {
                        let Some((pe, pn)) = stack.pop() else { break };
                        e = pe;
                        n = pn;
                    }
                } else {
                    let e2 = self.cw_edge(e, n);
                    if self.edge_flags[e2].available && e1 != e2 {
                        stack.push((e2, n));
                    }
                    e = e1;
                }
            }
        }
    }

    pub(crate) fn initialize_edge_data(&mut self, both_sides: bool, left: bool) {
        let num_edges = self.get_num_edges();
        let vd_id = self.vd_id;
        if num_edges < self.e_data.len() {
            if self.e_data_vd_id != vd_id {
                self.e_data_vd_id = vd_id;
                self.e_data[..num_edges]
                    .iter_mut()
                    .for_each(|d| d.pending = true);
            }
        } else {
            let fresh = EdgeData {
                pending: true,
                px: 0.0,
                py: 0.0,
                d: 0.0,
                ccw: false,
            };
            self.e_data.resize(num_edges + 1, fresh);
            self.e_data[..num_edges]
                .iter_mut()
                .for_each(|d| d.pending = true);
            self.e_data_vd_id = vd_id;
        }
        if self.edge_flags.len() <= num_edges {
            self.edge_flags.resize(num_edges + 1, EdgeFlags::default());
        }
        self.active_edges.clear();
        self.max_bbox_pnt = self.pnts.len() - 2;
        self.min_bbox_pnt = 1;

        if both_sides {
            for e in FIRST_REGULAR_EDGE..num_edges {
                self.edge_flags[e].available = true;
                if self.is_offset_edge(e) {
                    self.active_edges.push(e);
                }
            }
            self.clear_bbox_edge_flags();
        } else {
            for e in FIRST_REGULAR_EDGE..num_edges {
                self.edge_flags[e].available = true;
            }
            self.clear_bbox_edge_flags();
            for e in FIRST_REGULAR_EDGE..num_edges {
                if !self.edge_flags[e].available || !self.is_offset_edge(e) {
                    continue;
                }
                let edge = self.edges[e];
                let (site, stype) = if edge.ltype == SiteType::Pnt {
                    (edge.rgt, edge.rtype)
                } else {
                    (edge.lft, edge.ltype)
                };
                if stype != SiteType::Seg {
                    continue;
                }
                let r1 = self.nodes[edge.n1].r;
                let r2 = self.nodes[edge.n2].r;
                let n = if r2 > r1 || r1.is_nan() || r2.is_nan() {
                    if r2 <= ZERO { NIL } else { edge.n2 }
                } else if r1 <= ZERO {
                    NIL
                } else {
                    edge.n1
                };
                if n == NIL {
                    continue;
                }
                let on_side = self.is_node_on_seg_side(site, n, left);
                self.scan_offset_edges(e, n, on_side);
                let e1 = self.ccw_edge(e, n);
                if e1 != NIL && e1 >= FIRST_REGULAR_EDGE {
                    self.scan_offset_edges(e1, n, on_side);
                }
                let e2 = self.cw_edge(e, n);
                if e2 != NIL && e2 >= FIRST_REGULAR_EDGE && e1 != e2 {
                    self.scan_offset_edges(e2, n, on_side);
                }
            }
            for e in FIRST_REGULAR_EDGE..num_edges {
                self.edge_flags[e].available = true;
            }
        }
        self.new_offset_needed = false;
    }

    fn clear_bbox_edge_flags(&mut self) {
        for e in 0..FIRST_REGULAR_EDGE {
            self.edge_flags[e].available = false;
        }
    }

    fn edge_radii_degenerate(&self, e: usize) -> bool {
        let edge = &self.edges[e];
        within_zero(R::of(self.nodes[edge.n1].r) - self.nodes[edge.n2].r)
    }

    fn mark_edge_degenerate(&mut self, e: usize) {
        self.edge_flags[e].degenerate = true;
        self.e_data[e].pending = false;
    }

    pub(crate) fn create_parabola_data(&mut self, e: usize) {
        debug_assert!(!self.new_offset_needed, "!e_data_init");
        let edge = self.edges[e];
        let (pnt, seg) = if edge.ltype == SiteType::Seg {
            debug_assert!(edge.rtype == SiteType::Pnt, "type == PNT");
            (edge.rgt, edge.lft)
        } else {
            debug_assert!(edge.ltype == SiteType::Pnt, "type == PNT");
            debug_assert!(edge.rtype == SiteType::Seg, "type == SEG");
            (edge.lft, edge.rgt)
        };
        if self.edge_radii_degenerate(e) {
            self.mark_edge_degenerate(e);
            return;
        }
        let n1 = self.nodes[edge.n1];
        let n2 = self.nodes[edge.n2];
        let p = self.pnts[pnt];
        let s = self.segs[seg];
        let mut a = s.a;
        let mut b = s.b;
        let mut d = (R::of(p.x) * a + R::of(p.y) * b + s.c).stored();
        let far = if n2.r < n1.r { n1 } else { n2 };
        let mut du = d;
        if within_zero(du) {
            du = s.signed_dist::<R>(far.x, far.y).stored();
            debug_assert!(ZERO < du || du < -ZERO, "!eq(du, ZERO)");
        }
        let flip = du < 0.0;
        if flip {
            a = -a;
            b = -b;
            d = -d;
        }
        let ad = (R::of(a) * d).stored();
        let fx = (R::of(p.x) - ad).stored();
        let bd = (R::of(b) * d).stored();
        let fy = (R::of(p.y) - bd).stored();
        let na = -a;
        let k = (-(R::of(p.x) * b) - R::of(na) * p.y).stored();
        let sgn_far = sign_of(R::of(b) * far.x + R::of(na) * far.y + k);

        let disc = ((R::of(d) + d - d) * d).stored();
        let x = if disc < 0.0 {
            (R::of(ad) * 0.5 + fx).stored()
        } else {
            (R::of(ad) + fx + R::of(f64::from(sgn_far)) * b * R::of(disc).sqrt()).stored()
        };
        let y = if disc < 0.0 {
            R::of(bd) * 0.5 + fy
        } else {
            R::of((R::of(bd) + fy).stored()) - R::of(f64::from(sgn_far)) * a * R::of(disc).sqrt()
        };
        let mut sgn = sign_of(R::of(x) * b + y * na + k);
        if sgn != sgn_far {
            sgn = -sgn;
        }

        self.e_data[e].pending = false;
        self.edge_flags[e].degenerate = false;
        let data = &mut self.e_data[e];
        data.px = fx;
        data.py = fy;
        data.d = if flip { -d } else { d };
        data.ccw = sgn > 0;
    }

    fn eval_degenerate_edge_pnt(&self, e: usize, t: f64) -> (f64, f64) {
        debug_assert!(!self.e_data[e].pending, "!GetEdgeDataInit(i)");
        let edge = &self.edges[e];
        let lo = self.nodes[edge.n1];
        let hi = self.nodes[edge.n2];
        if t == lo.r {
            return (lo.x, lo.y);
        }
        if t == hi.r {
            return (hi.x, hi.y);
        }
        let (mut x1, mut y1, mut t1, mut t2) = (lo.x, lo.y, lo.r, hi.r);
        if hi.r < lo.r {
            x1 = hi.x;
            y1 = hi.y;
            t1 = hi.r;
            t2 = lo.r;
        }
        debug_assert!(t1 <= t && t < t2, "(t1 <= t) && (t < t2)");
        let step = ((R::of(t2) - t1) * DEGENERATE_STEP).stored();
        let mut k: i32 = 1;
        let mut reach = (R::of(t1) + step).stored();
        while reach < t && k <= DEGENERATE_MAX_BUCKET {
            k += 1;
            reach = (R::of(f64::from(k)) * step + t1).stored();
        }
        let f = ((R::of(f64::from(k)) - 0.5) / DEGENERATE_SUBDIVISIONS).stored();
        (
            ((R::of(hi.x) - x1) * f + x1).stored(),
            ((R::of(hi.y) - y1) * f + y1).stored(),
        )
    }

    fn debug_assert_t_on_edge(&self, e: usize, t: f64) {
        let edge = &self.edges[e];
        let t1 = self.nodes[edge.n1].r;
        let t2 = self.nodes[edge.n2].r;
        if t2 > t1 {
            debug_assert!(t1 <= t && t <= t2, "(t1 <= t) && (t <= t2)");
        } else {
            debug_assert!(t2 <= t && t <= t1, "(t2 <= t) && (t <= t1)");
        }
    }

    fn edge_midpoint(&self, e: usize) -> (f64, f64) {
        let edge = &self.edges[e];
        let n1 = &self.nodes[edge.n1];
        let n2 = &self.nodes[edge.n2];
        (
            ((R::of(n2.x) + n1.x) * 0.5).stored(),
            ((R::of(n2.y) + n1.y) * 0.5).stored(),
        )
    }

    pub(crate) fn eval_parabola_pnt(&self, e: usize, seg: usize, t: f64) -> (f64, f64) {
        debug_assert!(!self.e_data[e].pending, "!GetEdgeDataInit(i)");
        self.debug_assert_t_on_edge(e, t);
        if self.edge_flags[e].degenerate {
            return self.edge_midpoint(e);
        }
        let s = &self.segs[seg];
        let data = &self.e_data[e];
        let mut a = s.a;
        let mut b = s.b;
        let mut d = data.d;
        let sgn = if data.ccw { 1.0 } else { -1.0 };
        if d < 0.0 {
            a = -a;
            b = -b;
            d = -d;
        }
        let disc = ((R::of(t) + t - d) * d).stored();
        let x = if disc < 0.0 {
            (R::of(d) * a * 0.5 + data.px).stored()
        } else {
            (R::of(a) * t + data.px + R::of(sgn) * b * R::of(disc).sqrt()).stored()
        };
        let y = if disc < 0.0 {
            (R::of(d) * b * 0.5 + data.py).stored()
        } else {
            (R::of((R::of(b) * t + data.py).stored()) - R::of(sgn) * a * R::of(disc).sqrt())
                .stored()
        };
        (x, y)
    }

    pub(crate) fn create_pp_bisector_data(&mut self, e: usize) {
        debug_assert!(!self.new_offset_needed, "!e_data_init");
        let edge = self.edges[e];
        debug_assert!(edge.ltype == SiteType::Pnt, "type == PNT");
        debug_assert!(edge.rtype == SiteType::Pnt, "type == PNT");
        if self.edge_radii_degenerate(e) {
            self.mark_edge_degenerate(e);
            return;
        }
        let n1 = self.nodes[edge.n1];
        let n2 = self.nodes[edge.n2];
        let p1 = self.pnts[edge.lft];
        let p2 = self.pnts[edge.rgt];
        let dx = (R::of(p1.x) - p2.x).stored();
        let dy = (R::of(p1.y) - p2.y).stored();
        let half_dist_sq = ((R::of(dx) * dx + R::of(dy) * dy) * 0.25).stored();
        let (vx, vy) = node_direction::<R>(&n1, &n2);
        let lgth = (R::of(vy) * vy + R::of(vx) * vx).sqrt().stored();
        if within_zero(lgth) {
            self.mark_edge_degenerate(e);
            return;
        }
        self.e_data[e].pending = false;
        self.edge_flags[e].degenerate = false;
        let data = &mut self.e_data[e];
        data.px = (R::of(vx) / lgth).stored();
        data.py = (R::of(vy) / lgth).stored();
        data.d = half_dist_sq;
    }

    pub(crate) fn eval_pp_bisector_pnt(&self, e: usize, t: f64) -> (f64, f64) {
        debug_assert!(!self.e_data[e].pending, "!GetEdgeDataInit(i)");
        self.debug_assert_t_on_edge(e, t);
        if self.edge_flags[e].degenerate {
            return self.edge_midpoint(e);
        }
        let edge = &self.edges[e];
        let n1 = &self.nodes[edge.n1];
        let n2 = &self.nodes[edge.n2];
        let lo = if n1.r > n2.r { n2 } else { n1 };
        let data = &self.e_data[e];
        let mut s0 = (R::of(lo.r) * lo.r - data.d).stored();
        if s0 < 0.0 {
            s0 = 0.0;
        }
        let s1 = R::of(t) * t - data.d;
        let s1 = if s1 < 0.0 { 0.0 } else { s1.stored() };
        let delta = R::of(R::of(s1).sqrt().stored()) - R::of(s0).sqrt();
        (
            (R::of(data.px) * delta + lo.x).stored(),
            (delta * data.py + lo.y).stored(),
        )
    }

    pub(crate) fn create_ss_bisector_data(&mut self, e: usize) {
        debug_assert!(!self.new_offset_needed, "!e_data_init");
        if self.edge_radii_degenerate(e) {
            self.mark_edge_degenerate(e);
            return;
        }
        let edge = self.edges[e];
        let n1 = self.nodes[edge.n1];
        let n2 = self.nodes[edge.n2];
        let (vx, vy) = node_direction::<R>(&n1, &n2);
        self.e_data[e].pending = false;
        self.edge_flags[e].degenerate = false;
        self.e_data[e].px = vx;
        self.e_data[e].py = vy;
    }

    pub(crate) fn eval_ss_bisector_pnt(&self, e: usize, t: f64) -> (f64, f64) {
        debug_assert!(!self.e_data[e].pending, "!GetEdgeDataInit(i)");
        self.debug_assert_t_on_edge(e, t);
        if self.edge_flags[e].degenerate {
            return self.eval_degenerate_edge_pnt(e, t);
        }
        let edge = &self.edges[e];
        let n1 = &self.nodes[edge.n1];
        let n2 = &self.nodes[edge.n2];
        let diff = (R::of(n1.r) - n2.r).stored();
        debug_assert!(ZERO < diff || diff < -ZERO, "!eq(t1-t2, ZERO)");
        let (f, base) = if n2.r < n1.r {
            ((R::of(t) - n2.r) / diff, n2)
        } else {
            ((R::of(t) - n1.r) / (R::of(n2.r) - n1.r), n1)
        };
        let data = &self.e_data[e];
        (
            (R::of(data.px) * f + base.x).stored(),
            (f * data.py + base.y).stored(),
        )
    }
}

/// The vector from the node of smaller clearance to the other.
fn node_direction<R: Real>(n1: &Node, n2: &Node) -> (f64, f64) {
    if n2.r < n1.r {
        ((R::of(n1.x) - n2.x).stored(), (R::of(n1.y) - n2.y).stored())
    } else {
        ((R::of(n2.x) - n1.x).stored(), (R::of(n2.y) - n1.y).stored())
    }
}

fn within_zero<R: Real>(v: R) -> bool {
    (R::of(-ZERO)..=R::of(ZERO)).contains(&v)
}

fn sign_of<R: Real>(v: R) -> i32 {
    if v > 0.0 {
        1
    } else if v < 0.0 {
        -1
    } else {
        0
    }
}
