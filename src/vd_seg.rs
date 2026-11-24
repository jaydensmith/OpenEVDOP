use crate::consts::*;
use crate::state::*;

const BETWEEN_EPS_RELAX_FACTOR: f64 = 10.0;

fn eps_sign(v: f64) -> i32 {
    if ZERO < v {
        1
    } else if -ZERO <= v {
        0
    } else {
        -1
    }
}

impl Engine {
    pub(crate) fn get_vd_ptr(&self, i: usize, t: SiteType) -> usize {
        match t {
            SiteType::Pnt => self.pnts[i].vd_edge,
            SiteType::Seg => self.segs[i].vd_edge,
            SiteType::Arc => self.arcs[i].vd_edge,
            SiteType::NoSite => NIL,
        }
    }

    pub(crate) fn is_lft_site(&self, e: usize, i: usize, t: SiteType) -> bool {
        self.edges[e].lft == i && self.edges[e].ltype == t
    }

    pub(crate) fn is_rgt_site(&self, e: usize, i: usize, t: SiteType) -> bool {
        self.edges[e].rgt == i && self.edges[e].rtype == t
    }

    pub(crate) fn is_site_of_edge(&self, e: usize, i: usize, t: SiteType) -> bool {
        self.is_lft_site(e, i, t) || self.is_rgt_site(e, i, t)
    }

    /// Node of `e` from which a CCW walk around the cell of site (i, t) starts.
    fn cell_start_node(&self, e: usize, i: usize, t: SiteType) -> usize {
        if self.is_rgt_site(e, i, t) {
            self.edges[e].n2
        } else {
            self.edges[e].n1
        }
    }

    fn other_site(&self, e: usize, i: usize, t: SiteType) -> (usize, SiteType) {
        let edge = &self.edges[e];
        if edge.lft == i && edge.ltype == t {
            (edge.rgt, edge.rtype)
        } else {
            debug_assert!(
                edge.rgt == i && edge.rtype == t,
                "(edges[e].rgt == i) && (edges[e].rtype == t)"
            );
            (edge.lft, edge.ltype)
        }
    }

    /// |a x + b y + c| when (x, y) projects into seg's extent widened by `eps`, `large` otherwise.
    pub(crate) fn seg_clearance(&self, seg: usize, x: f64, y: f64, eps: f64) -> f64 {
        let s = &self.segs[seg];
        let p = &self.pnts[s.i1];
        let t = (x - p.x) * s.b - (y - p.y) * s.a;
        if s.lgth + eps < t || t < -eps {
            self.large
        } else {
            s.signed_dist(x, y).abs()
        }
    }

    /// False also for NaN clearances.
    fn node_in_seg_conflict(&self, seg: usize, k: usize) -> bool {
        let node = &self.nodes[k];
        self.seg_clearance(seg, node.x, node.y, self.between_eps) - node.r <= self.between_eps
    }

    fn relax_between_eps(&mut self) -> bool {
        if self.between_eps < 0.0 {
            self.between_eps = 0.0;
            return true;
        }
        if self.between_eps == 0.0 {
            self.between_eps = SMALL;
            return true;
        }
        self.between_eps *= BETWEEN_EPS_RELAX_FACTOR;
        self.between_eps <= ZERO_MAX
    }

    fn push_check_pnt(&mut self, site: usize, stype: SiteType) {
        if stype == SiteType::Pnt && !self.pnts[site].has_deg2 {
            self.check_pnts.push(CheckPnt {
                site,
                stype: SiteType::Pnt,
            });
        }
    }

    pub(crate) fn scan_seg_start_node_rec(
        &self,
        mut e: usize,
        k_stop: usize,
        i1: usize,
        seg: usize,
        k_from: usize,
    ) -> (usize, f64) {
        debug_assert!(
            !(self.is_lft_site(e, i1, SiteType::Pnt) && self.is_rgt_site(e, i1, SiteType::Pnt)),
            "!(IsLftSite(e, i1, PNT) && IsRgtSite(e, i1, PNT))"
        );
        let mut k = self.other_node(e, k_from);
        loop {
            if k == k_stop {
                return (NIL, f64::MAX);
            }
            let node = self.nodes[k];
            if node.site != NIL && node.deg2 {
                let e1 = self.ccw_edge(e, k);
                debug_assert!(self.cw_edge(e, k) == e1, "GetCWEdge(e, k1) == e1");
                e = e1;
                k = self.other_node(e, k);
                continue;
            }
            if node.site == NIL {
                return (k, self.seg_clearance(seg, node.x, node.y, SMALL) - node.r);
            }
            let e1 = self.ccw_edge(e, k);
            let e2 = self.ccw_edge(e1, k);
            let (na, da) = self.scan_seg_start_node_rec(e1, k_stop, i1, seg, k);
            let (nb, db) = self.scan_seg_start_node_rec(e2, k_stop, i1, seg, k);
            return if db <= da { (nb, db) } else { (na, da) };
        }
    }

    /// Node of pnt i1's cell most in conflict with seg; the flag reports that it also lies on i2's cell.
    pub(crate) fn find_seg_start_node(&self, i1: usize, i2: usize, seg: usize) -> (usize, bool) {
        let s = self.segs[seg];
        debug_assert!(
            s.i1 == i1 || s.i2 == i1,
            "IsSegStartPnt(j, i1) || IsSegEndPnt(j, i1)"
        );
        debug_assert!(
            s.i1 == i2 || s.i2 == i2,
            "IsSegStartPnt(j, i2) || IsSegEndPnt(j, i2)"
        );
        let e0 = self.pnts[i1].vd_edge;
        debug_assert!(
            self.is_site_of_edge(e0, i1, SiteType::Pnt),
            "IsLftSite(e, i1, PNT) || IsRgtSite(e, i1, PNT)"
        );
        debug_assert!(!(s.lgth <= ZERO && s.lgth >= -ZERO), "!eq(lgth, ZERO)");

        let shared_by_both = |e: usize| {
            self.is_site_of_edge(e, i1, SiteType::Pnt) && self.is_site_of_edge(e, i2, SiteType::Pnt)
        };
        let mut same_node = false;
        let mut best = NIL;
        let mut best_dist = f64::MAX;
        let mut k = self.cell_start_node(e0, i1, SiteType::Pnt);
        let mut prev_shared = shared_by_both(e0);
        let mut e = e0;
        loop {
            e = self.ccw_edge(e, k);
            let cur_shared = shared_by_both(e);
            let node = self.nodes[k];
            if node.site == NIL {
                let dist = self.seg_clearance(seg, node.x, node.y, SMALL) - node.r;
                if dist < best_dist {
                    best = k;
                    best_dist = dist;
                    if prev_shared || cur_shared {
                        same_node = true;
                    }
                }
            }
            k = self.other_node(e, k);
            prev_shared = cur_shared;
            if e == e0 {
                break;
            }
        }

        if (best == NIL || best_dist >= 0.0) && self.pnts[i1].has_deg2 {
            same_node = false;
            let kd = self.pnts[i1].deg2_node;
            debug_assert!(self.nodes[kd].site == i1, "GetNodeSite(k) == i1");
            debug_assert!(self.nodes[kd].deg2, "IsDeg2Node(k)");
            let ea = self.nodes[kd].edge;
            let eb = self.ccw_edge(ea, kd);
            debug_assert!(self.cw_edge(ea, kd) == eb, "GetCWEdge(e, k) == e1");
            let (na, da) = self.scan_seg_start_node_rec(ea, kd, i1, seg, kd);
            let (nb, db) = self.scan_seg_start_node_rec(eb, kd, i1, seg, kd);
            best = if db <= da { nb } else { na };
        }
        (best, same_node)
    }

    pub(crate) fn find_seg_start_node_at_pnt(&mut self, i: usize, t: SiteType, i0: usize) -> bool {
        debug_assert!(t == SiteType::Pnt, "(t == PNT) && InPntsList(i)");
        let e0 = self.get_vd_ptr(i, t);
        debug_assert!(
            self.is_site_of_edge(e0, i, t),
            "IsLftSite(e, i, t) || IsRgtSite(e, i, t)"
        );

        let mut k = self.cell_start_node(e0, i, t);
        let mut e = e0;
        loop {
            if self.nodes[k].status != NodeStatus::Processed {
                return true;
            }
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e0 {
                break;
            }
        }

        let s = self.segs[i0];
        if s.i1 == i || s.i2 == i {
            return self.split_cell_at_seg_endpoint(i, t);
        }

        let p = self.pnts[i];
        let (mut a, mut b) = (s.a, s.b);
        let d = p.x * a + p.y * b + s.c;
        if (-ZERO..=ZERO).contains(&d) {
            if self.check_site_intersections(i0, SiteType::Seg) {
                self.restart = true;
            }
            return false;
        }
        if d <= 0.0 {
            b = -b;
        } else {
            a = -a;
        }
        let c = -(b * p.x) - a * p.y;
        let dir = Point { x: -a, y: b };
        let side = |node: &Node| eps_sign(node.x * b + node.y * a + c);

        let (mut k, k_prev) = if self.is_rgt_site(e0, i, t) {
            (self.edges[e0].n2, self.edges[e0].n1)
        } else {
            (self.edges[e0].n1, self.edges[e0].n2)
        };
        let mut r_prev = self.nodes[k_prev].r;
        let mut side_prev = side(&self.nodes[k_prev]);
        let mut e = e0;
        loop {
            let r_k = self.nodes[k].r;
            let side_k = side(&self.nodes[k]);
            if side_k == -side_prev
                && ZERO < r_k
                && ZERO < r_prev
                && let Some((q, r)) = self.pnt_ray_cntr(i, dir, e)
                && self.is_pnt_between_nodes(e, q.x, q.y, false, self.between_eps)
            {
                let n = self.store_node(q.x, q.y, r);
                self.split_edge(e, n, true);
                self.nodes[n].status = NodeStatus::Dummy;
                return true;
            }
            side_prev = side_k;
            r_prev = r_k;
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e0 {
                return false;
            }
        }
    }

    /// Splits the cell edge ending at the node that coincides with seg endpoint i by a degree-2 node there.
    fn split_cell_at_seg_endpoint(&mut self, i: usize, t: SiteType) -> bool {
        let e0 = self.get_vd_ptr(i, t);
        let mut k = self.cell_start_node(e0, i, t);
        let mut e = e0;
        loop {
            if self.nodes[k].site == i {
                let Node { x, y, r, .. } = self.nodes[k];
                debug_assert!(r == 0.0, "r2 == 0.0");
                let n = self.store_node(x, y, 0.0);
                self.split_edge(e, n, false);
                self.nodes[n].site = i;
                self.nodes[n].status = NodeStatus::Visited;
                let pnt = &mut self.pnts[i];
                pnt.deg2_node = n;
                pnt.vd_edge = e;
                pnt.has_deg2 = true;
                return true;
            }
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e0 {
                return false;
            }
        }
    }

    pub(crate) fn break_loop_at_pnt_pnt_edges(&mut self, i0: usize) -> bool {
        let num_loop = self.loop_stack.len();
        if num_loop <= self.loop_stack_start {
            return true;
        }
        let s = self.segs[i0];
        for j in self.loop_stack_start..num_loop {
            let e = self.loop_stack[j];
            let edge = self.edges[e];
            let site = match edge.ltype {
                SiteType::Seg if edge.rtype == SiteType::Pnt => edge.rgt,
                SiteType::Pnt if edge.rtype == SiteType::Pnt => {
                    let pl = self.pnts[edge.lft];
                    let pr = self.pnts[edge.rgt];
                    let dl = self.seg_clearance(i0, pl.x, pl.y, self.between_eps);
                    let dr = self.seg_clearance(i0, pr.x, pr.y, self.between_eps);
                    if dr < dl { edge.rgt } else { edge.lft }
                }
                SiteType::Pnt => edge.lft,
                _ => continue,
            };

            let q = self.pnts[site];
            let (a, b, c) = if s.i1 == site {
                (s.a, s.b, s.c)
            } else if s.i2 == site {
                (-s.a, -s.b, -s.c)
            } else {
                let d = q.y * s.b + q.x * s.a + s.c;
                if (-ZERO..=ZERO).contains(&d) {
                    if self.check_site_intersections(i0, SiteType::Seg) {
                        self.restart = true;
                    }
                    return false;
                }
                let (a1, b1) = if d < 0.0 { (-s.a, -s.b) } else { (s.a, s.b) };
                let (na, nb) = (b1, -a1);
                (na, nb, -(na * q.x) - q.y * nb)
            };

            let n1 = self.nodes[edge.n1];
            let n2 = self.nodes[edge.n2];
            let side2 = eps_sign(n2.x * a + n2.y * b + c);
            let side1 = eps_sign(n1.x * a + n1.y * b + c);
            if side1 == -side2
                && ZERO < n1.r
                && ZERO < n2.r
                && n1.status == NodeStatus::Deleted
                && n2.status == NodeStatus::Deleted
                && let Some((p, r)) = self.pnt_ray_cntr(site, Point { x: -b, y: a }, e)
                && self.is_pnt_between_nodes(e, p.x, p.y, false, self.between_eps)
            {
                let n = self.store_node(p.x, p.y, r);
                self.split_edge(e, n, true);
                self.nodes[n].status = NodeStatus::Dummy;
                return true;
            }
        }
        false
    }

    pub(crate) fn mark_seg_conflict_tree(&mut self, mut e: usize, mut n: usize, seg: usize) {
        loop {
            debug_assert!(
                self.edges[e].n1 == n || self.edges[e].n2 == n,
                "(edges[e].n1 == n) || (edges[e].n2 == n)"
            );
            let mut k = self.other_node(e, n);
            loop {
                if self.nodes[k].status != NodeStatus::Unchecked {
                    return;
                }
                if !self.nodes[k].deg2 {
                    break;
                }
                let e1 = self.ccw_edge(e, k);
                debug_assert!(self.cw_edge(e, k) == e1, "GetCWEdge(e, k) == e1");
                let k1 = self.other_node(e1, k);
                if self.nodes[k].site != NIL || !self.node_in_seg_conflict(seg, k) {
                    self.nodes[k].status = NodeStatus::Visited;
                    return;
                }
                self.nodes[k].status = NodeStatus::Processed;
                if self.nodes[k1].status == NodeStatus::Processed {
                    let edge1 = self.edges[e1];
                    self.push_check_pnt(edge1.lft, edge1.ltype);
                    self.push_check_pnt(edge1.rgt, edge1.rtype);
                    return;
                }
                if self.nodes[k1].status != NodeStatus::Unchecked {
                    return;
                }
                e = e1;
                k = k1;
            }

            let site = self.nodes[k].site;
            if self.segs[seg].i1 == site
                || self.segs[seg].i2 == site
                || !self.node_in_seg_conflict(seg, k)
            {
                self.nodes[k].status = NodeStatus::Visited;
                return;
            }
            self.nodes[k].status = NodeStatus::Processed;

            let e1 = self.ccw_edge(e, k);
            let e2 = self.cw_edge(e, k);
            let (i1, t1) = if self.edges[e1].n2 == k {
                (self.edges[e1].rgt, self.edges[e1].rtype)
            } else {
                debug_assert!(k == self.edges[e1].n1, "k == GetStartNode(e1)");
                (self.edges[e1].lft, self.edges[e1].ltype)
            };
            debug_assert!(
                self.is_site_of_edge(e2, i1, t1),
                "(IsLftSite(e2, i1, t1)) || (IsRgtSite(e2, i1, t1))"
            );
            let k1 = self.other_node(e1, k);
            let k2 = self.other_node(e2, k);
            if self.nodes[k1].status == NodeStatus::Processed {
                let (s1, ts1) = self.other_site(e1, i1, t1);
                self.push_check_pnt(s1, ts1);
                if self.nodes[k2].status == NodeStatus::Processed {
                    self.push_check_pnt(i1, t1);
                    let (s2, ts2) = self.other_site(e2, i1, t1);
                    self.push_check_pnt(s2, ts2);
                }
            } else if self.nodes[k2].status == NodeStatus::Processed {
                let (s2, ts2) = self.other_site(e2, i1, t1);
                self.push_check_pnt(s2, ts2);
            }
            if self.nodes[k1].status == NodeStatus::Unchecked {
                self.mark_seg_conflict_tree(e1, k, seg);
            }
            if self.nodes[k2].status != NodeStatus::Unchecked {
                return;
            }
            e = e2;
            n = k;
        }
    }

    /// Returns false when the insertion has to be abandoned (restart or input intersections).
    fn break_traced_loop(&mut self, seg: usize, n1: usize, end: usize) -> bool {
        self.set_loop_stack_start(n1, end);
        self.between_eps = ZERO;
        loop {
            if self.break_loop_at_pnt_pnt_edges(seg) {
                break;
            }
            if self.restart {
                return false;
            }
            if self.relax_between_eps() {
                continue;
            }
            if self.check_site_intersections(seg, SiteType::Seg) {
                return false;
            }
            if !self.desperate_split_loop_max_dist(seg) {
                return false;
            }
            break;
        }
        self.loop_stack.clear();
        self.loop_stack_start = 0;
        true
    }

    /// Creates the new edge of seg from `nr` (end of `er`) to `nl` (end of `el`) and links it in at both nodes.
    fn store_seg_edge(&mut self, seg: usize, er: usize, nr: usize, el: usize, nl: usize) -> usize {
        let (i, t) = if self.edges[er].n2 == nr {
            (self.edges[er].rgt, self.edges[er].rtype)
        } else {
            debug_assert!(self.edges[er].n1 == nr, "GetStartNode(e1r) == n1r");
            (self.edges[er].lft, self.edges[er].ltype)
        };
        let e0 = self.store_edge(nr, nl, seg, i, SiteType::Seg, t);
        self.set_vd_ptr(i, t, e0);
        self.set_ccw_edge(er, nr, e0);
        self.set_cw_edge(e0, nr, er);
        self.set_cw_edge(el, nl, e0);
        self.set_ccw_edge(e0, nl, el);
        let edummy = self.cw_edge(er, nr);
        if edummy != NIL {
            self.set_cw_edge(edummy, nr, e0);
            self.set_ccw_edge(e0, nr, edummy);
        }
        let edummy = self.ccw_edge(el, nl);
        if edummy != NIL {
            self.set_ccw_edge(edummy, nl, e0);
            self.set_cw_edge(e0, nl, edummy);
        }
        e0
    }

    pub(crate) fn insert_seg_into_vd(&mut self, seg: usize) {
        self.check_pnts.clear();
        let i1 = self.segs[seg].i1;
        let i2 = self.segs[seg].i2;

        let (mut n1, mut same_node) = self.find_seg_start_node(i1, i2, seg);
        if n1 == NIL {
            if self.check_site_intersections(seg, SiteType::Seg) {
                return;
            }
            n1 = self.desperate_split_pnt_edge(i1, seg);
            same_node = false;
            debug_assert!(n1 != NIL, "n1 != NIL");
        }
        self.nodes[n1].status = NodeStatus::Unchecked;
        let mut n2 = n1;
        if !same_node {
            n2 = self.find_seg_start_node(i2, i1, seg).0;
            if n2 == NIL {
                if self.check_site_intersections(seg, SiteType::Seg) {
                    return;
                }
                n2 = self.desperate_split_pnt_edge(i2, seg);
                debug_assert!(n2 != NIL, "n2 != NIL");
            }
            self.nodes[n2].status = NodeStatus::Unchecked;
        }

        self.between_eps = -SMALL;
        loop {
            let deg2 = self.nodes[n1].deg2;
            let e1 = self.nodes[n1].edge;
            let e2 = self.ccw_edge(e1, n1);
            let e3 = if deg2 { NIL } else { self.ccw_edge(e2, n1) };
            self.nodes[n1].status = NodeStatus::Processed;
            if !deg2 {
                debug_assert!(
                    {
                        let (ed1, ed2) = (self.edges[e1], self.edges[e2]);
                        let matches = |i: usize, t: SiteType| {
                            (i == ed1.lft && t == ed1.ltype) || (i == ed1.rgt && t == ed1.rtype)
                        };
                        let (i3, t3) = if matches(ed2.lft, ed2.ltype) {
                            (ed2.rgt, ed2.rtype)
                        } else {
                            (ed2.lft, ed2.ltype)
                        };
                        !((i3 == ed1.lft && t3 == ed1.ltype) && (i3 == ed1.rgt && t3 == ed1.rtype))
                    },
                    "!(((i3 == i1) && (t3 == t1)) && ((i3 == i2) && (t3 == t2)))"
                );
            }
            self.mark_seg_conflict_tree(e1, n1, seg);
            self.mark_seg_conflict_tree(e2, n1, seg);
            if !deg2 {
                self.mark_seg_conflict_tree(e3, n1, seg);
            }

            if self.nodes[n2].status != NodeStatus::Processed {
                if self.relax_between_eps() {
                    self.check_pnts.clear();
                    self.reset_node_status_rec(e1, n1);
                    self.reset_node_status_rec(e2, n1);
                    if !self.nodes[n1].deg2 {
                        self.reset_node_status_rec(e3, n1);
                    }
                } else {
                    if self.check_site_intersections(seg, SiteType::Seg) {
                        return;
                    }
                    self.nodes[n2].status = NodeStatus::Unchecked;
                    self.desperate_used = true;
                    self.troubles = true;
                    n2 = n1;
                }
            }
            if self.nodes[n2].status == NodeStatus::Processed {
                break;
            }
        }

        let mut j = 0;
        while j < self.check_pnts.len() {
            let cp = self.check_pnts[j];
            self.between_eps = ZERO;
            loop {
                if self.find_seg_start_node_at_pnt(cp.site, cp.stype, seg) {
                    break;
                }
                if self.restart {
                    return;
                }
                if self.relax_between_eps() {
                    continue;
                }
                if self.check_site_intersections(seg, SiteType::Seg) {
                    return;
                }
                self.desperate_split_cell_max_dist(cp.site, cp.stype, seg);
                break;
            }
            j += 1;
        }

        self.check_pnts.clear();
        let mut e1 = self.nodes[n1].edge;
        self.loop_stack.clear();
        self.loop_stack_start = 0;
        self.nodes[n1].status = NodeStatus::Deleted;
        let mut end = NIL;
        while self.trace_loop_rec(e1, n1, &mut end) {
            if !self.break_traced_loop(seg, n1, end) {
                return;
            }
            e1 = self.nodes[n1].edge;
            self.mark_loop_nodes_rec(e1, n1);
            if end == n1 {
                let e2 = self.ccw_edge(e1, n1);
                self.mark_loop_nodes_rec(e2, n1);
                let e3 = self.ccw_edge(e2, n1);
                if e3 != e1 {
                    self.mark_loop_nodes_rec(e3, n1);
                }
            }
            self.nodes[n1].status = NodeStatus::Deleted;
        }

        let e1 = self.nodes[n1].edge;
        let mut e2 = self.ccw_edge(e1, n1);
        self.loop_stack.clear();
        self.loop_stack_start = 0;
        while self.trace_loop_rec(e2, n1, &mut end) {
            if !self.break_traced_loop(seg, n1, end) {
                return;
            }
            e2 = self.ccw_edge(e1, n1);
            self.mark_loop_nodes_rec(e2, n1);
            if end == n1 {
                let e3 = self.ccw_edge(e2, n1);
                if e3 != e1 {
                    self.mark_loop_nodes_rec(e3, n1);
                }
            }
        }

        let e2 = self.ccw_edge(e1, n1);
        let mut e3 = self.ccw_edge(e2, n1);
        self.loop_stack.clear();
        self.loop_stack_start = 0;
        if e3 != e1 {
            while self.trace_loop_rec(e3, n1, &mut end) {
                if !self.break_traced_loop(seg, n1, end) {
                    return;
                }
                e3 = self.ccw_edge(e2, n1);
                self.mark_loop_nodes_rec(e3, n1);
            }
        }

        let e1 = self.nodes[n1].edge;
        let e2 = self.ccw_edge(e1, n1);
        let last_edge = if self.nodes[n1].deg2 {
            let (e1l, e1r, n1l, n1r) = self.insert_new_edges_rec(seg, SiteType::Seg, n1, e1);
            if self.restart {
                return;
            }
            let (e2l, e2r, n2l, n2r) = self.insert_new_edges_rec(seg, SiteType::Seg, n1, e2);
            if self.restart {
                return;
            }
            self.store_seg_edge(seg, e1r, n1r, e2l, n2l);
            self.store_seg_edge(seg, e2r, n2r, e1l, n1l)
        } else {
            let e3 = self.ccw_edge(e2, n1);
            let (e1l, e1r, n1l, n1r) = self.insert_new_edges_rec(seg, SiteType::Seg, n1, e1);
            if self.restart {
                return;
            }
            let (e2l, e2r, n2l, n2r) = self.insert_new_edges_rec(seg, SiteType::Seg, n1, e2);
            if self.restart {
                return;
            }
            let (e3l, e3r, n3l, n3r) = self.insert_new_edges_rec(seg, SiteType::Seg, n1, e3);
            if self.restart {
                return;
            }
            self.store_seg_edge(seg, e1r, n1r, e2l, n2l);
            self.store_seg_edge(seg, e2r, n2r, e3l, n3l);
            self.store_seg_edge(seg, e3r, n3r, e1l, n1l)
        };

        self.del_node(n1);
        self.segs[seg].vd_edge = last_edge;
        self.remove_dummy_nodes();

        let e_start = self.segs[seg].vd_edge;
        debug_assert!(
            self.is_site_of_edge(e_start, seg, SiteType::Seg),
            "IsLftSite(e, i, SEG) || IsRgtSite(e, i, SEG)"
        );
        let mut k = self.cell_start_node(e_start, seg, SiteType::Seg);
        let mut e = e_start;
        loop {
            let (site, stype) = self.other_site(e, seg, SiteType::Seg);
            if stype == SiteType::Pnt {
                self.insert_deg2_node_seg(e, seg, site);
            }
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e_start {
                break;
            }
        }
        debug_assert!(
            self.pnts[self.segs[seg].i1].vd_edge != NIL,
            "GetVDPtr(i1, PNT) != NIL"
        );
        debug_assert!(
            self.pnts[self.segs[seg].i2].vd_edge != NIL,
            "GetVDPtr(i2, PNT) != NIL"
        );
    }
}
