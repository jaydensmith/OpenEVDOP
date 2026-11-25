use crate::consts::*;
use crate::real::Real;
use crate::state::*;

/// (e_l, e_r, n_l, n_r): outermost new edges and nodes created along the conflict boundary.
type NewEdges = (usize, usize, usize, usize);

const NO_NEW_EDGES: NewEdges = (NIL, NIL, NIL, NIL);

fn sign_open_lower<R: Real>(v: R) -> i32 {
    if v > ZERO {
        1
    } else if v > -ZERO {
        0
    } else {
        -1
    }
}

fn sign_closed_lower<R: Real>(v: R) -> i32 {
    if v > ZERO || v != v {
        1
    } else if v >= -ZERO {
        0
    } else {
        -1
    }
}

impl<R: Real> Engine<R> {
    pub(crate) fn reset_node_status_rec(&mut self, mut e: usize, mut n: usize) {
        loop {
            let k = self.other_node(e, n);
            let status = self.nodes[k].status;
            if status != NodeStatus::Visited && status != NodeStatus::Processed {
                return;
            }
            self.nodes[k].status = NodeStatus::Unchecked;
            if self.nodes[k].deg2 {
                let e1 = self.ccw_edge(e, k);
                debug_assert!(self.ccw_edge(e1, k) == e, "GetCCWEdge(e1, k) == e");
                e = e1;
            } else {
                let e1 = self.ccw_edge(e, k);
                let e2 = self.cw_edge(e, k);
                self.reset_node_status_rec(e1, k);
                e = e2;
            }
            n = k;
        }
    }

    pub(crate) fn mark_loop_nodes_rec(&mut self, mut e: usize, mut n: usize) {
        loop {
            let k = self.other_node(e, n);
            if self.nodes[k].status != NodeStatus::Deleted {
                return;
            }
            self.nodes[k].status = NodeStatus::Processed;
            if self.nodes[k].deg2 {
                let e1 = self.ccw_edge(e, k);
                debug_assert!(self.ccw_edge(e1, k) == e, "GetCCWEdge(e1, k) == e");
                e = e1;
            } else {
                let e1 = self.ccw_edge(e, k);
                let e2 = self.cw_edge(e, k);
                self.mark_loop_nodes_rec(e1, k);
                e = e2;
            }
            n = k;
        }
    }

    fn finish_deg2_node_at_pnt(&mut self, e: usize, pnt: usize, n: usize) {
        self.nodes[n].status = NodeStatus::Unchecked;
        self.nodes[n].site = pnt;
        let p = &mut self.pnts[pnt];
        p.deg2_node = n;
        p.vd_edge = e;
        p.has_deg2 = true;
    }

    pub(crate) fn insert_deg2_node_seg(&mut self, e: usize, seg: usize, pnt: usize) {
        let s = self.segs[seg];
        if s.i1 == pnt || s.i2 == pnt {
            if self.pnts[pnt].has_deg2 {
                return;
            }
            let p = self.pnts[pnt];
            let n = self.store_node(p.x, p.y, 0.0);
            self.split_edge(e, n, false);
            self.finish_deg2_node_at_pnt(e, pnt, n);
            return;
        }

        let edge = self.edges[e];
        let n1 = self.nodes[edge.n1];
        let n2 = self.nodes[edge.n2];
        let p = self.pnts[pnt];
        let u = self.pnts[s.i1];
        let dir_x = s.b;
        let dir_y = -s.a;
        let t0 = (R::of(p.x) - u.x) * dir_x + (R::of(p.y) - u.y) * dir_y;
        let s1 = (R::of(n1.x) - u.x) * dir_x + (R::of(n1.y) - u.y) * dir_y - t0;
        let s2 = (R::of(n2.x) - u.x) * dir_x + (R::of(n2.y) - u.y) * dir_y - t0;
        let prod = sign_open_lower(s1) * sign_closed_lower(s2);
        if prod >= 0 && (prod != 0 || n1.r <= ZERO || n2.r <= ZERO) {
            return;
        }
        let h = (R::of(s.a) * p.x + R::of(p.y) * s.b + s.c) * 0.5;
        if prod >= 0 && (h > ZERO || h < -ZERO) {
            return;
        }
        let r = if h < 0.0 { -h } else { h }.stored();
        let n = self.store_node(
            (R::of(p.x) - h * s.a).stored(),
            (R::of(p.y) - h * s.b).stored(),
            r,
        );
        self.split_edge(e, n, false);
    }

    pub(crate) fn set_loop_stack_start(&mut self, n_prev: usize, n: usize) {
        let found = self
            .loop_stack
            .iter()
            .position(|&e| self.edges[e].n1 == n || self.edges[e].n2 == n);
        debug_assert!(found.is_some(), "1 == 0");
        if let Some(j) = found {
            self.loop_stack_start = if n_prev != n { j + 1 } else { j };
        }
    }

    pub(crate) fn trace_loop_rec(&mut self, mut e: usize, n: usize, end_node: &mut usize) -> bool {
        let saved = self.loop_stack.len();
        let mut k = self.other_node(e, n);
        self.loop_stack.push(e);
        while self.nodes[k].deg2 && self.nodes[k].status == NodeStatus::Processed {
            self.nodes[k].status = NodeStatus::Deleted;
            let e1 = self.ccw_edge(e, k);
            debug_assert!(self.ccw_edge(e1, k) == e, "GetCCWEdge(e1, k) == e");
            let k1 = self.other_node(e1, k);
            self.loop_stack.push(e1);
            k = k1;
            e = e1;
        }
        match self.nodes[k].status {
            NodeStatus::Deleted => {
                *end_node = k;
                return true;
            }
            NodeStatus::Processed => {}
            _ => {
                self.loop_stack.truncate(saved);
                return false;
            }
        }
        self.nodes[k].status = NodeStatus::Deleted;
        let e1 = self.ccw_edge(e, k);
        if self.trace_loop_rec(e1, k, end_node) {
            return true;
        }
        let e2 = self.cw_edge(e, k);
        if self.trace_loop_rec(e2, k, end_node) {
            return true;
        }
        self.loop_stack.truncate(saved);
        false
    }

    /// The candidate is meaningful even when the bool is false: the solvers record a coincident point on rejected paths.
    fn compute_vd_node(
        &mut self,
        e: usize,
        site: usize,
        site_type: SiteType,
    ) -> (bool, NodeCandidate) {
        let mut out = NodeCandidate::new();
        if site_type != SiteType::Seg {
            return (false, out);
        }
        let edge = self.edges[e];
        let (lft, rgt) = (edge.lft, edge.rgt);
        let ok = match (edge.ltype, edge.rtype) {
            (SiteType::Pnt, SiteType::Pnt) => self.seg_pnt_pnt_cntr(site, lft, rgt, e, &mut out),
            (SiteType::Pnt, SiteType::Seg) => self.seg_seg_pnt_cntr(rgt, site, lft, e, &mut out),
            (SiteType::Seg, SiteType::Pnt) => self.seg_seg_pnt_cntr(lft, site, rgt, e, &mut out),
            (SiteType::Seg, SiteType::Seg) => self.seg_seg_seg_cntr(lft, rgt, site, e, &mut out),
            _ => false,
        };
        (ok, out)
    }

    pub(crate) fn insert_new_edges_rec(
        &mut self,
        i: usize,
        t: SiteType,
        n: usize,
        e: usize,
    ) -> NewEdges {
        debug_assert!(
            t == SiteType::Seg || t == SiteType::Arc,
            "(t == SEG) || (t == ARC)"
        );
        let k = self.other_node(e, n);
        if self.nodes[k].status == NodeStatus::Deleted {
            let out;
            if !self.nodes[k].deg2 {
                let e1 = self.ccw_edge(e, k);
                let e2 = self.cw_edge(e, k);
                let first = self.insert_new_edges_rec(i, t, k, e1);
                if self.restart {
                    return NO_NEW_EDGES;
                }
                let second = self.insert_new_edges_rec(i, t, k, e2);
                if self.restart {
                    return NO_NEW_EDGES;
                }
                let (first_e_l, e1r, first_n_l, n1r) = first;
                let (e2l, second_e_r, n2l, second_n_r) = second;
                let e1r_edge = self.edges[e1r];
                let (i0, t0) = if e1r_edge.n2 == n1r {
                    (e1r_edge.rgt, e1r_edge.rtype)
                } else {
                    debug_assert!(e1r_edge.n1 == n1r, "GetStartNode(e1r) == n1r");
                    (e1r_edge.lft, e1r_edge.ltype)
                };
                let e0 = self.store_edge(n1r, n2l, i, i0, t, t0);
                self.set_ccw_edge(e1r, n1r, e0);
                self.set_cw_edge(e0, n1r, e1r);
                self.set_cw_edge(e2l, n2l, e0);
                self.set_ccw_edge(e0, n2l, e2l);
                let dummy = self.cw_edge(e1r, n1r);
                if dummy != NIL {
                    self.set_cw_edge(dummy, n1r, e0);
                    self.set_ccw_edge(e0, n1r, dummy);
                }
                let dummy = self.ccw_edge(e2l, n2l);
                if dummy != NIL {
                    self.set_ccw_edge(dummy, n2l, e0);
                    self.set_cw_edge(e0, n2l, dummy);
                }
                self.set_vd_ptr(i0, t0, e0);
                out = (first_e_l, second_e_r, first_n_l, second_n_r);
            } else {
                let e1 = self.ccw_edge(e, k);
                out = self.insert_new_edges_rec(i, t, k, e1);
                if self.restart {
                    return NO_NEW_EDGES;
                }
            }
            self.del_edge(e);
            self.del_node(k);
            return out;
        }

        if self.nodes[k].status != NodeStatus::Dummy {
            self.nodes[k].status = NodeStatus::Unchecked;
        }
        let (ok, cand) = self.compute_vd_node(e, i, t);
        let (center, r) = if ok {
            if cand.ambiguous && self.check_site_intersections(i, t) {
                return NO_NEW_EDGES;
            }
            (cand.center, cand.r)
        } else {
            if self.check_site_intersections(i, t) {
                return NO_NEW_EDGES;
            }
            match self.desperate_node(e, i, t) {
                Some(node) => node,
                None => return NO_NEW_EDGES,
            }
        };
        let nn = self.store_node(center.x, center.y, r);
        if cand.coincident != NIL {
            self.nodes[nn].site = cand.coincident;
        }
        let edge = &mut self.edges[e];
        if edge.n1 == n {
            edge.n1 = nn;
            edge.s_cw = NIL;
            edge.s_ccw = NIL;
        } else {
            edge.n2 = nn;
            edge.e_cw = NIL;
            edge.e_ccw = NIL;
        }
        self.nodes[nn].edge = e;
        (e, e, nn, nn)
    }
}
