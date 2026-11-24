use crate::consts::*;
use crate::state::*;

impl Engine {
    /// Finds the node of `near_pnt`'s Voronoi cell that the new point at (x, y) violates most.
    pub(crate) fn find_most_violated_node_pnt(&mut self, x: f64, y: f64, near_pnt: usize) -> usize {
        let e0 = self.pnts[near_pnt].vd_edge;
        debug_assert!(
            self.edges[e0].lft == near_pnt || self.edges[e0].rgt == near_pnt,
            "IsLftPnt(e, i) || IsRgtPnt(e, i)"
        );
        let first_node = |engine: &Self| {
            if engine.edges[e0].rgt == near_pnt {
                engine.edges[e0].n2
            } else {
                engine.edges[e0].n1
            }
        };

        let mut best = NIL;
        let mut best_dist = 1.0;
        let mut best_r2 = 0.0;
        let mut k = first_node(self);
        let mut e = e0;
        loop {
            let node = self.nodes[k];
            let dx = x - node.x;
            let dy = y - node.y;
            let dist = dy * dy + dx * dx;
            if dist < node.r && dist * best_r2 < node.r * best_dist {
                best = k;
                best_dist = dist;
                best_r2 = node.r;
            }
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e0 {
                break;
            }
        }
        if best != NIL {
            return best;
        }

        self.troubles = true;
        if self.check_site_intersections(near_pnt, SiteType::Pnt) {
            return NIL;
        }
        let mut best_violation = -f64::MAX;
        let mut k = first_node(self);
        let mut e = e0;
        loop {
            let node = self.nodes[k];
            let dx = x - node.x;
            let dy = y - node.y;
            let violation = node.r.sqrt() - (dy * dy + dx * dx).sqrt();
            if best_violation < violation {
                best_violation = violation;
                best = k;
            }
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e0 {
                break;
            }
        }
        best
    }

    /// Splices a new PNT-PNT edge between boundary ends (e1r, n1r) and (e2l, n2l) of the conflict region.
    fn link_new_pnt_edge(
        &mut self,
        i: usize,
        e1r: usize,
        n1r: usize,
        e2l: usize,
        n2l: usize,
    ) -> usize {
        let i0 = if self.edges[e1r].n2 == n1r {
            self.edges[e1r].rgt
        } else {
            debug_assert!(self.edges[e1r].n1 == n1r, "GetStartNode(e1r) == n1r");
            self.edges[e1r].lft
        };
        let e0 = self.store_edge(n1r, n2l, i, i0, SiteType::Pnt, SiteType::Pnt);
        self.set_ccw_edge(e1r, n1r, e0);
        self.set_cw_edge(e0, n1r, e1r);
        self.set_cw_edge(e2l, n2l, e0);
        self.set_ccw_edge(e0, n2l, e2l);
        let edummy = self.cw_edge(e1r, n1r);
        if edummy != NIL {
            self.set_cw_edge(edummy, n1r, e0);
            self.set_ccw_edge(e0, n1r, edummy);
        }
        let edummy = self.ccw_edge(e2l, n2l);
        if edummy != NIL {
            self.set_ccw_edge(edummy, n2l, e0);
            self.set_cw_edge(e0, n2l, edummy);
        }
        self.pnts[i0].vd_edge = e0;
        e0
    }

    /// Deletes `i`'s conflict region beyond node `n` along `e`; returns its boundary (e_l, e_r, n_l, n_r).
    pub(crate) fn insert_pnt_rec(
        &mut self,
        i: usize,
        n: usize,
        e: usize,
    ) -> (usize, usize, usize, usize) {
        const ABORTED: (usize, usize, usize, usize) = (NIL, NIL, NIL, NIL);
        let k = self.other_node(e, n);
        if self.nodes[k].status == NodeStatus::Unchecked {
            let e1 = self.ccw_edge(e, k);
            let e2 = self.cw_edge(e, k);
            let i1 = if self.edges[e1].n2 == k {
                self.edges[e1].rgt
            } else {
                debug_assert!(k == self.edges[e1].n1, "k == GetStartNode(e1)");
                self.edges[e1].lft
            };
            debug_assert!(
                self.edges[e2].lft == i1 || self.edges[e2].rgt == i1,
                "(IsLftPnt(e2, i1)) || (IsRgtPnt(e2, i1))"
            );
            if i1 != NIL && !self.pnts[i1].visited {
                let k1 = self.other_node(e1, k);
                let k2 = self.other_node(e2, k);
                if self.nodes[k1].status != NodeStatus::Deleted
                    && self.nodes[k2].status != NodeStatus::Deleted
                {
                    let node = self.nodes[k];
                    let dx = self.pnts[i].x - node.x;
                    let dy = self.pnts[i].y - node.y;
                    if dy * dy + dx * dx < node.r {
                        self.pnts[i1].visited = true;
                        self.pnt_visited.push(i1);
                        self.nodes[k].status = NodeStatus::Deleted;
                        debug_assert!(e1 != NIL, "e1 != NIL");
                        debug_assert!(e2 != NIL, "e2 != NIL");
                        let (e_l, e1r, n_l, n1r) = self.insert_pnt_rec(i, k, e1);
                        if self.restart {
                            return ABORTED;
                        }
                        let (e2l, e_r, n2l, n_r) = self.insert_pnt_rec(i, k, e2);
                        if self.restart {
                            return ABORTED;
                        }
                        self.link_new_pnt_edge(i, e1r, n1r, e2l, n2l);
                        self.del_edge(e);
                        self.del_node(k);
                        return (e_l, e_r, n_l, n_r);
                    }
                    self.nodes[k].status = NodeStatus::Visited;
                    self.nodes_checked.push(k);
                }
            }
        }

        let Edge { lft, rgt, .. } = self.edges[e];
        let (ok, center, r2) = self.pnt_pnt_pnt_node(i, lft, rgt);
        if !ok {
            self.troubles = true;
            if self.check_site_intersections(i, SiteType::Pnt)
                || self.check_site_intersections(lft, SiteType::Pnt)
                || self.check_site_intersections(rgt, SiteType::Pnt)
            {
                return ABORTED;
            }
        }
        let nn = self.store_node(center.x, center.y, r2);
        self.nodes[nn].status = NodeStatus::Visited;
        self.nodes_checked.push(nn);
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

    pub(crate) fn insert_pnt_into_vd(&mut self, i: usize) {
        let mut ring_search = false;
        let (mut j, mut dist_sq) = if self.grid_active {
            self.insert_pnt_into_grid(i)
        } else {
            self.find_nearest_pnt(i, false)
        };
        if self.grid_active && dist_sq > ZERO_SQ {
            ring_search = true;
            (j, dist_sq) = self.find_nearest_pnt_in_grid(i, j, dist_sq);
        }

        if dist_sq > ZERO_SQ || dist_sq < -ZERO_SQ {
            let n = self.find_most_violated_node_pnt(self.pnts[i].x, self.pnts[i].y, j);
            if self.restart {
                return;
            }
            let e1 = self.nodes[n].edge;
            let e2 = self.ccw_edge(e1, n);
            let e3 = self.ccw_edge(e2, n);
            self.nodes[n].status = NodeStatus::Deleted;
            let i1 = self.edges[e1].lft;
            let i2 = self.edges[e1].rgt;
            let mut i3 = self.edges[e2].lft;
            if i3 == i1 || i3 == i2 {
                i3 = self.edges[e2].rgt;
            }
            debug_assert!(i3 != i1 && i3 != i2, "(i3 != i1) && (i3 != i2)");
            self.pnt_visited.extend([i1, i2, i3]);
            self.pnts[i1].visited = true;
            self.pnts[i2].visited = true;
            self.pnts[i3].visited = true;

            let (e1l, e1r, n1l, n1r) = self.insert_pnt_rec(i, n, e1);
            if self.restart {
                return;
            }
            let (e2l, e2r, n2l, n2r) = self.insert_pnt_rec(i, n, e2);
            if self.restart {
                return;
            }
            let (e3l, e3r, n3l, n3r) = self.insert_pnt_rec(i, n, e3);
            if self.restart {
                return;
            }
            self.link_new_pnt_edge(i, e1r, n1r, e2l, n2l);
            self.link_new_pnt_edge(i, e2r, n2r, e3l, n3l);
            let e0 = self.link_new_pnt_edge(i, e3r, n3r, e1l, n1l);
            self.del_node(n);
            self.pnts[i].vd_edge = e0;

            for &p in &self.pnt_visited {
                self.pnts[p].visited = false;
            }
            self.pnt_visited.clear();
            for &k in &self.nodes_checked {
                self.nodes[k].status = NodeStatus::Unchecked;
            }
            self.nodes_checked.clear();
            return;
        }

        if ring_search {
            self.delete_pnt_from_grid(i);
        }
        if self.pnts[j].x == self.pnts[i].x && self.pnts[j].y == self.pnts[i].y {
            self.pnts[i].deleted = true;
            if !self.segs.is_empty() || self.num_arcs > 0 {
                self.redirect_pnt(i, j);
            }
            return;
        }
        self.merge_pnts(i, j);
    }
}
