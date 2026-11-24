use crate::consts::*;
use crate::state::*;

impl Engine {
    pub(crate) fn ccw_edge(&self, e: usize, n: usize) -> usize {
        let edge = &self.edges[e];
        debug_assert!(
            edge.n1 == n || edge.n2 == n,
            "(edges[e].n1 == n) || (edges[e].n2 == n)"
        );
        if edge.n1 == n { edge.s_ccw } else { edge.e_ccw }
    }

    pub(crate) fn cw_edge(&self, e: usize, n: usize) -> usize {
        let edge = &self.edges[e];
        debug_assert!(
            edge.n1 == n || edge.n2 == n,
            "(edges[e].n1 == n) || (edges[e].n2 == n)"
        );
        if edge.n1 == n { edge.s_cw } else { edge.e_cw }
    }

    pub(crate) fn set_ccw_edge(&mut self, e: usize, n: usize, value: usize) {
        let edge = &mut self.edges[e];
        debug_assert!(
            edge.n1 == n || edge.n2 == n,
            "(edges[e].n1 == n) || (edges[e].n2 == n)"
        );
        if edge.n1 == n {
            edge.s_ccw = value;
        } else {
            edge.e_ccw = value;
        }
    }

    pub(crate) fn set_cw_edge(&mut self, e: usize, n: usize, value: usize) {
        let edge = &mut self.edges[e];
        debug_assert!(
            edge.n1 == n || edge.n2 == n,
            "(edges[e].n1 == n) || (edges[e].n2 == n)"
        );
        if edge.n1 == n {
            edge.s_cw = value;
        } else {
            edge.e_cw = value;
        }
    }

    pub(crate) fn other_node(&self, e: usize, n: usize) -> usize {
        let edge = &self.edges[e];
        if edge.n1 == n { edge.n2 } else { edge.n1 }
    }

    pub(crate) fn set_vd_ptr(&mut self, site: usize, site_type: SiteType, e: usize) {
        match site_type {
            SiteType::Pnt => self.pnts[site].vd_edge = e,
            SiteType::Seg => self.segs[site].vd_edge = e,
            SiteType::Arc => self.arcs[site].vd_edge = e,
            SiteType::NoSite => {}
        }
    }

    pub(crate) fn init_nodes(&mut self) {
        self.nodes.clear();
        self.free_nodes.clear();
    }

    pub(crate) fn init_edges(&mut self) {
        self.edges.clear();
        self.free_edges.clear();
    }

    pub(crate) fn nodes_r2_to_r(&mut self) {
        for node in &mut self.nodes {
            debug_assert!(node.r >= 0.0, "nodes[i].r2 >= 0.0");
            node.r = node.r.sqrt();
        }
        self.nodes_store_r2 = false;
    }

    pub(crate) fn del_node(&mut self, n: usize) {
        self.free_nodes.push(n);
        self.nodes[n].status = NodeStatus::Deleted;
    }

    pub(crate) fn del_edge(&mut self, e: usize) {
        let edge = &mut self.edges[e];
        edge.lft = NIL;
        edge.rgt = NIL;
        edge.deleted = true;
        self.free_edges.push(e);
    }

    pub(crate) fn store_node(&mut self, x: f64, y: f64, r: f64) -> usize {
        let node = Node {
            x,
            y,
            r,
            status: NodeStatus::Unchecked,
            edge: NIL,
            site: NIL,
            deg2: false,
        };
        match self.free_nodes.pop() {
            Some(i) => {
                self.nodes[i] = node;
                i
            }
            None => {
                self.nodes.push(node);
                self.nodes.len() - 1
            }
        }
    }

    pub(crate) fn store_edge(
        &mut self,
        n1: usize,
        n2: usize,
        lft: usize,
        rgt: usize,
        ltype: SiteType,
        rtype: SiteType,
    ) -> usize {
        let edge = Edge {
            n1,
            n2,
            lft,
            rgt,
            ltype,
            rtype,
            s_ccw: NIL,
            s_cw: NIL,
            e_ccw: NIL,
            e_cw: NIL,
            deleted: false,
        };
        match self.free_edges.pop() {
            Some(i) => {
                self.edges[i] = edge;
                i
            }
            None => {
                self.edges.push(edge);
                self.edges.len() - 1
            }
        }
    }

    /// Builds the initial VD of the four bounding-box corners (pnts 0, 1, N-2, N-1) and pnt 2.
    pub(crate) fn init_vd(&mut self) {
        debug_assert!(self.edges.is_empty(), "num_edges == 0");
        debug_assert!(self.nodes.is_empty(), "num_nodes == 0");
        debug_assert!(self.pnts.len() > 4, "num_pnts > 4");
        self.vd_id += 1;
        self.vd_computed = false;
        self.init_nodes();
        self.init_edges();

        let i = self.pnts.len() - 2;
        let j = self.pnts.len() - 1;
        debug_assert!(
            self.pnts[0].x == self.pnts[1].x,
            "pnts[0].p.x == pnts[1].p.x"
        );
        debug_assert!(self.pnts[0].y < self.pnts[1].y, "pnts[0].p.y < pnts[1].p.y");
        debug_assert!(
            self.pnts[i].x == self.pnts[j].x,
            "pnts[i].p.x == pnts[j].p.x"
        );
        debug_assert!(self.pnts[i].y < self.pnts[j].y, "pnts[i].p.y < pnts[j].p.y");
        debug_assert!(self.pnts[0].x < self.pnts[i].x, "pnts[0].p.x < pnts[i].p.x");
        debug_assert!(
            self.pnts[0].y == self.pnts[i].y,
            "pnts[0].p.y == pnts[i].p.y"
        );
        debug_assert!(self.pnts[1].x < self.pnts[j].x, "pnts[1].p.x < pnts[j].p.x");
        debug_assert!(
            self.pnts[1].y == self.pnts[j].y,
            "pnts[1].p.y == pnts[j].p.y"
        );
        debug_assert!(self.pnts[0].x < self.pnts[2].x, "pnts[0].p.x < pnts[2].p.x");
        debug_assert!(self.pnts[j].x > self.pnts[2].x, "pnts[j].p.x > pnts[2].p.x");
        debug_assert!(self.pnts[0].y < self.pnts[2].y, "pnts[0].p.y < pnts[2].p.y");
        debug_assert!(self.pnts[j].y > self.pnts[2].y, "pnts[j].p.y > pnts[2].p.y");

        let (_, c0, r0) = self.pnt_pnt_pnt_node(0, 1, 2);
        let (_, c1, r1) = self.pnt_pnt_pnt_node(1, 2, j);
        let (_, c2, r2) = self.pnt_pnt_pnt_node(2, i, j);
        let (_, c3, r3) = self.pnt_pnt_pnt_node(0, 2, i);
        let far = FAR_NODE_DIST2.sqrt();

        self.store_node(c0.x - far, c0.y, 0.0);
        self.store_node(c1.x, c1.y + far, 0.0);
        self.store_node(c2.x + far, c2.y, 0.0);
        self.store_node(c3.x, c3.y - far, 0.0);
        self.store_node(c0.x, c0.y, r0);
        self.store_node(c1.x, c1.y, r1);
        self.store_node(c2.x, c2.y, r2);
        self.store_node(c3.x, c3.y, r3);

        let pnt = SiteType::Pnt;
        self.store_edge(0, 4, 1, 0, pnt, pnt);
        self.store_edge(1, 5, j, 1, pnt, pnt);
        self.store_edge(2, 6, i, j, pnt, pnt);
        self.store_edge(3, 7, 0, i, pnt, pnt);
        self.store_edge(4, 5, 1, 2, pnt, pnt);
        self.store_edge(5, 6, j, 2, pnt, pnt);
        self.store_edge(6, 7, i, 2, pnt, pnt);
        self.store_edge(7, 4, 0, 2, pnt, pnt);
        self.store_edge(0, 1, NIL, 1, SiteType::NoSite, pnt);
        self.store_edge(1, 2, NIL, j, SiteType::NoSite, pnt);
        self.store_edge(2, 3, NIL, i, SiteType::NoSite, pnt);
        self.store_edge(3, 0, NIL, 0, SiteType::NoSite, pnt);

        const LINKS: [[usize; 4]; 12] = [
            [8, 11, 7, 4],
            [9, 8, 4, 5],
            [10, 9, 5, 6],
            [11, 10, 6, 7],
            [0, 7, 5, 1],
            [1, 4, 6, 2],
            [2, 5, 7, 3],
            [3, 6, 4, 0],
            [11, 0, 1, 9],
            [8, 1, 2, 10],
            [9, 2, 3, 11],
            [10, 3, 0, 8],
        ];
        for (edge, [s_ccw, s_cw, e_ccw, e_cw]) in self.edges.iter_mut().zip(LINKS) {
            edge.s_ccw = s_ccw;
            edge.s_cw = s_cw;
            edge.e_ccw = e_ccw;
            edge.e_cw = e_cw;
        }

        self.pnts[0].vd_edge = 11;
        self.pnts[1].vd_edge = 8;
        self.pnts[2].vd_edge = 4;
        self.pnts[i].vd_edge = 9;
        self.pnts[j].vd_edge = 10;
        for (k, node) in self.nodes.iter_mut().enumerate() {
            node.edge = k % 4;
        }

        self.insert_pnt_into_grid(2);
    }

    /// Splits edge `e = (n1, n2)` at node `n` into `e = (n1, n)` and a new edge `(n, n2)`.
    pub(crate) fn split_edge(&mut self, e: usize, n: usize, record_dummy: bool) {
        let Edge {
            lft,
            rgt,
            ltype,
            rtype,
            n2,
            ..
        } = self.edges[e];
        let e1 = self.store_edge(n, n2, lft, rgt, ltype, rtype);
        self.nodes[n].edge = e;
        if self.nodes[n2].edge == e {
            self.nodes[n2].edge = e1;
        }

        let e0 = self.ccw_edge(e, n2);
        self.set_ccw_edge(e1, n2, e0);
        self.set_cw_edge(e0, n2, e1);
        let e0 = self.cw_edge(e, n2);
        self.set_cw_edge(e1, n2, e0);
        self.set_ccw_edge(e0, n2, e1);

        self.edges[e].n2 = n;
        self.set_cw_edge(e, n, e1);
        self.set_ccw_edge(e, n, e1);
        self.set_cw_edge(e1, n, e);
        self.set_ccw_edge(e1, n, e);

        if record_dummy {
            self.dummy_nodes.push(DummyRec { e, e1, n });
            self.nodes[n].status = NodeStatus::Dummy;
        }
        self.nodes[n].deg2 = true;
    }

    pub(crate) fn remove_dummy_nodes(&mut self) {
        for k in 0..self.dummy_nodes.len() {
            let DummyRec { e, e1, n } = self.dummy_nodes[k];
            debug_assert!(self.nodes[n].status == NodeStatus::Dummy, "IsNodeDummy(n)");
            let (e1, e2) = if self.edges[e].n1 == n {
                (e1, e)
            } else {
                (e, e1)
            };
            debug_assert!(self.edges[e1].n2 == n, "GetEndNode(e1) == n");
            debug_assert!(self.edges[e2].n1 == n, "GetStartNode(e2) == n");
            debug_assert!(self.ccw_edge(e1, n) == e2, "GetCCWEdge(e1, n) == e2");
            debug_assert!(self.cw_edge(e1, n) == e2, "GetCWEdge(e1, n) == e2");
            debug_assert!(self.ccw_edge(e2, n) == e1, "GetCCWEdge(e2, n) == e1");
            debug_assert!(self.cw_edge(e2, n) == e1, "GetCWEdge(e2, n) == e1");

            let n2 = self.edges[e2].n2;
            self.edges[e1].n2 = n2;
            if self.nodes[n2].edge == e2 {
                self.nodes[n2].edge = e1;
            }
            let e0 = self.ccw_edge(e2, n2);
            self.set_ccw_edge(e1, n2, e0);
            self.set_cw_edge(e0, n2, e1);
            let e0 = self.cw_edge(e2, n2);
            self.set_cw_edge(e1, n2, e0);
            self.set_ccw_edge(e0, n2, e1);

            let Edge {
                lft,
                rgt,
                ltype,
                rtype,
                ..
            } = self.edges[e2];
            self.set_vd_ptr(lft, ltype, e1);
            self.set_vd_ptr(rgt, rtype, e1);
            self.del_node(n);
            self.del_edge(e2);
        }
        self.dummy_nodes.clear();
    }

    /// Deletes `e` and every dangling edge chain that its removal leaves behind.
    pub(crate) fn del_edge_rec(&mut self, mut e: usize) {
        loop {
            let n1 = self.edges[e].n1;
            let e1 = self.ccw_edge(e, n1);
            let e2 = self.cw_edge(e, n1);
            let mut dangling_start = NIL;
            if e1 == NIL {
                self.del_node(n1);
            } else {
                debug_assert!(e2 != NIL, "e2 != NIL");
                self.nodes[n1].edge = e1;
                if e1 == e2 {
                    self.set_cw_edge(e1, n1, NIL);
                    self.set_ccw_edge(e1, n1, NIL);
                    dangling_start = e1;
                } else {
                    self.set_cw_edge(e1, n1, e2);
                    self.set_ccw_edge(e2, n1, e1);
                    self.nodes[n1].deg2 = true;
                }
            }

            let n2 = self.edges[e].n2;
            let e1 = self.ccw_edge(e, n2);
            let e2 = self.cw_edge(e, n2);
            let mut dangling_end = NIL;
            if e1 == NIL {
                self.del_node(n2);
            } else {
                debug_assert!(e2 != NIL, "e2 != NIL");
                self.nodes[n2].edge = e1;
                if e1 != e2 {
                    self.set_cw_edge(e1, n2, e2);
                    self.set_ccw_edge(e2, n2, e1);
                    self.nodes[n2].deg2 = true;
                } else {
                    self.set_cw_edge(e1, n2, NIL);
                    self.set_ccw_edge(e1, n2, NIL);
                    dangling_end = e1;
                }
            }

            self.del_edge(e);
            if dangling_start != NIL && !self.edges[dangling_start].deleted {
                self.del_edge_rec(dangling_start);
            }
            if dangling_end == NIL || self.edges[dangling_end].deleted {
                return;
            }
            e = dangling_end;
        }
    }

    /// Removes the VD edges between consecutive approximation sites of the same arc.
    pub(crate) fn recover_arc_vd(&mut self) {
        if self.num_apx_arcs == 0 {
            return;
        }
        for pnt in &mut self.pnts {
            pnt.vd_edge = NIL;
        }
        for seg in &mut self.segs {
            seg.vd_edge = NIL;
        }

        let mut k = 0;
        while k < self.edges.len() {
            let edge = self.edges[k];
            if !edge.deleted && edge.lft != NIL && edge.rgt != NIL {
                let arc_l = self.apx_arc_of_site(edge.lft, edge.ltype);
                let arc_r = self.apx_arc_of_site(edge.rgt, edge.rtype);
                if arc_l != NIL && arc_l == arc_r {
                    let r1 = self.nodes[edge.n1].r;
                    let r2 = self.nodes[edge.n2].r;
                    if ((-ZERO..=ZERO).contains(&r1) && !(r2 <= ZERO))
                        || ((-ZERO..=ZERO).contains(&r2) && ZERO < r1)
                    {
                        self.del_edge_rec(k);
                    }
                }
            }
            k += 1;
        }

        let mut last = self.edges.len() - 1;
        while self.edges[last].deleted {
            last -= 1;
        }
        self.edges.truncate(last + 1);
        self.free_edges.retain(|&e| e <= last);

        for k in 0..self.edges.len() {
            let Edge {
                lft,
                rgt,
                ltype,
                rtype,
                deleted,
                ..
            } = self.edges[k];
            if deleted {
                continue;
            }
            if lft != NIL {
                self.set_vd_ptr_with_apx_arc(lft, ltype, k);
            }
            if rgt != NIL {
                self.set_vd_ptr_with_apx_arc(rgt, rtype, k);
            }
        }
    }

    pub(crate) fn apx_arc_of_site(&self, site: usize, site_type: SiteType) -> usize {
        match site_type {
            SiteType::Pnt => self.pnts[site].apx_arc,
            SiteType::Seg => self.segs[site].apx_arc,
            SiteType::Arc | SiteType::NoSite => NIL,
        }
    }

    fn set_vd_ptr_with_apx_arc(&mut self, site: usize, site_type: SiteType, e: usize) {
        self.set_vd_ptr(site, site_type, e);
        let arc = self.apx_arc_of_site(site, site_type);
        if arc != NIL {
            self.arcs[arc].vd_edge = e;
        }
    }
}
