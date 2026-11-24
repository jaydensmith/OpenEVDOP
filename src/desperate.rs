use crate::consts::*;
use crate::error::OffsetError;
use crate::state::*;

const NO_DIST_YET: f64 = -1.0;

impl Engine {
    /// `None` after aborting the computation.
    pub(crate) fn desperate_node(
        &mut self,
        e: usize,
        i: usize,
        t: SiteType,
    ) -> Option<(Point, f64)> {
        self.troubles = true;
        let edge = self.edges[e];
        let (lft, rgt) = (edge.lft, edge.rgt);
        let pnt = |v: &Engine, k: usize| Point {
            x: v.pnts[k].x,
            y: v.pnts[k].y,
        };
        let (center, r) = match (edge.ltype, edge.rtype, t) {
            (SiteType::Pnt, SiteType::Pnt, SiteType::Seg) => {
                self.seg_pnt_pnt_desperate(i, pnt(self, lft), pnt(self, rgt))
            }
            (SiteType::Pnt, SiteType::Seg, SiteType::Seg) => {
                self.seg_seg_pnt_desperate(rgt, i, pnt(self, lft))
            }
            (SiteType::Seg, SiteType::Pnt, SiteType::Seg) => {
                self.seg_seg_pnt_desperate(lft, i, pnt(self, rgt))
            }
            (SiteType::Seg, SiteType::Seg, SiteType::Seg) => {
                self.seg_seg_seg_desperate(lft, rgt, i)
            }
            _ => {
                self.abort_with(OffsetError::NumericalFailure);
                return None;
            }
        };
        self.desperate_used = true;
        debug_assert!(r >= -ZERO, "ge(*r2, ZERO)");
        Some((center, r))
    }

    fn clearance_beyond_seg(&self, node: usize, i: usize, s: usize) -> f64 {
        let n = &self.nodes[node];
        if n.site == i {
            return 0.0;
        }
        self.seg_clearance(s, n.x, n.y, SMALL) - n.r
    }

    pub(crate) fn desperate_split_pnt_edge(&mut self, i: usize, s: usize) -> usize {
        let e = self.pnts[i].vd_edge;
        let n_end = self.edges[e].n2;
        let n_start = self.edges[e].n1;
        let d_end = self.clearance_beyond_seg(n_end, i, s);
        let d_start = self.clearance_beyond_seg(n_start, i, s);
        let chosen = self.nodes[if d_end < d_start { n_end } else { n_start }];
        let n = self.store_node(chosen.x, chosen.y, chosen.r);
        self.split_edge(e, n, false);
        self.desperate_used = true;
        n
    }

    pub(crate) fn desperate_split_cell_max_dist(&mut self, i: usize, t: SiteType, s: usize) {
        self.troubles = true;
        let e0 = self.get_vd_ptr(i, t);
        let edge = self.edges[e0];
        debug_assert!(
            (edge.lft == i && edge.ltype == t) || (edge.rgt == i && edge.rtype == t),
            "IsLftSite(e, i, t) || IsRgtSite(e, i, t)"
        );
        let seg = self.segs[s];
        let (mut k, start) = if edge.rgt == i && edge.rtype == t {
            (edge.n2, edge.n1)
        } else {
            (edge.n1, edge.n2)
        };
        let mut best = self.nodes[start];
        let mut best_dist = seg.signed_dist(best.x, best.y).abs();
        let mut best_edge = e0;
        let mut e = e0;
        loop {
            let node = self.nodes[k];
            let dist = seg.signed_dist(node.x, node.y).abs();
            if best_dist < dist {
                best_edge = e;
                best = node;
                best_dist = dist;
            }
            e = self.ccw_edge(e, k);
            k = self.other_node(e, k);
            if e == e0 {
                break;
            }
        }
        let n = self.store_node(best.x, best.y, best.r);
        self.split_edge(best_edge, n, true);
        self.desperate_used = true;
    }

    /// False after aborting the computation because no loop node qualifies.
    pub(crate) fn desperate_split_loop_max_dist(&mut self, s: usize) -> bool {
        self.troubles = true;
        let seg = self.segs[s];
        let start = self.pnts[seg.i1];
        let mut best = (start.x, start.y, 0.0);
        let mut best_dist = NO_DIST_YET;
        let mut best_edge = NIL;
        for j in self.loop_stack_start..self.loop_stack.len() {
            let e = self.loop_stack[j];
            let edge = &self.edges[e];
            for node in [self.nodes[edge.n2], self.nodes[edge.n1]] {
                let dist = seg.signed_dist(node.x, node.y).abs();
                if best_dist < dist {
                    best = (node.x, node.y, node.r);
                    best_dist = dist;
                    best_edge = e;
                }
            }
        }
        if best_edge == NIL {
            self.abort_with(OffsetError::NumericalFailure);
            return false;
        }
        let n = self.store_node(best.0, best.1, best.2);
        self.split_edge(best_edge, n, true);
        self.desperate_used = true;
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_vd() -> Engine {
        let mut v = Engine::new();
        let sq = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        for i in 0..4 {
            let (a, b) = (sq[i], sq[(i + 1) % 4]);
            v.api_add_seg(a.0, a.1, b.0, b.1);
        }
        v.api_compute_vd(true, 0.001).expect("square VD");
        v
    }

    #[test]
    fn empty_conflict_loop_aborts_instead_of_splitting_a_missing_edge() {
        let mut v = square_vd();
        v.loop_stack.clear();
        v.loop_stack_start = 0;
        let edges_before = v.edges.len();
        assert!(!v.desperate_split_loop_max_dist(0));
        assert_eq!(v.failure, Some(OffsetError::NumericalFailure));
        assert!(v.restart, "callers unwind through the restart checks");
        assert_eq!(v.edges.len(), edges_before, "no edge was split");
    }

    #[test]
    fn unsupported_site_combination_aborts_instead_of_exiting() {
        let mut v = square_vd();
        let e = v
            .edges
            .iter()
            .position(|e| !e.deleted && e.ltype == SiteType::Pnt && e.rtype == SiteType::Pnt)
            .expect("a point-point edge");
        assert_eq!(v.desperate_node(e, 2, SiteType::Pnt), None);
        assert_eq!(v.failure, Some(OffsetError::NumericalFailure));
        assert!(v.restart);
    }

    #[test]
    fn driver_reports_an_abort_as_an_error_instead_of_restarting() {
        let mut v = square_vd();
        v.abort_with(OffsetError::NumericalFailure);
        assert_eq!(
            v.api_compute_vd(true, 0.001),
            Err(OffsetError::NumericalFailure)
        );
        assert_eq!(
            v.restart_count, 1,
            "exactly one pass ran before the error surfaced"
        );
    }
}
