use crate::consts::*;
use crate::real::Real;
use crate::state::*;

const FIRST_INNER_PNT: usize = 2;

#[derive(Clone, Copy, PartialEq, Eq)]
enum BisectorKind {
    Straight,
    Parabola,
    PntPnt,
}

impl<R: Real> Engine<R> {
    pub(crate) fn offset_isolated_pnts(&mut self, t: f64, dt: f64) {
        let end = self.pnts.len().saturating_sub(2);
        for i in FIRST_INNER_PNT..end {
            let p = self.pnts[i];
            if p.has_deg2 || p.deleted || p.vd_edge == NIL {
                continue;
            }
            let e0 = p.vd_edge;
            let edge = self.edges[e0];
            let mut n = if edge.lft == i && edge.ltype == SiteType::Pnt {
                edge.n2
            } else {
                debug_assert!(
                    edge.rgt == i && edge.rtype == SiteType::Pnt,
                    "IsRgtSite(e, i, PNT)"
                );
                edge.n1
            };
            let mut r_min = self.nodes[n].r;
            let mut e = e0;
            loop {
                n = self.other_node(e, n);
                let r = self.nodes[n].r;
                if r <= r_min {
                    r_min = r;
                }
                e = self.ccw_edge(e, n);
                if e == e0 {
                    break;
                }
            }
            let mut tt = t;
            while tt <= r_min {
                let start = self.offset_elems.len();
                self.offset_loops.push(OffsetLoop {
                    start,
                    end: start + 1,
                    dist: (R::of(tt) / self.scale_factor).stored(),
                });
                self.offset_elems.push(OffsetElem {
                    site: i,
                    kind: OffsetKind::PntArc,
                    x: (R::of(tt) + p.x).stored(),
                    y: p.y,
                });
                self.offset_elems.push(OffsetElem {
                    site: i,
                    kind: OffsetKind::PntArc,
                    x: (R::of(p.x) - tt).stored(),
                    y: p.y,
                });
                tt = (R::of(tt) + dt).stored();
            }
        }
    }

    fn bisector_kind(&self, e: usize) -> BisectorKind {
        let edge = &self.edges[e];
        let pnt_is_seg_endpoint = |pnt: usize, seg: usize| {
            let s = &self.segs[seg];
            s.i1 == pnt || s.i2 == pnt
        };
        if edge.ltype == SiteType::Pnt {
            if edge.rtype == SiteType::Pnt {
                return BisectorKind::PntPnt;
            }
            debug_assert!(edge.rtype == SiteType::Seg, "t_rgt == SEG");
            if !pnt_is_seg_endpoint(edge.lft, edge.rgt) {
                return BisectorKind::Parabola;
            }
        } else if edge.rtype == SiteType::Pnt {
            debug_assert!(edge.ltype == SiteType::Seg, "t_lft == SEG");
            if !pnt_is_seg_endpoint(edge.rgt, edge.lft) {
                return BisectorKind::Parabola;
            }
        }
        BisectorKind::Straight
    }

    fn trace_offset_loop(&mut self, e0: usize, t: f64) {
        let loop_idx = self.offset_loops.len();
        let start = self.offset_elems.len();
        self.offset_loops.push(OffsetLoop {
            start,
            end: start,
            dist: (R::of(t) / self.scale_factor).stored(),
        });

        let mut e = e0;
        loop {
            let edge = self.edges[e];
            let kind = self.bisector_kind(e);
            if self.e_data[e].pending {
                match kind {
                    BisectorKind::Straight => self.create_ss_bisector_data(e),
                    BisectorKind::Parabola => self.create_parabola_data(e),
                    BisectorKind::PntPnt => self.create_pp_bisector_data(e),
                }
            }
            self.edge_flags[e].available = false;
            let (x, y) = match kind {
                BisectorKind::Straight => self.eval_ss_bisector_pnt(e, t),
                BisectorKind::Parabola => {
                    let seg = if edge.ltype == SiteType::Seg {
                        edge.lft
                    } else {
                        edge.rgt
                    };
                    self.eval_parabola_pnt(e, seg, t)
                }
                BisectorKind::PntPnt => self.eval_pp_bisector_pnt(e, t),
            };

            let r1 = self.nodes[edge.n1].r;
            let r2 = self.nodes[edge.n2].r;
            let (site, stype, mut n) = if r2 <= t {
                debug_assert!(t < r1, "t < t1");
                (edge.lft, edge.ltype, edge.n1)
            } else {
                (edge.rgt, edge.rtype, edge.n2)
            };
            let apx = self.apx_arc_of_site(site, stype);
            let elem = if apx == NIL {
                let kind = match stype {
                    SiteType::Seg => OffsetKind::Line,
                    SiteType::Pnt => OffsetKind::PntArc,
                    SiteType::Arc | SiteType::NoSite => {
                        unreachable!("offset element from unsupported site type")
                    }
                };
                OffsetElem { site, kind, x, y }
            } else {
                let arc = &self.arcs[apx];
                let dx = (R::of(arc.cx) - x).stored();
                let dy = (R::of(arc.cy) - y).stored();
                let kind = if R::of(arc.r) * arc.r <= R::of(dx) * dx + R::of(dy) * dy {
                    OffsetKind::ApxArcOutside
                } else {
                    OffsetKind::ApxArcInside
                };
                OffsetElem {
                    site: apx,
                    kind,
                    x,
                    y,
                }
            };
            self.offset_elems.push(elem);

            loop {
                let next = self.ccw_edge(e, n);
                let next_edge = &self.edges[next];
                let nr1 = self.nodes[next_edge.n1].r;
                let nr2 = self.nodes[next_edge.n2].r;
                let r_min = if nr2 <= nr1 { nr2 } else { nr1 };
                let r_max = if nr1 <= nr2 { nr2 } else { nr1 };
                e = next;
                if r_min <= t && t < r_max {
                    break;
                }
                n = self.other_node(e, n);
            }
            if e == e0 {
                self.offset_loops[loop_idx].end = self.offset_elems.len() - 1;
                return;
            }
        }
    }

    pub(crate) fn compute_offset(&mut self, t: f64, dt: f64, left: bool, right: bool) -> bool {
        debug_assert!(self.vd_computed, "GetVDStatus()");
        if self.nodes_store_r2 {
            self.nodes_r2_to_r();
        }
        if !self.deg2_nodes_added {
            self.insert_pnt_pnt_mid_nodes();
        }
        if self.vd_id != self.offset_vd_id {
            self.new_offset_needed = true;
            self.offset_vd_id = self.vd_id;
        }

        let mut finished;
        if self.new_offset_needed {
            self.offset_t = (R::of(t) * self.scale_factor).stored();
            self.offset_dt = (R::of(self.scale_factor) * dt).stored();
            if self.max_offset <= self.offset_t || self.offset_t <= ZERO {
                return true;
            }
            finished = false;
            let both_sides = !left && !right;
            self.initialize_edge_data(both_sides, left);
            if self.offset_dt <= ZERO {
                self.offset_dt = self.max_offset;
            }
            if self.isolated_pnts && both_sides {
                self.offset_isolated_pnts(self.offset_t, self.offset_dt);
            }
        } else {
            self.offset_t = (R::of(self.offset_t) + self.offset_dt).stored();
            finished = self.max_offset <= self.offset_t;
        }

        while !finished {
            let mut j = 0;
            while j < self.active_edges.len() {
                let e = self.active_edges[j];
                if !self.edge_flags[e].available {
                    self.edge_flags[e].available = true;
                    j += 1;
                    continue;
                }
                let edge = &self.edges[e];
                let r1 = self.nodes[edge.n1].r;
                let r2 = self.nodes[edge.n2].r;
                let r_min = if r1 < r2 { r1 } else { r2 };
                if self.offset_t < r_min {
                    j += 1;
                    continue;
                }
                let r_max = if r2 < r1 { r1 } else { r2 };
                if self.offset_t < r_max {
                    self.trace_offset_loop(e, self.offset_t);
                    self.edge_flags[e].available = true;
                    j += 1;
                } else {
                    self.active_edges.swap_remove(j);
                }
            }
            self.offset_t = (R::of(self.offset_t) + self.offset_dt).stored();
            finished = self.max_offset <= self.offset_t;
        }
        finished
    }
}
