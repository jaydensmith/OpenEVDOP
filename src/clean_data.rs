use crate::data::{FIRST_INPUT_PNT, NUM_BBOX_RGT_PNTS};
use crate::misc::quicksort;
use crate::state::*;

fn pnt_compare(a: &Pnt, b: &Pnt) -> i32 {
    if b.x > a.x {
        -1
    } else if b.x < a.x {
        1
    } else if b.y > a.y {
        -1
    } else if b.y < a.y {
        1
    } else {
        0
    }
}

impl Engine {
    fn num_arc_records(&self) -> usize {
        if self.num_arcs == 0 {
            self.num_apx_arcs
        } else {
            self.num_arcs
        }
    }

    fn build_reverse_index(&mut self, n_segs: usize, n_arcs: usize, store_orig_index: bool) {
        let n = self.pnts.len();
        self.reverse_ind.clear();
        self.reverse_ind.resize(
            n,
            ReverseInd {
                isite: NIL,
                site_type: ReverseType::IsoPnt,
                next: NIL,
            },
        );
        if store_orig_index {
            for (i, p) in self.pnts.iter_mut().enumerate() {
                p.vd_edge = i;
            }
        }
        for j in 0..n_segs {
            let (i1, i2) = (self.segs[j].i1, self.segs[j].i2);
            self.add_reverse_ref(i1, j, ReverseType::SegStart);
            self.add_reverse_ref(i2, j, ReverseType::SegEnd);
        }
        for j in 0..n_arcs {
            let (i1, i2) = (self.arcs[j].i1, self.arcs[j].i2);
            self.add_reverse_ref(i1, j, ReverseType::ArcStart);
            self.add_reverse_ref(i2, j, ReverseType::ArcEnd);
        }
        self.reverse_index_exists = true;
    }

    fn add_reverse_ref(&mut self, k: usize, isite: usize, site_type: ReverseType) {
        if self.reverse_ind[k].site_type == ReverseType::IsoPnt {
            debug_assert!(
                self.reverse_ind[k].next == NIL,
                "reverse_ind[k].next == NIL"
            );
            self.reverse_ind[k].isite = isite;
            self.reverse_ind[k].site_type = site_type;
        } else {
            let entry = ReverseInd {
                isite,
                site_type,
                next: self.reverse_ind[k].next,
            };
            self.reverse_ind.push(entry);
            self.reverse_ind[k].next = self.reverse_ind.len() - 1;
        }
    }

    pub(crate) fn free_reverse_index(&mut self) {
        self.reverse_ind.clear();
        self.reverse_index_exists = false;
    }

    fn set_site_endpoint(&mut self, r: ReverseInd, pnt: usize) {
        match r.site_type {
            ReverseType::SegStart => self.segs[r.isite].i1 = pnt,
            ReverseType::SegEnd => self.segs[r.isite].i2 = pnt,
            ReverseType::ArcStart => self.arcs[r.isite].i1 = pnt,
            ReverseType::ArcEnd => self.arcs[r.isite].i2 = pnt,
            ReverseType::IsoPnt => {}
        }
    }

    /// Makes every segment/arc endpoint that refers to `i_del` refer to `i_old` instead.
    pub(crate) fn redirect_pnt(&mut self, i_del: usize, i_old: usize) {
        if i_del == i_old {
            return;
        }
        if self.data_cleaned {
            self.restart = true;
            return;
        }
        if !self.reverse_index_exists {
            let n_arcs = self.num_arc_records();
            self.build_reverse_index(self.segs.len(), n_arcs, false);
        }
        let mut k = i_del;
        while k != NIL {
            let r = self.reverse_ind[k];
            self.set_site_endpoint(r, i_old);
            k = r.next;
        }
    }

    /// First call: turns the placeholder points into the four bbox corners. Later calls (after a
    /// restart) sort the points, drop duplicates, rewire all endpoints and drop degenerate segments.
    pub(crate) fn clean_data(&mut self) {
        let n_arcs = self.num_arc_records();
        if !self.bbox_rgt_added {
            debug_assert!(self.first_pnt_stored, "bbox_lft_added");
            debug_assert!(!self.bbox_added, "!bbox_added");
            debug_assert!(self.pnts[2].x == self.pnts[0].x && self.pnts[2].x == self.pnts[1].x);
            debug_assert!(self.pnts[2].y == self.pnts[0].y && self.pnts[2].y == self.pnts[1].y);
            self.pnts[0].x = -f64::MAX;
            self.pnts[0].y = -f64::MAX;
            self.pnts[1].x = -f64::MAX;
            self.pnts[1].y = f64::MAX;
            self.store_pnt(f64::MAX, -f64::MAX);
            self.store_pnt(f64::MAX, f64::MAX);
            self.bbox_rgt_added = true;
            self.data_cleaned = false;
            return;
        }

        self.isolated_pnts = false;
        if !self.segs.is_empty() || n_arcs > 0 {
            self.build_reverse_index(self.segs.len(), n_arcs, true);
        }
        quicksort(&mut self.pnts, pnt_compare);
        self.data_cleaned = true;

        let n = self.pnts.len();
        let mut kept = 0;
        if self.segs.is_empty() && n_arcs == 0 {
            for i in 1..n {
                let p = self.pnts[i];
                if !p.deleted && (self.pnts[kept].x != p.x || self.pnts[kept].y != p.y) {
                    kept += 1;
                    self.pnts[kept] = p;
                }
            }
            self.pnts.truncate(kept + 1);
            self.isolated_pnts = true;
            return;
        }

        let last_interior = n - NUM_BBOX_RGT_PNTS;
        for i in 0..n {
            let head = self.pnts[i].vd_edge;
            self.pnts[i].vd_edge = NIL;
            let p = self.pnts[i];
            if p.deleted {
                continue;
            }
            if self.pnts[kept].x != p.x || self.pnts[kept].y != p.y {
                kept += 1;
                self.pnts[kept] = p;
            }
            let mut k = head;
            while k != NIL {
                let r = self.reverse_ind[k];
                if r.site_type == ReverseType::IsoPnt {
                    if i >= FIRST_INPUT_PNT && i < last_interior {
                        self.isolated_pnts = true;
                    }
                } else {
                    self.set_site_endpoint(r, kept);
                }
                k = r.next;
            }
        }
        self.pnts.truncate(kept + 1);

        let mut j = 0;
        while j < self.segs.len() {
            if self.segs[j].i1 == self.segs[j].i2 {
                self.segs.swap_remove(j);
            } else {
                j += 1;
            }
        }
        self.free_reverse_index();
    }
}
