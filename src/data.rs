use crate::consts::*;
use crate::real::Real;
use crate::state::*;
use std::marker::PhantomData;

/// Index at which the first stored input point lands; slots 0 and 1 are the left bbox corners.
pub(crate) const FIRST_INPUT_PNT: usize = 2;
pub(crate) const NUM_BBOX_RGT_PNTS: usize = 2;
const MIN_SCALED_EXTENT: f64 = 0.5;

impl<R: Real> Engine<R> {
    pub(crate) fn new() -> Self {
        Self {
            pnts: Vec::new(),
            segs: Vec::new(),
            arcs: Vec::new(),
            num_arcs: 0,
            num_apx_arcs: 0,
            scale_factor: 1.0,
            shift_x: 0.0,
            shift_y: 0.0,
            data_scaled: false,
            first_pnt_stored: false,
            bbox_added: false,
            bbox_rgt_added: false,
            data_cleaned: false,
            bb_min_x: 0.0,
            bb_min_y: 0.0,
            bb_max_x: 0.0,
            bb_max_y: 0.0,
            min_bbox_pnt: NIL,
            max_bbox_pnt: NIL,
            max_offset: 0.0,
            large: (R::of(f64::MAX).sqrt() / LARGE_DIVISOR).stored(),
            isolated_pnts: false,
            failure: None,

            reverse_ind: Vec::new(),
            reverse_index_exists: false,

            restart: false,
            restart_count: 0,
            troubles: false,
            desperate_used: false,

            nodes: Vec::new(),
            edges: Vec::new(),
            free_nodes: Vec::new(),
            free_edges: Vec::new(),
            dummy_nodes: Vec::new(),
            nodes_store_r2: true,
            vd_computed: false,
            deg2_nodes_added: false,
            vd_id: 0,
            between_eps: 0.0,
            pnt_visited: Vec::new(),
            nodes_checked: Vec::new(),
            check_pnts: Vec::new(),
            loop_stack: Vec::new(),
            loop_stack_start: 0,

            grid_active: false,
            grid: Vec::new(),
            grid_x_coords: Vec::new(),
            grid_y_coords: Vec::new(),
            grid_pnt_list: Vec::new(),
            grid_first_n_pnts: Vec::new(),
            grid_first_counter: 0,
            grid_n_initial: 0,
            grid_nx: 0,
            grid_ny: 0,
            grid_n: 0,
            grid_min_x: 0.0,
            grid_min_y: 0.0,
            grid_max_x: 0.0,
            grid_max_y: 0.0,
            grid_delta_x: 0.0,
            grid_delta_y: 0.0,
            grid_cell_w: 0.0,
            grid_cell_h: 0.0,
            grid_aspect: 0.0,
            grid_inv_aspect: 0.0,

            tree: Vec::new(),
            tree_root: NIL,
            tree_xmin: 0.0,
            tree_ymin: 0.0,
            tree_xmax: 0.0,
            tree_ymax: 0.0,

            heap: Vec::new(),
            heap_top_pending: false,
            heap_top_valid: false,

            e_data: Vec::new(),
            e_data_vd_id: 0,
            edge_flags: Vec::new(),
            active_edges: Vec::new(),
            offset_loops: Vec::new(),
            offset_elems: Vec::new(),
            offset_t: 0.0,
            offset_dt: 0.0,
            offset_vd_id: 0,
            new_offset_needed: true,
            arithmetic: PhantomData,
        }
    }

    pub(crate) fn store_pnt(&mut self, x: f64, y: f64) -> usize {
        let pnt = Pnt {
            x,
            y,
            vd_edge: NIL,
            deg2_node: NIL,
            visited: false,
            has_deg2: false,
            deleted: false,
            apx_arc: NIL,
        };
        if !self.first_pnt_stored {
            self.pnts.clear();
            self.pnts.extend([pnt, pnt]);
            self.first_pnt_stored = true;
        }
        self.pnts.push(pnt);
        self.pnts.len() - 1
    }

    pub(crate) fn store_seg(&mut self, i1: usize, i2: usize) -> usize {
        self.segs.push(Seg {
            i1,
            i2,
            vd_flag: false,
            a: 0.0,
            b: 0.0,
            c: 0.0,
            lgth: 0.0,
            vd_edge: NIL,
            apx_arc: NIL,
        });
        self.segs.len() - 1
    }

    pub(crate) fn store_arc(&mut self, i1: usize, i2: usize, cx: f64, cy: f64, ccw: bool) -> usize {
        let p1 = self.pnts[i1];
        let p2 = self.pnts[i2];
        let r1 = dist::<R>((R::of(p1.x) - cx).stored(), (R::of(p1.y) - cy).stored()).stored();
        let r2 = dist::<R>((R::of(p2.x) - cx).stored(), (R::of(p2.y) - cy).stored()).stored();
        let arc = Arc {
            i1,
            i2,
            cx,
            cy,
            r: ((R::of(r2) + r1) * 0.5).stored(),
            ccw,
            vd_edge: NIL,
        };
        if self.num_arcs < self.arcs.len() {
            self.arcs[self.num_arcs] = arc;
        } else {
            self.arcs.push(arc);
        }
        self.num_arcs += 1;
        self.num_arcs - 1
    }

    pub(crate) fn merge_pnts(&mut self, i1: usize, i2: usize) {
        let x = ((R::of(self.pnts[i2].x) + self.pnts[i1].x) * 0.5).stored();
        let y = ((R::of(self.pnts[i2].y) + self.pnts[i1].y) * 0.5).stored();
        self.pnts[i1].x = x;
        self.pnts[i1].y = y;
        self.pnts[i2].x = x;
        self.pnts[i2].y = y;
        self.restart = true;
    }

    pub(crate) fn delete_seg(&mut self, i: usize) {
        self.segs.swap_remove(i);
        self.restart = true;
    }

    pub(crate) fn delete_arc(&mut self, i: usize) {
        let last = self.num_arcs - 1;
        self.arcs[i] = self.arcs[last];
        self.num_arcs = last;
        if self.num_arcs < self.num_apx_arcs {
            self.num_apx_arcs = self.num_arcs;
        }
        self.restart = true;
    }

    pub(crate) fn reset_site_vd_data(&mut self) {
        for p in &mut self.pnts {
            p.vd_edge = NIL;
            p.deg2_node = NIL;
            p.visited = false;
            p.has_deg2 = false;
        }
        for s in &mut self.segs {
            s.vd_flag = false;
            s.vd_edge = NIL;
        }
        let n = if self.num_arcs == 0 {
            self.num_apx_arcs
        } else {
            self.num_arcs
        };
        for a in &mut self.arcs[..n] {
            a.vd_edge = NIL;
        }
    }

    pub(crate) fn compute_bbox(&mut self) {
        if self.pnts.is_empty() {
            return;
        }
        debug_assert!(!self.bbox_added, "!bbox_added");
        debug_assert!(!self.bbox_rgt_added, "!bbox_rgt_added");
        self.bb_min_x = self.pnts[0].x;
        self.bb_max_x = self.bb_min_x;
        self.bb_min_y = self.pnts[0].y;
        self.bb_max_y = self.bb_min_y;
        for p in &self.pnts[1..] {
            if self.bb_min_x <= p.x {
                if self.bb_max_x < p.x {
                    self.bb_max_x = p.x;
                }
            } else {
                self.bb_min_x = p.x;
            }
            if self.bb_min_y <= p.y {
                if self.bb_max_y < p.y {
                    self.bb_max_y = p.y;
                }
            } else {
                self.bb_min_y = p.y;
            }
        }
    }

    pub(crate) fn scale_data(&mut self) {
        let dx = (R::of(self.bb_max_x) - self.bb_min_x).stored();
        let dy = (R::of(self.bb_max_y) - self.bb_min_y).stored();
        let mut d = if dy < dx { dx } else { dy };
        if !(ZERO < d) {
            d = ZERO;
        }
        let mut s: f64 = 1.0;
        self.scale_factor = 1.0;
        if ONE_PLUS_ZERO < d {
            loop {
                s = (R::of(s) * 0.5).stored();
                if !(R::of(s) * d > ONE_PLUS_ZERO) {
                    break;
                }
            }
            self.data_scaled = true;
        } else if d < MIN_SCALED_EXTENT {
            loop {
                s = (R::of(s) + s).stored();
                if !(R::of(s) * d < MIN_SCALED_EXTENT) {
                    break;
                }
            }
            self.data_scaled = true;
        } else if self.bb_max_x > ONE_PLUS_ZERO
            || self.bb_max_y > ONE_PLUS_ZERO
            || self.bb_min_x < -ZERO
            || self.bb_min_y < -ZERO
        {
            self.data_scaled = true;
        }

        if !self.data_scaled {
            self.shift_x = 0.0;
            self.shift_y = 0.0;
            self.scale_factor = 1.0;
            return;
        }

        self.scale_factor = s;
        let (min_x, min_y) = (self.bb_min_x, self.bb_min_y);
        let n_arcs = self.num_arcs;
        for p in &mut self.pnts {
            p.x = ((R::of(p.x) - min_x) * s).stored();
            p.y = ((R::of(p.y) - min_y) * s).stored();
        }
        for a in &mut self.arcs[..n_arcs] {
            a.cx = ((R::of(a.cx) - min_x) * s).stored();
            a.cy = ((R::of(a.cy) - min_y) * s).stored();
            a.r = (R::of(a.r) * s).stored();
        }
        self.shift_x = min_x;
        self.shift_y = min_y;
        self.bb_max_x = ((R::of(self.bb_max_x) - min_x) * s).stored();
        self.bb_max_y = ((R::of(self.bb_max_y) - min_y) * s).stored();
        self.bb_min_x = 0.0;
        self.bb_min_y = 0.0;
    }

    pub(crate) fn add_bbox(&mut self) {
        debug_assert!(!self.bbox_added, "!bbox_added");
        debug_assert!(self.first_pnt_stored, "bbox_lft_added");
        debug_assert!(self.bbox_rgt_added, "bbox_rgt_added");
        let factor = (R::of(f64::from(BBOX_BOUND)) + BBOX_MARGIN_FRACTION).stored();
        let dx = (R::of(self.bb_max_x) - self.bb_min_x).stored();
        let mut mx = if ZERO < dx || dx < -ZERO {
            (R::of(factor) * dx).stored()
        } else {
            factor
        };
        let dy = (R::of(self.bb_max_y) - self.bb_min_y).stored();
        let mut my = if dy > ZERO || dy < -ZERO || dy.is_nan() {
            (R::of(dy) * factor).stored()
        } else {
            factor
        };
        if R::of(my) * 0.5 <= mx {
            if R::of(mx) * 0.5 > my {
                my = (R::of(mx) * 0.5).stored();
            }
        } else {
            mx = (R::of(my) * 0.5).stored();
        }

        let n = self.pnts.len();
        debug_assert!(self.pnts[0].x == -f64::MAX && self.pnts[0].y == -f64::MAX);
        self.pnts[0].x = (R::of(self.bb_min_x) - mx).stored();
        self.pnts[0].y = (R::of(self.bb_min_y) - my).stored();
        debug_assert!(self.pnts[1].x == -f64::MAX && self.pnts[1].y == f64::MAX);
        self.pnts[1].x = (R::of(self.bb_min_x) - mx).stored();
        self.pnts[1].y = (R::of(self.bb_max_y) + my).stored();
        debug_assert!(self.pnts[n - 2].x == f64::MAX && self.pnts[n - 2].y == -f64::MAX);
        self.pnts[n - 2].x = (R::of(self.bb_max_x) + mx).stored();
        self.pnts[n - 2].y = (R::of(self.bb_min_y) - my).stored();
        debug_assert!(self.pnts[n - 1].x == f64::MAX && self.pnts[n - 1].y == f64::MAX);
        self.pnts[n - 1].x = (R::of(self.bb_max_x) + mx).stored();
        self.pnts[n - 1].y = (R::of(self.bb_max_y) + my).stored();

        self.max_offset = ((R::of(mx) * mx + R::of(my) * my).sqrt() * MAX_OFFSET_FRACTION).stored();
        self.bbox_added = true;
    }

    pub(crate) fn find_isolated_pnts(&mut self) {
        debug_assert!(self.bbox_added, "bbox_added");
        let mut i = FIRST_INPUT_PNT;
        while self.pnts[i].has_deg2 || self.pnts[i].deleted {
            i += 1;
        }
        self.isolated_pnts = i < self.pnts.len() - NUM_BBOX_RGT_PNTS;
    }

    pub(crate) fn discard_apx_arc_refs(&mut self) {
        if self.num_apx_arcs > 0 {
            for p in &mut self.pnts {
                p.apx_arc = NIL;
            }
            for s in &mut self.segs {
                s.apx_arc = NIL;
            }
            self.num_apx_arcs = 0;
        }
    }

    pub(crate) fn split_seg(&mut self, i_seg: usize, i_pnt: usize) {
        let i2 = self.segs[i_seg].i2;
        let new_seg = self.store_seg(i_pnt, i2);
        self.segs[i_seg].i2 = i_pnt;
        self.pnts[i_pnt].apx_arc = NIL;
        self.segs[i_seg].apx_arc = NIL;
        self.segs[new_seg].apx_arc = NIL;
        self.restart = true;
    }
}

pub(crate) fn dist<R: Real>(dx: f64, dy: f64) -> R {
    (R::of(dy) * dy + R::of(dx) * dx).sqrt()
}
