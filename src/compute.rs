use crate::consts::*;
use crate::data::NUM_BBOX_RGT_PNTS;
use crate::error::OffsetError;
use crate::state::*;

/// Points 0..=2 (two left bbox corners and the first input point) seed the initial VD.
const FIRST_INSERTED_PNT: usize = 3;

impl Engine {
    pub(crate) fn reset_vd(&mut self) {
        self.nodes.clear();
        self.edges.clear();
        self.free_nodes.clear();
        self.free_edges.clear();
        self.dummy_nodes.clear();
        self.vd_computed = false;
        self.deg2_nodes_added = false;
        self.nodes_store_r2 = true;
    }

    pub(crate) fn reset_vd_data(&mut self) {
        self.pnt_visited.clear();
        self.check_pnts.clear();
        self.loop_stack.clear();
        self.nodes_checked.clear();
    }

    pub(crate) fn free_insert_buffers(&mut self) {
        self.heap.clear();
        self.tree.clear();
        self.tree_root = NIL;
    }

    /// Abandons the computation: callers unwind through the existing restart checks and the
    /// driver reports `err` instead of restarting.
    pub(crate) fn abort_with(&mut self, err: OffsetError) {
        self.failure = Some(err);
        self.restart = true;
    }

    /// One pass of the Voronoi diagram computation; returns early with `self.restart` set when the input
    /// had to be repaired, in which case the caller runs it again.
    pub(crate) fn compute_vd(
        &mut self,
        recover_arcs: bool,
        apx_abs: f64,
    ) -> Result<(), OffsetError> {
        if self.pnts.is_empty() {
            return Ok(());
        }
        if self.restart {
            self.restart_count += 1;
            if self.restart_count > MAX_RESTARTS {
                return Err(OffsetError::TooManyRestarts);
            }
            self.reset_vd();
            self.reset_vd_data();
            self.reset_site_vd_data();
        } else {
            self.restart_count = 0;
            self.data_cleaned = false;
        }
        self.restart = false;

        if !self.bbox_added {
            self.compute_bbox();
            self.scale_data();
        }
        if self.num_arcs > 0 && apx_abs > 0.0 {
            self.approx_arcs(apx_abs);
        }
        self.clean_data();
        if self.num_arcs > 0 {
            return Err(OffsetError::ArcsNotApproximated);
        }
        if !self.segs.is_empty() {
            self.check_segs();
        }
        let has_segs = !self.segs.is_empty();
        if self.restart {
            self.clean_data();
            self.restart = false;
        }
        if !self.bbox_added {
            self.add_bbox();
        }
        self.setup_grid();

        self.init_vd();
        let max_n_pnts = self.pnts.len() - NUM_BBOX_RGT_PNTS;
        debug_assert!(max_n_pnts >= FIRST_INSERTED_PNT, "max_n_pnts >= 3");
        for i in FIRST_INSERTED_PNT..max_n_pnts {
            self.insert_pnt_into_vd(i);
            if self.restart {
                return Ok(());
            }
        }
        self.free_grid_buffers();
        self.free_insert_buffers();
        self.free_reverse_index();
        if has_segs {
            self.nodes_r2_to_r();
            self.insert_pnt_pnt_mid_nodes();
        }

        let mut i = 0;
        while i < self.segs.len() {
            self.insert_seg_into_vd(i);
            if self.restart {
                return Ok(());
            }
            i += 1;
        }

        self.insert_seg_endpoint_nodes();
        if self.segs.is_empty() {
            self.isolated_pnts = true;
        } else if !self.data_cleaned {
            self.find_isolated_pnts();
        }
        if recover_arcs {
            self.recover_arc_vd();
        } else {
            self.discard_apx_arc_refs();
        }
        self.vd_computed = true;
        Ok(())
    }
}
