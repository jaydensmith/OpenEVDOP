use crate::error::OffsetError;

pub(crate) const NIL: usize = usize::MAX;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SiteType {
    Seg = 0,
    Arc = 1,
    Pnt = 2,
    NoSite = 8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NodeStatus {
    Unchecked = 0,
    Visited = 1,
    Deleted = 2,
    Processed = 3,
    Dummy = 4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReverseType {
    IsoPnt = 0,
    SegStart = 1,
    SegEnd = 2,
    ArcStart = 3,
    ArcEnd = 4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OffsetKind {
    Line = 0,
    PntArc = 2,
    ApxArcInside = 6,
    ApxArcOutside = 7,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Point {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NodeCandidate {
    pub center: Point,
    pub r: f64,
    pub ambiguous: bool,
    pub coincident: usize,
}

impl NodeCandidate {
    pub fn new() -> Self {
        Self {
            center: Point { x: 0.0, y: 0.0 },
            r: 0.0,
            ambiguous: false,
            coincident: NIL,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Pnt {
    pub x: f64,
    pub y: f64,
    pub vd_edge: usize,
    pub deg2_node: usize,
    pub visited: bool,
    pub has_deg2: bool,
    pub deleted: bool,
    pub apx_arc: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Seg {
    pub i1: usize,
    pub i2: usize,
    pub vd_flag: bool,
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub lgth: f64,
    pub vd_edge: usize,
    pub apx_arc: usize,
}

impl Seg {
    pub fn signed_dist(&self, x: f64, y: f64) -> f64 {
        self.a * x + self.b * y + self.c
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Arc {
    pub i1: usize,
    pub i2: usize,
    pub cx: f64,
    pub cy: f64,
    pub r: f64,
    pub ccw: bool,
    pub vd_edge: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ReverseInd {
    pub isite: usize,
    pub site_type: ReverseType,
    pub next: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Node {
    pub x: f64,
    pub y: f64,
    /// Squared clearance while `nodes_store_r2` is set, plain clearance afterwards.
    pub r: f64,
    pub status: NodeStatus,
    pub edge: usize,
    pub site: usize,
    pub deg2: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Edge {
    pub n1: usize,
    pub n2: usize,
    pub lft: usize,
    pub rgt: usize,
    pub ltype: SiteType,
    pub rtype: SiteType,
    pub s_ccw: usize,
    pub s_cw: usize,
    pub e_ccw: usize,
    pub e_cw: usize,
    pub deleted: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EdgeData {
    pub pending: bool,
    pub px: f64,
    pub py: f64,
    pub d: f64,
    pub ccw: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct EdgeFlags {
    pub flag0: bool,
    pub available: bool,
    pub degenerate: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct GridEntry {
    pub pnt: usize,
    pub next: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct HeapEntry {
    pub key: f64,
    pub value: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TreeNode {
    pub pnt: usize,
    pub child: [usize; 4],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct DummyRec {
    pub e: usize,
    pub e1: usize,
    pub n: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CheckPnt {
    pub site: usize,
    pub stype: SiteType,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OffsetLoop {
    pub start: usize,
    pub end: usize,
    pub dist: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct OffsetElem {
    pub site: usize,
    pub kind: OffsetKind,
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct Engine {
    // input data
    pub pnts: Vec<Pnt>,
    pub segs: Vec<Seg>,
    pub arcs: Vec<Arc>,
    pub num_arcs: usize,
    pub num_apx_arcs: usize,
    pub scale_factor: f64,
    pub shift_x: f64,
    pub shift_y: f64,
    pub data_scaled: bool,
    pub first_pnt_stored: bool,
    pub bbox_added: bool,
    pub bbox_rgt_added: bool,
    pub data_cleaned: bool,
    pub bb_min_x: f64,
    pub bb_min_y: f64,
    pub bb_max_x: f64,
    pub bb_max_y: f64,
    pub min_bbox_pnt: usize,
    pub max_bbox_pnt: usize,
    pub max_offset: f64,
    pub large: f64,

    // input cleaning
    pub reverse_ind: Vec<ReverseInd>,
    pub reverse_index_exists: bool,

    // driver state
    pub restart: bool,
    pub restart_count: usize,
    pub failure: Option<OffsetError>,
    pub troubles: bool,
    pub desperate_used: bool,
    pub isolated_pnts: bool,

    // Voronoi diagram
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
    pub free_nodes: Vec<usize>,
    pub free_edges: Vec<usize>,
    pub dummy_nodes: Vec<DummyRec>,
    pub nodes_store_r2: bool,
    pub vd_computed: bool,
    pub deg2_nodes_added: bool,
    pub vd_id: i32,
    pub between_eps: f64,
    pub pnt_visited: Vec<usize>,
    pub nodes_checked: Vec<usize>,
    pub check_pnts: Vec<CheckPnt>,
    pub loop_stack: Vec<usize>,
    pub loop_stack_start: usize,

    // grid for nearest point queries
    pub grid_active: bool,
    pub grid: Vec<usize>,
    pub grid_x_coords: Vec<f64>,
    pub grid_y_coords: Vec<f64>,
    pub grid_pnt_list: Vec<GridEntry>,
    pub grid_first_n_pnts: Vec<usize>,
    pub grid_first_counter: usize,
    pub grid_n_initial: usize,
    pub grid_nx: usize,
    pub grid_ny: usize,
    pub grid_n: usize,
    pub grid_min_x: f64,
    pub grid_min_y: f64,
    pub grid_max_x: f64,
    pub grid_max_y: f64,
    pub grid_delta_x: f64,
    pub grid_delta_y: f64,
    pub grid_cell_w: f64,
    pub grid_cell_h: f64,
    pub grid_aspect: f64,
    pub grid_inv_aspect: f64,

    // quadtree for nearest point queries
    pub tree: Vec<TreeNode>,
    pub tree_root: usize,
    pub tree_xmin: f64,
    pub tree_ymin: f64,
    pub tree_xmax: f64,
    pub tree_ymax: f64,

    // heap
    pub heap: Vec<HeapEntry>,
    pub heap_top_pending: bool,
    pub heap_top_valid: bool,

    // offsets
    pub e_data: Vec<EdgeData>,
    pub e_data_vd_id: i32,
    pub edge_flags: Vec<EdgeFlags>,
    pub active_edges: Vec<usize>,
    pub offset_loops: Vec<OffsetLoop>,
    pub offset_elems: Vec<OffsetElem>,
    pub offset_t: f64,
    pub offset_dt: f64,
    pub offset_vd_id: i32,
    pub new_offset_needed: bool,
}
