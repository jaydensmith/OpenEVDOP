use crate::consts::*;
use crate::state::*;

const TREE_SEED_MAX_DEPTH: i32 = 2;

impl Engine {
    pub(crate) fn heap_reset(&mut self) {
        self.heap.clear();
        self.heap_top_pending = false;
        self.heap_top_valid = false;
    }

    /// Inserts `i` as a quadtree leaf unless duplicate; a heap reset limits the search to the descent path.
    pub(crate) fn find_nearest_pnt(
        &mut self,
        i: usize,
        reset_heap_after_descent: bool,
    ) -> (usize, f64) {
        let px = self.pnts[i].x;
        let py = self.pnts[i].y;
        debug_assert!(px > self.tree_xmin, "p.x > tree_min_x");
        debug_assert!(py > self.tree_ymin, "p.y > tree_min_y");
        debug_assert!(px < self.tree_xmax, "p.x < tree_max_x");
        debug_assert!(py < self.tree_ymax, "p.y < tree_max_y");
        debug_assert!(self.tree_root == 0, "next == 0");

        let mut min_ind = NIL;
        let mut min_dist = f64::MAX;
        let mut root_insert = NIL;
        let mut child_insert = NIL;
        let mut first_descent = true;
        let mut next = self.tree_root;
        self.heap_reset();
        loop {
            let mut node;
            let mut quadrant;
            loop {
                node = next;
                let split = self.tree[node].pnt;
                let dx = px - self.pnts[split].x;
                let dy = py - self.pnts[split].y;
                let dx2 = dx * dx;
                let dy2 = dy * dy;
                let dist = dy2 + dx2;
                if dist < min_dist {
                    min_ind = split;
                    min_dist = dist;
                }
                let child = self.tree[node].child;
                let (across_y, across_x, down) = if dx >= 0.0 {
                    if dy >= 0.0 { (2, 1, 3) } else { (3, 0, 2) }
                } else if dy >= 0.0 {
                    (0, 3, 1)
                } else {
                    (1, 2, 0)
                };
                if child[across_y] != NIL && dy2 <= min_dist {
                    self.heap_push(dy2, child[across_y]);
                }
                if child[across_x] != NIL && dx2 <= min_dist {
                    self.heap_push(dx2, child[across_x]);
                }
                next = child[down];
                quadrant = down;
                if next == NIL {
                    break;
                }
            }
            if first_descent {
                first_descent = false;
                root_insert = node;
                child_insert = quadrant;
                if reset_heap_after_descent {
                    self.heap_reset();
                }
            }
            while next == NIL {
                match self.heap_pop_min() {
                    Some((key, candidate)) => {
                        if key <= min_dist {
                            next = candidate;
                        }
                    }
                    None => break,
                }
            }
            if next == NIL {
                break;
            }
        }

        if min_dist > ZERO_SQ {
            debug_assert!(root_insert != NIL, "root_insert != NIL");
            let leaf = self.tree.len();
            self.tree.push(TreeNode {
                pnt: i,
                child: [NIL; 4],
            });
            self.tree[root_insert].child[child_insert] = leaf;
        }
        debug_assert!(min_ind != NIL, "min_ind != NIL");
        (min_ind, min_dist)
    }

    /// Seeds the tree with points nearest the centres of recursive box splits, swap-removing them from `list`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn build_tree_rec(
        &mut self,
        mut xmin: f64,
        mut ymin: f64,
        xmax: f64,
        ymax: f64,
        list: &mut [usize],
        count: &mut usize,
        mut depth: i32,
    ) {
        loop {
            let depth1 = depth + 1;
            let cx = (xmin + xmax) * 0.5;
            let cy = (ymin + ymax) * 0.5;
            if *count < 1 {
                return;
            }
            let mut best = NIL;
            let mut best_dist = f64::MAX;
            for (n, &pi) in list[..*count].iter().enumerate() {
                let p = self.pnts[pi];
                if xmin <= p.x && p.x <= xmax && ymin <= p.y && p.y <= ymax {
                    let dist = (cx - p.x) * (cx - p.x) + (cy - p.y) * (cy - p.y);
                    if dist < best_dist {
                        best = n;
                        best_dist = dist;
                    }
                }
            }
            if best == NIL {
                return;
            }
            let pnt = list[best];
            if self.tree_root == NIL {
                self.tree_root = self.tree.len();
                self.tree.push(TreeNode {
                    pnt,
                    child: [NIL; 4],
                });
                debug_assert!(self.tree_root == 0, "root_of_tree == 0");
            } else {
                self.find_nearest_pnt(pnt, true);
            }
            let px = self.pnts[pnt].x;
            let py = self.pnts[pnt].y;
            *count -= 1;
            list[best] = list[*count];
            if depth1 > TREE_SEED_MAX_DEPTH {
                return;
            }
            self.build_tree_rec(xmin, ymin, px, py, list, count, depth1);
            self.build_tree_rec(xmin, py, px, ymax, list, count, depth1);
            self.build_tree_rec(px, ymin, xmax, py, list, count, depth1);
            xmin = px;
            ymin = py;
            depth = depth1;
        }
    }

    /// Returns the nearest neighbour and squared distance of `grid_first_n_pnts[n_initial]`.
    pub(crate) fn build_tree(&mut self, n_initial: usize) -> (usize, f64) {
        self.tree.clear();
        self.tree_root = NIL;
        let mut margin_x = (self.bb_max_x - self.bb_min_x) * BBOX_MARGIN_FRACTION;
        let mut margin_y = (self.bb_max_y - self.bb_min_y) * BBOX_MARGIN_FRACTION;
        if margin_x <= GRID_MIN_MARGIN {
            margin_x = GRID_MIN_MARGIN;
        }
        if margin_y <= GRID_MIN_MARGIN {
            margin_y = GRID_MIN_MARGIN;
        }
        self.tree_xmin = self.bb_min_x - margin_x;
        self.tree_xmax = self.bb_max_x + margin_x;
        self.tree_ymin = self.bb_min_y - margin_y;
        self.tree_ymax = self.bb_max_y + margin_y;
        debug_assert!(self.tree_xmax - self.tree_xmin > ZERO, "gt(delta_x, ZERO)");
        debug_assert!(self.tree_ymax - self.tree_ymin > ZERO, "gt(delta_y, ZERO)");
        debug_assert!(n_initial >= 1, "N_Initial >= 1");

        let mut list = std::mem::take(&mut self.grid_first_n_pnts);
        let mut count = n_initial;
        let (xmin, ymin, xmax, ymax) = (
            self.tree_xmin,
            self.tree_ymin,
            self.tree_xmax,
            self.tree_ymax,
        );
        self.build_tree_rec(xmin, ymin, xmax, ymax, &mut list, &mut count, 0);
        debug_assert!(self.tree_root != NIL, "root_of_tree != NIL");
        for &pnt in &list[..count] {
            self.find_nearest_pnt(pnt, true);
        }
        let query = list[n_initial];
        self.grid_first_n_pnts = list;
        self.find_nearest_pnt(query, false)
    }
}
