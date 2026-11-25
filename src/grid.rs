use crate::consts::*;
use crate::real::Real;
use crate::state::*;

const GRID_SAMPLE_MIN: usize = 1500;
const GRID_SAMPLE_DOUBLINGS: usize = 3;

impl<R: Real> Engine<R> {
    pub(crate) fn free_grid_buffers(&mut self) {
        self.grid.clear();
        self.grid_x_coords.clear();
        self.grid_y_coords.clear();
        self.grid_pnt_list.clear();
        self.grid_first_n_pnts.clear();
        self.grid_n = 0;
        self.grid_nx = 0;
        self.grid_ny = 0;
        self.grid_first_counter = 0;
    }

    fn init_grid_cells(&mut self) {
        debug_assert!(self.grid_n > 0, "N > 0");
        debug_assert!(self.grid_nx * self.grid_ny <= self.grid_n, "N_x * N_y <= N");
        self.grid.clear();
        self.grid.resize(self.grid_n, NIL);
    }

    fn init_grid_coords(&mut self) {
        let (nx, ny) = (self.grid_nx, self.grid_ny);
        self.grid_x_coords.clear();
        self.grid_x_coords.extend(
            (0..nx).map(|i| (R::of(i as f64) * self.grid_cell_w + self.grid_min_x).stored()),
        );
        self.grid_x_coords.push(self.grid_max_x);
        self.grid_y_coords.clear();
        self.grid_y_coords.extend(
            (0..ny).map(|j| (R::of(j as f64) * self.grid_cell_h + self.grid_min_y).stored()),
        );
        self.grid_y_coords.push(self.grid_max_y);
    }

    fn init_grid(&mut self, number: usize) {
        let s = R::of(number as f64).sqrt().stored();
        self.grid_nx = (R::of(self.grid_aspect) * s).truncate().max(1) as usize;
        self.grid_ny = (R::of(s) * self.grid_inv_aspect).truncate().max(1) as usize;
        self.grid_n = self.grid_ny * self.grid_nx;
        self.grid_cell_w = (R::of(self.grid_delta_x) / self.grid_nx as f64).stored();
        self.grid_cell_h = (R::of(self.grid_delta_y) / self.grid_ny as f64).stored();
        self.init_grid_cells();
        self.init_grid_coords();
        self.grid_pnt_list.clear();
        self.grid_pnt_list.reserve(self.pnts.len());
        self.grid_active = true;
    }

    pub(crate) fn setup_grid(&mut self) {
        let mut margin_x = ((R::of(self.bb_max_x) - self.bb_min_x) * BBOX_MARGIN_FRACTION).stored();
        let mut margin_y = ((R::of(self.bb_max_y) - self.bb_min_y) * BBOX_MARGIN_FRACTION).stored();
        if margin_x <= GRID_MIN_MARGIN {
            margin_x = GRID_MIN_MARGIN;
        }
        if margin_y <= GRID_MIN_MARGIN {
            margin_y = GRID_MIN_MARGIN;
        }
        self.grid_min_x = (R::of(self.bb_min_x) - margin_x).stored();
        self.grid_max_x = (R::of(self.bb_max_x) + margin_x).stored();
        self.grid_min_y = (R::of(self.bb_min_y) - margin_y).stored();
        self.grid_max_y = (R::of(self.bb_max_y) + margin_y).stored();
        self.grid_delta_x = (R::of(self.grid_max_x) - self.grid_min_x).stored();
        self.grid_delta_y = (R::of(self.grid_max_y) - self.grid_min_y).stored();
        debug_assert!(self.grid_delta_x > ZERO, "gt(delta_x, ZERO)");
        debug_assert!(self.grid_delta_y > ZERO, "gt(delta_y, ZERO)");

        self.grid_aspect = (R::of(self.grid_delta_x) / self.grid_delta_y).stored();
        if GRID_MIN_ASPECT <= self.grid_aspect {
            if self.grid_aspect <= GRID_MAX_ASPECT {
                self.grid_inv_aspect = (R::of(1.0) / self.grid_aspect).stored();
            } else {
                self.grid_aspect = GRID_MAX_ASPECT;
                self.grid_inv_aspect = GRID_MIN_ASPECT;
            }
        } else {
            self.grid_aspect = GRID_MIN_ASPECT;
            self.grid_inv_aspect = GRID_MAX_ASPECT;
        }

        let mut n = (R::of(self.pnts.len() as f64) * GRID_SAMPLE_FRACTION).truncate() as usize;
        let mut doublings = 0;
        while n < GRID_SAMPLE_MIN && doublings < GRID_SAMPLE_DOUBLINGS {
            n += n;
            doublings += 1;
        }
        self.grid_first_n_pnts.clear();
        if n >= GRID_SAMPLE_MIN {
            self.grid_n_initial = n;
            self.grid_first_n_pnts.reserve(n + 2);
        } else {
            self.grid_n_initial = NIL;
            n = self.pnts.len();
        }
        self.init_grid(2 * n);
        self.grid_first_counter = 0;
    }

    fn grid_cell_of(&self, x: f64, y: f64) -> usize {
        debug_assert!(x > self.grid_min_x, "p.x > grid_min_x");
        let mut i = ((R::of(x) - self.grid_min_x) / self.grid_cell_w).truncate() as isize;
        if self.grid_x_coords[i as usize] > x {
            i -= 1;
        } else if self.grid_x_coords[i as usize + 1] < x {
            i += 1;
        }
        debug_assert!(
            i >= 0 && (i as usize) < self.grid_nx,
            "(i >= 0) && (i < N_x)"
        );
        debug_assert!(y > self.grid_min_y, "p.y > grid_min_y");
        let mut j = ((R::of(y) - self.grid_min_y) / self.grid_cell_h).truncate() as isize;
        if self.grid_y_coords[j as usize] > y {
            j -= 1;
        } else if self.grid_y_coords[j as usize + 1] < y {
            j += 1;
        }
        debug_assert!(
            j >= 0 && (j as usize) < self.grid_ny,
            "(j >= 0) && (j < N_y)"
        );
        self.grid_ny * i as usize + j as usize
    }

    fn grid_cell_indices(&self, x: f64, y: f64) -> (usize, usize) {
        let k = self.grid_cell_of(x, y);
        (k / self.grid_ny, k % self.grid_ny)
    }

    fn push_grid_entry(&mut self, index: usize, k: usize) {
        self.grid_pnt_list.push(GridEntry {
            pnt: index,
            next: self.grid[k],
        });
        self.grid[k] = self.grid_pnt_list.len() - 1;
    }

    fn scan_grid_cell(&self, k: usize, x: f64, y: f64, best: &mut f64, nearest: &mut usize) {
        let mut e = self.grid[k];
        while e != NIL {
            let q = self.grid_pnt_list[e].pnt;
            let dx = (R::of(x) - self.pnts[q].x).stored();
            let dy = (R::of(y) - self.pnts[q].y).stored();
            let d = (R::of(dy) * dy + R::of(dx) * dx).stored();
            if d < *best {
                *best = d;
                *nearest = q;
            }
            e = self.grid_pnt_list[e].next;
        }
    }

    pub(crate) fn delete_pnt_from_grid(&mut self, index: usize) {
        let k = self.grid_cell_of(self.pnts[index].x, self.pnts[index].y);
        debug_assert!(
            self.grid_pnt_list[self.grid[k]].pnt == index,
            "pnt_list[grid[k]].pnt == index"
        );
        self.grid[k] = self.grid_pnt_list[self.grid[k]].next;
    }

    /// Returns the nearest grid point in the cell of `index` and its squared distance, then
    /// inserts `index` unless it coincides with that point.
    pub(crate) fn insert_pnt_into_grid(&mut self, index: usize) -> (usize, f64) {
        let (x, y) = (self.pnts[index].x, self.pnts[index].y);
        let k = self.grid_cell_of(x, y);
        let mut best = f64::MAX;
        let mut nearest = NIL;
        self.scan_grid_cell(k, x, y, &mut best, &mut nearest);
        if best <= ZERO_SQ {
            return (nearest, best);
        }
        self.push_grid_entry(index, k);

        let n_initial = self.grid_n_initial;
        if n_initial == NIL || self.grid_first_counter > n_initial {
            return (nearest, best);
        }
        self.grid_first_n_pnts.push(index);
        if self.grid_first_counter < n_initial {
            self.grid_first_counter += 1;
            return (nearest, best);
        }

        let num = self.grid.iter().filter(|&&head| head != NIL).count();
        debug_assert!(num > 0, "num > 0");
        let ratio = (R::of(n_initial as f64) / num as f64).stored();
        if ratio <= GRID_REGRID_RATIO {
            let density = if GRID_MIN_DENSITY > ratio {
                GRID_MIN_DENSITY
            } else if GRID_MAX_DENSITY >= ratio {
                ratio
            } else {
                GRID_MAX_DENSITY
            };
            self.init_grid((R::of(density) * self.pnts.len() as f64).truncate() as usize);
            for m in 0..=n_initial {
                let pm = self.grid_first_n_pnts[m];
                let km = self.grid_cell_of(self.pnts[pm].x, self.pnts[pm].y);
                self.push_grid_entry(pm, km);
            }
        } else {
            self.grid_active = false;
            (nearest, best) = self.build_tree(n_initial);
            self.free_grid_buffers();
        }
        self.grid_first_counter += 1;
        (nearest, best)
    }

    /// Expanding-ring search around the cell of `index`, starting from `candidate` at squared
    /// distance `dist_sq`; returns the nearest point and its squared distance.
    pub(crate) fn find_nearest_pnt_in_grid(
        &self,
        index: usize,
        candidate: usize,
        dist_sq: f64,
    ) -> (usize, f64) {
        let (x, y) = (self.pnts[index].x, self.pnts[index].y);
        let (i, j) = self.grid_cell_indices(x, y);
        let mut best = dist_sq;
        let mut nearest = candidate;
        let (mut i_lo, mut i_hi, mut j_lo, mut j_hi) = (i, i, j, j);
        let (mut lft, mut rgt, mut bot, mut top) = (true, true, true, true);
        loop {
            if lft || rgt {
                if lft {
                    if i_lo < 1 || sq::<R>((R::of(x) - self.grid_x_coords[i_lo]).stored()) >= best {
                        lft = false;
                    } else {
                        i_lo -= 1;
                        for jj in j_lo..=j_hi {
                            self.scan_grid_cell(
                                self.grid_ny * i_lo + jj,
                                x,
                                y,
                                &mut best,
                                &mut nearest,
                            );
                        }
                    }
                }
                if rgt {
                    if i_hi + 1 >= self.grid_nx
                        || sq::<R>((R::of(x) - self.grid_x_coords[i_hi + 1]).stored()) >= best
                    {
                        rgt = false;
                    } else {
                        i_hi += 1;
                        for jj in j_lo..=j_hi {
                            self.scan_grid_cell(
                                self.grid_ny * i_hi + jj,
                                x,
                                y,
                                &mut best,
                                &mut nearest,
                            );
                        }
                    }
                }
            } else if !top && !bot {
                debug_assert!(nearest != NIL, "i_min != NIL");
                debug_assert!(best < f64::MAX, "dist_min < DBL_MAX");
                return (nearest, best);
            }
            if bot {
                if j_lo < 1 || sq::<R>((R::of(y) - self.grid_y_coords[j_lo]).stored()) >= best {
                    bot = false;
                } else {
                    j_lo -= 1;
                    for ii in i_lo..=i_hi {
                        self.scan_grid_cell(
                            self.grid_ny * ii + j_lo,
                            x,
                            y,
                            &mut best,
                            &mut nearest,
                        );
                    }
                }
            }
            if top {
                if j_hi + 1 >= self.grid_ny
                    || sq::<R>((R::of(y) - self.grid_y_coords[j_hi + 1]).stored()) >= best
                {
                    top = false;
                } else {
                    j_hi += 1;
                    for ii in i_lo..=i_hi {
                        self.scan_grid_cell(
                            self.grid_ny * ii + j_hi,
                            x,
                            y,
                            &mut best,
                            &mut nearest,
                        );
                    }
                }
            }
        }
    }
}

fn sq<R: Real>(v: f64) -> R {
    R::of(v) * v
}
