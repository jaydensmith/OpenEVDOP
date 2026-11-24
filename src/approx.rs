use crate::consts::*;
use crate::state::*;
use std::f64::consts::TAU;

impl Engine {
    /// Replaces every arc by a circumscribed polygonal chain whose deviation from the arc is at
    /// most `apx_abs` (in input units); the arcs stay addressable as `num_apx_arcs` records.
    pub(crate) fn approx_arcs(&mut self, apx_abs: f64) {
        self.check_arcs();
        debug_assert!(
            apx_abs > 0.0,
            "(apx_absolute > 0.0) || (apx_relative > 0.0)"
        );
        let mut apx_abs = apx_abs;
        if self.data_scaled {
            debug_assert!(self.scale_factor > ZERO, "gt(scale_factor, ZERO)");
            apx_abs *= self.scale_factor;
        }

        let mut vertices: Vec<(f64, f64)> = Vec::new();
        for i in 0..self.num_arcs {
            let arc = self.arcs[i];
            let p1 = self.pnts[arc.i1];
            let p2 = self.pnts[arc.i2];
            let (cx, cy) = (arc.cx, arc.cy);
            let r = arc.r;
            debug_assert!(r > ZERO, "gt(radius, ZERO)");
            let (start, end) = if arc.ccw {
                (arc.i1, arc.i2)
            } else {
                (arc.i2, arc.i1)
            };

            let mut a1 = (p1.y - cy).atan2(p1.x - cx);
            let mut a2 = (p2.y - cy).atan2(p2.x - cx);
            if 0.0 > a1 {
                a1 += TAU;
            }
            if 0.0 > a2 {
                a2 += TAU;
            }
            if a1 > a2 {
                a2 += TAU;
            }

            let max_incr = (r / (r + apx_abs)).acos();
            debug_assert!(max_incr > ZERO, "gt(max_incr, ZERO)");
            let sweep = a2 - a1;
            let pieces = sweep / (max_incr + max_incr);
            let mut number = pieces as i32;
            if f64::from(number) != pieces {
                number += 1;
            }
            debug_assert!(number > 0, "number > 0");
            let mut incr = sweep / f64::from(number);
            let half = incr * 0.5;
            let radius = r / half.cos();
            let mut angle = if arc.ccw {
                half + a1
            } else {
                incr = -incr;
                a2 - half
            };

            vertices.clear();
            for _ in 0..number {
                vertices.push((angle.cos() * radius + cx, angle.sin() * radius + cy));
                angle += incr;
            }

            let mut prev = start;
            for &(x, y) in &vertices {
                let v = self.store_pnt(x, y);
                self.pnts[v].apx_arc = i;
                if x <= self.bb_max_x {
                    if x < self.bb_min_x {
                        self.bb_min_x = x;
                    }
                } else {
                    self.bb_max_x = x;
                }
                if y <= self.bb_max_y {
                    if y < self.bb_min_y {
                        self.bb_min_y = y;
                    }
                } else {
                    self.bb_max_y = y;
                }
                let s = self.store_seg(prev, v);
                self.segs[s].apx_arc = i;
                prev = v;
            }
            let s = self.store_seg(prev, end);
            self.segs[s].apx_arc = i;
        }

        self.num_apx_arcs = self.num_arcs;
        self.num_arcs = 0;
    }
}
