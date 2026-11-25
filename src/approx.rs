use crate::consts::*;
use crate::real::Real;
use crate::state::*;
use std::f64::consts::TAU;
use x87_extended::double::{acos, cos, sin};

impl<R: Real> Engine<R> {
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
            apx_abs = (R::of(apx_abs) * self.scale_factor).stored();
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

            let angle_of = |p: Pnt| {
                let (y, x) = ((R::of(p.y) - cy).stored(), (R::of(p.x) - cx).stored());
                R::of(y).atan2(R::of(x)).stored()
            };
            let mut a1 = angle_of(p1);
            let mut a2 = angle_of(p2);
            if 0.0 > a1 {
                a1 = (R::of(a1) + TAU).stored();
            }
            if 0.0 > a2 {
                a2 = (R::of(a2) + TAU).stored();
            }
            if a1 > a2 {
                a2 = (R::of(a2) + TAU).stored();
            }

            let max_incr = acos((R::of(r) / (R::of(r) + apx_abs)).stored());
            debug_assert!(max_incr > ZERO, "gt(max_incr, ZERO)");
            let sweep = (R::of(a2) - a1).stored();
            let pieces = (R::of(sweep) / (R::of(max_incr) + max_incr)).stored();
            let mut number = R::of(pieces).truncate();
            if f64::from(number) != pieces {
                number += 1;
            }
            debug_assert!(number > 0, "number > 0");
            let mut incr = (R::of(sweep) / f64::from(number)).stored();
            let half = (R::of(incr) * 0.5).stored();
            let radius = (R::of(r) / cos(half)).stored();
            let mut angle = if arc.ccw {
                (R::of(half) + a1).stored()
            } else {
                incr = -incr;
                (R::of(a2) - half).stored()
            };

            vertices.clear();
            for _ in 0..number {
                vertices.push((
                    (R::of(cos(angle)) * radius + cx).stored(),
                    (R::of(sin(angle)) * radius + cy).stored(),
                ));
                angle = (R::of(angle) + incr).stored();
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
