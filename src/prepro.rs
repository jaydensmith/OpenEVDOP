use crate::consts::*;
use crate::data::dist;
use crate::state::*;

fn sign_zero(v: f64) -> i32 {
    if ZERO < v {
        1
    } else if v >= -ZERO {
        0
    } else {
        -1
    }
}

impl Engine {
    /// Computes the normalized line equation and length of every segment; zero-length segments
    /// are removed (merging their endpoints if they differ).
    pub(crate) fn check_segs(&mut self) {
        let mut i = 0;
        while i < self.segs.len() {
            let (i1, i2) = (self.segs[i].i1, self.segs[i].i2);
            let p1 = self.pnts[i1];
            let p2 = self.pnts[i2];
            let ux = p2.x - p1.x;
            let uy = p2.y - p1.y;
            let lgth = (uy * uy + ux * ux).sqrt();
            if ZERO < lgth {
                let ux = ux / lgth;
                let uy = uy / lgth;
                let mut c = uy * p1.x - p1.y * ux;
                let a = -uy;
                let b = ux;
                let d = a * p2.x + p2.y * ux + c;
                if ZERO < d || d < -ZERO {
                    c = ((uy * p2.x - p2.y * ux) + c) * 0.5;
                }
                let seg = &mut self.segs[i];
                seg.a = a;
                seg.b = b;
                seg.c = c;
                seg.lgth = lgth;
                i += 1;
            } else {
                if p2.x != p1.x || p2.y != p1.y {
                    self.merge_pnts(i1, i2);
                }
                self.delete_seg(i);
            }
        }
    }

    /// Validates every arc: degenerate chords merge their endpoints, inconsistent radii move the
    /// centre onto the chord bisector, and arcs with vanishing radius become segments.
    pub(crate) fn check_arcs(&mut self) {
        let mut i = 0;
        while i < self.num_arcs {
            let arc = self.arcs[i];
            let (x1, y1) = (self.pnts[arc.i1].x, self.pnts[arc.i1].y);
            let (x2, y2) = (self.pnts[arc.i2].x, self.pnts[arc.i2].y);
            let (cx, cy) = (arc.cx, arc.cy);
            let x12 = x1 - x2;
            let y12 = y1 - y2;
            let orient = (cx - x1) * y12 + x12 * (y1 - cy);
            let (u, v) = if 0.0 > orient {
                (x12, y12)
            } else {
                (x2 - x1, y2 - y1)
            };
            let s_orient = sign_zero(orient);
            let chord = (v * v + u * u).sqrt();
            let half = chord * 0.5;

            if !(ZERO < chord) {
                if x2 != x1 || y2 != y1 {
                    self.merge_pnts(arc.i1, arc.i2);
                }
                self.delete_arc(i);
                continue;
            }

            let d1 = dist(x1 - cx, y1 - cy);
            let d2 = dist(x2 - cx, y2 - cy);
            let mut r = (d1 + d2) * 0.5;
            if !(half <= r) {
                r = half;
            }
            if !(ZERO < r) {
                if arc.ccw {
                    self.store_seg(arc.i1, arc.i2);
                } else {
                    self.store_seg(arc.i2, arc.i1);
                }
                self.delete_arc(i);
                i += 1;
                continue;
            }

            let diff = d1 - d2;
            if SMALL < diff || diff < -SMALL {
                let mid_x = (x2 + x1) * 0.5;
                let mid_y = (y2 + y1) * 0.5;
                let w = r * r - half * half;
                let k = if 0.0 < w { w.sqrt() / chord } else { 0.0 };
                let new_cx = -v * k + mid_x;
                let new_cy = k * u + mid_y;
                let orient2 = (new_cx - x1) * y12 + x12 * (y1 - new_cy);
                debug_assert!(
                    !(ZERO < orient2 || orient2 < -ZERO)
                        || s_orient == 0
                        || sign_zero(orient2) == s_orient,
                    "(SignEps(orient, ZERO) == 0) || (s_orient == 0) || (SignEps(orient, ZERO) == s_orient)"
                );
                let d1 = dist(x1 - new_cx, y1 - new_cy);
                let d2 = dist(x2 - new_cx, y2 - new_cy);
                r = (d2 + d1) * 0.5;
                self.arcs[i].cx = new_cx;
                self.arcs[i].cy = new_cy;
            }
            debug_assert!(ZERO < r, "gt(radius, ZERO)");
            self.arcs[i].r = r;
            i += 1;
        }
    }
}
