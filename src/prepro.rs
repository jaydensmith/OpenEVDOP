use crate::consts::*;
use crate::data::dist;
use crate::real::Real;
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

impl<R: Real> Engine<R> {
    /// Computes the normalized line equation and length of every segment; zero-length segments
    /// are removed (merging their endpoints if they differ).
    pub(crate) fn check_segs(&mut self) {
        let mut i = 0;
        while i < self.segs.len() {
            let (i1, i2) = (self.segs[i].i1, self.segs[i].i2);
            let p1 = self.pnts[i1];
            let p2 = self.pnts[i2];
            let ux = (R::of(p2.x) - p1.x).stored();
            let uy = (R::of(p2.y) - p1.y).stored();
            let lgth = (R::of(uy) * uy + R::of(ux) * ux).sqrt().stored();
            if ZERO < lgth {
                let ux = (R::of(ux) / lgth).stored();
                let uy = (R::of(uy) / lgth).stored();
                let mut c = (R::of(uy) * p1.x - R::of(p1.y) * ux).stored();
                let a = -uy;
                let b = ux;
                let d = (R::of(a) * p2.x + R::of(p2.y) * ux + c).stored();
                if ZERO < d || d < -ZERO {
                    c = ((R::of(uy) * p2.x - R::of(p2.y) * ux + c) * 0.5).stored();
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
            let x12 = (R::of(x1) - x2).stored();
            let y12 = (R::of(y1) - y2).stored();
            let orient = ((R::of(cx) - x1) * y12 + R::of(x12) * (R::of(y1) - cy)).stored();
            let (u, v) = if 0.0 > orient {
                (x12, y12)
            } else {
                ((R::of(x2) - x1).stored(), (R::of(y2) - y1).stored())
            };
            let s_orient = sign_zero(orient);
            let chord = (R::of(v) * v + R::of(u) * u).sqrt().stored();
            let half = (R::of(chord) * 0.5).stored();

            if !(ZERO < chord) {
                if x2 != x1 || y2 != y1 {
                    self.merge_pnts(arc.i1, arc.i2);
                }
                self.delete_arc(i);
                continue;
            }

            let d1 = dist::<R>((R::of(x1) - cx).stored(), (R::of(y1) - cy).stored()).stored();
            let d2 = dist::<R>((R::of(x2) - cx).stored(), (R::of(y2) - cy).stored()).stored();
            let mut r = ((R::of(d1) + d2) * 0.5).stored();
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

            let diff = (R::of(d1) - d2).stored();
            if SMALL < diff || diff < -SMALL {
                let mid_x = ((R::of(x2) + x1) * 0.5).stored();
                let mid_y = ((R::of(y2) + y1) * 0.5).stored();
                let w = (R::of(r) * r - R::of(half) * half).stored();
                let k = if 0.0 < w {
                    (R::of(w).sqrt() / chord).stored()
                } else {
                    0.0
                };
                let new_cx = (R::of(-v) * k + mid_x).stored();
                let new_cy = (R::of(k) * u + mid_y).stored();
                let orient2 =
                    ((R::of(new_cx) - x1) * y12 + R::of(x12) * (R::of(y1) - new_cy)).stored();
                debug_assert!(
                    !(ZERO < orient2 || orient2 < -ZERO)
                        || s_orient == 0
                        || sign_zero(orient2) == s_orient,
                    "(SignEps(orient, ZERO) == 0) || (s_orient == 0) || (SignEps(orient, ZERO) == s_orient)"
                );
                let d1 = dist::<R>((R::of(x1) - new_cx).stored(), (R::of(y1) - new_cy).stored())
                    .stored();
                let d2 = dist::<R>((R::of(x2) - new_cx).stored(), (R::of(y2) - new_cy).stored())
                    .stored();
                r = ((R::of(d2) + d1) * 0.5).stored();
                self.arcs[i].cx = new_cx;
                self.arcs[i].cy = new_cy;
            }
            debug_assert!(ZERO < r, "gt(radius, ZERO)");
            self.arcs[i].r = r;
            i += 1;
        }
    }
}
