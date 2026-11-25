use crate::consts::*;
use crate::numerics::{dist2, solve_2x2};
use crate::real::Real;
use crate::state::*;

pub(crate) fn degenerate_3pts_cntr<R: Real>(a: Point, b: Point, c: Point) -> (Point, f64) {
    let center = Point {
        x: ((R::of(a.x) + b.x) * 0.25 + R::of(c.x) * 0.5).stored(),
        y: ((R::of(b.y) + a.y) * 0.25 + R::of(c.y) * 0.5).stored(),
    };
    (center, dist2::<R>(a, center).stored())
}

pub(crate) fn pnt_pnt_cntr<R: Real>(a: Point, b: Point) -> (Point, f64) {
    let center = Point {
        x: ((R::of(b.x) + a.x) * 0.5).stored(),
        y: ((R::of(b.y) + a.y) * 0.5).stored(),
    };
    (center, dist2::<R>(a, center).stored())
}

/// An edge's midpoint and unit direction; the midpoint's x and the direction's y may still be
/// in the registers.
#[derive(Clone, Copy)]
struct Edge<R> {
    mid_x: R,
    mid_y: f64,
    u_x: f64,
    u_y: R,
}

/// Intersection of the perpendicular bisectors of edges `base` and `other`,
/// parametrised from the midpoint of `base`.
fn bisector_intersection<R: Real>(base: Edge<R>, other: Edge<R>) -> Option<Point> {
    let a = [
        [(-base.u_y).stored(), other.u_y.stored()],
        [base.u_x, -other.u_x],
    ];
    let rhs = [
        (other.mid_x - base.mid_x).stored(),
        (R::of(other.mid_y) - base.mid_y).stored(),
    ];
    let (x, _) = solve_2x2::<R>(a, rhs)?;
    let s = x[0];
    Some(Point {
        x: (base.mid_x + -base.u_y * s).stored(),
        y: (R::of(s) * base.u_x + base.mid_y).stored(),
    })
}

impl<R: Real> Engine<R> {
    pub(crate) fn pnt_pnt_pnt_node(&self, i: usize, j: usize, k: usize) -> (bool, Point, f64) {
        let mut idx = [i, j, k];
        idx.sort_unstable();
        let [i, j, k] = idx;
        let pi = self.pnt_point(i);
        if i == j || j == k {
            if i == k {
                return (true, pi, 0.0);
            }
            let (center, r2) = pnt_pnt_cntr::<R>(pi, self.pnt_point(k));
            return (true, center, r2);
        }
        let arc = self.pnts[i].apx_arc;
        if arc != NIL && arc == self.pnts[j].apx_arc && arc == self.pnts[k].apx_arc {
            let center = Point {
                x: self.arcs[arc].cx,
                y: self.arcs[arc].cy,
            };
            return (true, center, dist2::<R>(center, pi).stored());
        }
        circumcircle_3pts::<R>(pi, self.pnt_point(j), self.pnt_point(k))
    }

    pub(crate) fn pnt_point(&self, i: usize) -> Point {
        Point {
            x: self.pnts[i].x,
            y: self.pnts[i].y,
        }
    }
}

pub(crate) fn circumcircle_3pts<R: Real>(a: Point, b: Point, c: Point) -> (bool, Point, f64) {
    let (ab_x, ab_y) = ((R::of(b.x) - a.x).stored(), (R::of(b.y) - a.y).stored());
    let ab = (R::of(ab_x) * ab_x + R::of(ab_y) * ab_y).sqrt().stored();
    debug_assert!(ab >= 0.0, "ab >= 0.0");
    if ab <= ZERO {
        let (center, r2) = degenerate_3pts_cntr::<R>(a, b, c);
        return (false, center, r2);
    }
    let (bc_x, bc_y) = ((R::of(c.x) - b.x).stored(), (R::of(c.y) - b.y).stored());
    let bc = (R::of(bc_x) * bc_x + R::of(bc_y) * bc_y).sqrt().stored();
    debug_assert!(bc >= 0.0, "bc >= 0.0");
    if bc <= ZERO {
        let (center, r2) = degenerate_3pts_cntr::<R>(b, c, a);
        return (false, center, r2);
    }
    let (ca_x, ca_y) = ((R::of(a.x) - c.x).stored(), (R::of(a.y) - c.y).stored());
    let ca = (R::of(ca_y) * ca_y + R::of(ca_x) * ca_x).sqrt();
    debug_assert!(ca >= 0.0, "ca >= 0.0");
    if ca <= ZERO {
        let (center, r2) = degenerate_3pts_cntr::<R>(c, a, b);
        return (false, center, r2);
    }

    let (u_ab_x, u_ab_y) = ((R::of(ab_x) / ab).stored(), (R::of(ab_y) / ab).stored());
    let (u_bc_x, u_bc_y) = ((R::of(bc_x) / bc).stored(), R::of(bc_y) / bc);
    let (u_ca_x, u_ca_y) = ((R::of(ca_x) / ca).stored(), R::of(ca_y) / ca);
    let cos_b = (u_bc_y * u_ab_y + R::of(u_bc_x) * u_ab_x).stored().abs();
    let cos_c = (u_ca_y * u_bc_y + R::of(u_ca_x) * u_bc_x).stored().abs();
    let cos_a = (R::of(u_ab_y) * u_ca_y + R::of(u_ca_x) * u_ab_x).abs();
    let mut min_cos = cos_a;
    if !(cos_c > cos_b) {
        if !(min_cos <= cos_c) {
            min_cos = R::of(cos_c);
        }
    } else if !(min_cos <= cos_b) {
        min_cos = R::of(cos_b);
    }

    let midpoint = |p: Point, q: Point| ((R::of(p.x) + q.x) * 0.5, (R::of(p.y) + q.y) * 0.5);
    let edge = |(mid_x, mid_y): (R, R), u_x: f64, u_y: R| Edge {
        mid_x,
        mid_y: mid_y.stored(),
        u_x,
        u_y,
    };
    let stored = |(mid_x, mid_y): (R, R)| (R::of(mid_x.stored()), mid_y);
    let solved = if min_cos == cos_b {
        let e_ab = edge(midpoint(a, b), u_ab_x, R::of(u_ab_y));
        let e_ca = edge(stored(midpoint(c, a)), u_ca_x, u_ca_y);
        if !(R::of(ab) < ca) {
            bisector_intersection(e_ca, e_ab)
        } else {
            bisector_intersection(e_ab, e_ca)
        }
    } else if min_cos == cos_c {
        let e_ab = edge(stored(midpoint(a, b)), u_ab_x, R::of(u_ab_y));
        let e_bc = edge(stored(midpoint(c, b)), u_bc_x, u_bc_y);
        if !(ab < bc) {
            bisector_intersection(e_bc, e_ab)
        } else {
            bisector_intersection(e_ab, e_bc)
        }
    } else {
        let e_ca = edge(midpoint(c, a), u_ca_x, u_ca_y);
        let e_bc = edge(stored(midpoint(c, b)), u_bc_x, u_bc_y);
        if !(ca < bc) {
            bisector_intersection(e_bc, e_ca)
        } else {
            bisector_intersection(e_ca, e_bc)
        }
    };
    if let Some(center) = solved {
        return (true, center, dist2::<R>(a, center).stored());
    }

    let sum_x = (R::of(a.x) + b.x).stored();
    let sum_y = (R::of(a.y) + b.y).stored();
    let mx = (R::of(sum_x) * 0.5).stored();
    let my = (R::of(sum_y) * 0.5).stored();
    let dx = (R::of(mx) - c.x).stored();
    let dy = (R::of(my) - c.y).stored();
    let den = (R::of(dy) * u_ab_x + R::of(dx) * -u_ab_y).stored();
    let den = (R::of(den) + den).stored();
    if den > ZERO || den < -ZERO {
        let h2 = (R::of(ab) * ab * 0.25).stored();
        let t = ((R::of(h2) - (R::of(dy) * dy + R::of(dx) * dx)) / den).stored();
        let center = Point {
            x: (R::of(mx) + R::of(-u_ab_y) * t).stored(),
            y: (R::of(my) + R::of(t) * u_ab_x).stored(),
        };
        return (true, center, (R::of(t) * t + h2).stored());
    }
    let center = Point {
        x: ((R::of(sum_x) + c.x) / 3.0).stored(),
        y: ((R::of(sum_y) + c.y) / 3.0).stored(),
    };
    (false, center, dist2::<R>(a, center).stored())
}
