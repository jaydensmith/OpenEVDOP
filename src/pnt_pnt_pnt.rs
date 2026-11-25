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

/// Intersection of the perpendicular bisectors of edges `base` and `other`,
/// parametrised from the midpoint of `base`.
fn bisector_intersection<R: Real>(
    base_mid: Point,
    base_u: Point,
    other_mid: Point,
    other_u: Point,
) -> Option<Point> {
    let a = [[-base_u.y, other_u.y], [base_u.x, -other_u.x]];
    let rhs = [
        (R::of(other_mid.x) - base_mid.x).stored(),
        (R::of(other_mid.y) - base_mid.y).stored(),
    ];
    let (x, _) = solve_2x2::<R>(a, rhs)?;
    let s = x[0];
    Some(Point {
        x: (R::of(-base_u.y) * s + base_mid.x).stored(),
        y: (R::of(s) * base_u.x + base_mid.y).stored(),
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
    let ca = (R::of(ca_x) * ca_x + R::of(ca_y) * ca_y).sqrt().stored();
    debug_assert!(ca >= 0.0, "ca >= 0.0");
    if ca <= ZERO {
        let (center, r2) = degenerate_3pts_cntr::<R>(c, a, b);
        return (false, center, r2);
    }

    let u_ab = Point {
        x: (R::of(ab_x) / ab).stored(),
        y: (R::of(ab_y) / ab).stored(),
    };
    let u_bc = Point {
        x: (R::of(bc_x) / bc).stored(),
        y: (R::of(bc_y) / bc).stored(),
    };
    let u_ca = Point {
        x: (R::of(ca_x) / ca).stored(),
        y: (R::of(ca_y) / ca).stored(),
    };
    let cos_b = (R::of(u_bc.x) * u_ab.x + R::of(u_bc.y) * u_ab.y)
        .abs()
        .stored();
    let cos_c = (R::of(u_ca.x) * u_bc.x + R::of(u_ca.y) * u_bc.y)
        .abs()
        .stored();
    let cos_a = (R::of(u_ca.x) * u_ab.x + R::of(u_ab.y) * u_ca.y)
        .abs()
        .stored();
    let mut min_cos = cos_a;
    if cos_c <= cos_b {
        if cos_c < min_cos {
            min_cos = cos_c;
        }
    } else if cos_b < min_cos {
        min_cos = cos_b;
    }

    let m_ab = Point {
        x: ((R::of(a.x) + b.x) * 0.5).stored(),
        y: ((R::of(a.y) + b.y) * 0.5).stored(),
    };
    let m_bc = Point {
        x: ((R::of(c.x) + b.x) * 0.5).stored(),
        y: ((R::of(c.y) + b.y) * 0.5).stored(),
    };
    let m_ca = Point {
        x: ((R::of(c.x) + a.x) * 0.5).stored(),
        y: ((R::of(c.y) + a.y) * 0.5).stored(),
    };
    let solved = if cos_b == min_cos {
        if ca <= ab {
            bisector_intersection::<R>(m_ca, u_ca, m_ab, u_ab)
        } else {
            bisector_intersection::<R>(m_ab, u_ab, m_ca, u_ca)
        }
    } else if cos_c == min_cos {
        if bc <= ab {
            bisector_intersection::<R>(m_bc, u_bc, m_ab, u_ab)
        } else {
            bisector_intersection::<R>(m_ab, u_ab, m_bc, u_bc)
        }
    } else if bc <= ca {
        bisector_intersection::<R>(m_bc, u_bc, m_ca, u_ca)
    } else {
        bisector_intersection::<R>(m_ca, u_ca, m_bc, u_bc)
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
    let den = (R::of(dy) * u_ab.x + R::of(dx) * -u_ab.y).stored();
    let den = (R::of(den) + den).stored();
    if den > ZERO || den < -ZERO {
        let h2 = (R::of(ab) * ab * 0.25).stored();
        let t = ((R::of(h2) - (R::of(dy) * dy + R::of(dx) * dx)) / den).stored();
        let center = Point {
            x: (R::of(mx) + R::of(-u_ab.y) * t).stored(),
            y: (R::of(my) + R::of(t) * u_ab.x).stored(),
        };
        return (true, center, (R::of(t) * t + h2).stored());
    }
    let center = Point {
        x: ((R::of(sum_x) + c.x) / 3.0).stored(),
        y: ((R::of(sum_y) + c.y) / 3.0).stored(),
    };
    (false, center, dist2::<R>(a, center).stored())
}
