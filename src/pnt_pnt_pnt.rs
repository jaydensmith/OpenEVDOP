use crate::consts::*;
use crate::numerics::{dist2, solve_2x2};
use crate::state::*;

pub(crate) fn degenerate_3pts_cntr(a: Point, b: Point, c: Point) -> (Point, f64) {
    let center = Point {
        x: (a.x + b.x) * 0.25 + c.x * 0.5,
        y: (b.y + a.y) * 0.25 + c.y * 0.5,
    };
    (center, dist2(a, center))
}

pub(crate) fn pnt_pnt_cntr(a: Point, b: Point) -> (Point, f64) {
    let center = Point {
        x: (b.x + a.x) * 0.5,
        y: (b.y + a.y) * 0.5,
    };
    (center, dist2(a, center))
}

/// Intersection of the perpendicular bisectors of edges `base` and `other`,
/// parametrised from the midpoint of `base`.
fn bisector_intersection(
    base_mid: Point,
    base_u: Point,
    other_mid: Point,
    other_u: Point,
) -> Option<Point> {
    let a = [[-base_u.y, other_u.y], [base_u.x, -other_u.x]];
    let rhs = [other_mid.x - base_mid.x, other_mid.y - base_mid.y];
    let (x, _) = solve_2x2(a, rhs)?;
    let s = x[0];
    Some(Point {
        x: -base_u.y * s + base_mid.x,
        y: s * base_u.x + base_mid.y,
    })
}

impl Engine {
    pub(crate) fn pnt_pnt_pnt_node(&self, i: usize, j: usize, k: usize) -> (bool, Point, f64) {
        let mut idx = [i, j, k];
        idx.sort_unstable();
        let [i, j, k] = idx;
        let pi = self.pnt_point(i);
        if i == j || j == k {
            if i == k {
                return (true, pi, 0.0);
            }
            let (center, r2) = pnt_pnt_cntr(pi, self.pnt_point(k));
            return (true, center, r2);
        }
        let arc = self.pnts[i].apx_arc;
        if arc != NIL && arc == self.pnts[j].apx_arc && arc == self.pnts[k].apx_arc {
            let center = Point {
                x: self.arcs[arc].cx,
                y: self.arcs[arc].cy,
            };
            return (true, center, dist2(center, pi));
        }
        circumcircle_3pts(pi, self.pnt_point(j), self.pnt_point(k))
    }

    pub(crate) fn pnt_point(&self, i: usize) -> Point {
        Point {
            x: self.pnts[i].x,
            y: self.pnts[i].y,
        }
    }
}

pub(crate) fn circumcircle_3pts(a: Point, b: Point, c: Point) -> (bool, Point, f64) {
    let (ab_x, ab_y) = (b.x - a.x, b.y - a.y);
    let ab = (ab_x * ab_x + ab_y * ab_y).sqrt();
    debug_assert!(ab >= 0.0, "ab >= 0.0");
    if ab <= ZERO {
        let (center, r2) = degenerate_3pts_cntr(a, b, c);
        return (false, center, r2);
    }
    let (bc_x, bc_y) = (c.x - b.x, c.y - b.y);
    let bc = (bc_x * bc_x + bc_y * bc_y).sqrt();
    debug_assert!(bc >= 0.0, "bc >= 0.0");
    if bc <= ZERO {
        let (center, r2) = degenerate_3pts_cntr(b, c, a);
        return (false, center, r2);
    }
    let (ca_x, ca_y) = (a.x - c.x, a.y - c.y);
    let ca = (ca_x * ca_x + ca_y * ca_y).sqrt();
    debug_assert!(ca >= 0.0, "ca >= 0.0");
    if ca <= ZERO {
        let (center, r2) = degenerate_3pts_cntr(c, a, b);
        return (false, center, r2);
    }

    let u_ab = Point {
        x: ab_x / ab,
        y: ab_y / ab,
    };
    let u_bc = Point {
        x: bc_x / bc,
        y: bc_y / bc,
    };
    let u_ca = Point {
        x: ca_x / ca,
        y: ca_y / ca,
    };
    let cos_b = (u_bc.x * u_ab.x + u_bc.y * u_ab.y).abs();
    let cos_c = (u_ca.x * u_bc.x + u_ca.y * u_bc.y).abs();
    let cos_a = (u_ca.x * u_ab.x + u_ab.y * u_ca.y).abs();
    let mut min_cos = cos_a;
    if cos_c <= cos_b {
        if cos_c < min_cos {
            min_cos = cos_c;
        }
    } else if cos_b < min_cos {
        min_cos = cos_b;
    }

    let m_ab = Point {
        x: (a.x + b.x) * 0.5,
        y: (a.y + b.y) * 0.5,
    };
    let m_bc = Point {
        x: (c.x + b.x) * 0.5,
        y: (c.y + b.y) * 0.5,
    };
    let m_ca = Point {
        x: (c.x + a.x) * 0.5,
        y: (c.y + a.y) * 0.5,
    };
    let solved = if cos_b == min_cos {
        if ca <= ab {
            bisector_intersection(m_ca, u_ca, m_ab, u_ab)
        } else {
            bisector_intersection(m_ab, u_ab, m_ca, u_ca)
        }
    } else if cos_c == min_cos {
        if bc <= ab {
            bisector_intersection(m_bc, u_bc, m_ab, u_ab)
        } else {
            bisector_intersection(m_ab, u_ab, m_bc, u_bc)
        }
    } else if bc <= ca {
        bisector_intersection(m_bc, u_bc, m_ca, u_ca)
    } else {
        bisector_intersection(m_ca, u_ca, m_bc, u_bc)
    };
    if let Some(center) = solved {
        return (true, center, dist2(a, center));
    }

    let sum_x = a.x + b.x;
    let sum_y = a.y + b.y;
    let mx = sum_x * 0.5;
    let my = sum_y * 0.5;
    let dx = mx - c.x;
    let dy = my - c.y;
    let den = dy * u_ab.x + dx * -u_ab.y;
    let den = den + den;
    if den > ZERO || den < -ZERO {
        let h2 = ab * ab * 0.25;
        let t = (h2 - (dy * dy + dx * dx)) / den;
        let center = Point {
            x: mx + -u_ab.y * t,
            y: my + t * u_ab.x,
        };
        return (true, center, t * t + h2);
    }
    let center = Point {
        x: (sum_x + c.x) / 3.0,
        y: (sum_y + c.y) / 3.0,
    };
    (false, center, dist2(a, center))
}
