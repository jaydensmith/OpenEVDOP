use crate::state::*;

const CLEARANCE_TOLERANCE: f64 = 1e-9;

impl Engine {
    fn site_distance(&self, site: usize, site_type: SiteType, x: f64, y: f64) -> Option<f64> {
        match site_type {
            SiteType::Pnt => {
                let p = self.pnts[site];
                Some(((p.x - x).powi(2) + (p.y - y).powi(2)).sqrt())
            }
            SiteType::Seg => {
                let s = self.segs[site];
                let (p1, p2) = (self.pnts[s.i1], self.pnts[s.i2]);
                let (dx, dy) = (p2.x - p1.x, p2.y - p1.y);
                let t = (((x - p1.x) * dx + (y - p1.y) * dy) / (dx * dx + dy * dy)).clamp(0.0, 1.0);
                Some(((p1.x + t * dx - x).powi(2) + (p1.y + t * dy - y).powi(2)).sqrt())
            }
            SiteType::Arc | SiteType::NoSite => None,
        }
    }

    /// Every finite node must be equidistant (= its clearance) to the sites left and right of each incident edge.
    pub(crate) fn vd_violations(&self) -> Vec<String> {
        let mut out = Vec::new();
        for (e, edge) in self.edges.iter().enumerate() {
            if edge.deleted {
                continue;
            }
            for n in [edge.n1, edge.n2] {
                let node = self.nodes[n];
                if node.x.abs() > self.large || node.y.abs() > self.large {
                    continue;
                }
                let r = if self.nodes_store_r2 {
                    node.r.sqrt()
                } else {
                    node.r
                };
                for (site, t) in [(edge.lft, edge.ltype), (edge.rgt, edge.rtype)] {
                    if let Some(d) = self.site_distance(site, t, node.x, node.y)
                        && (d - r).abs() > CLEARANCE_TOLERANCE * r.max(1.0)
                    {
                        out.push(format!(
                            "edge {e} node {n} ({:.6},{:.6}) r={r:.9} but dist to {t:?} {site} = {d:.9}",
                            node.x, node.y
                        ));
                    }
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_engine() -> Engine {
        let mut v = Engine::new();
        let sq = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];
        for i in 0..4 {
            let (a, b) = (sq[i], sq[(i + 1) % 4]);
            v.api_add_seg(a.0, a.1, b.0, b.1);
        }
        v
    }

    #[test]
    fn square_vd_nodes_are_equidistant_to_their_sites() {
        let mut v = square_engine();
        v.api_compute_vd(true, 0.001).expect("compute_vd");
        let live_edges = v.edges.iter().filter(|e| !e.deleted).count();
        println!(
            "scale_factor={} shift=({}, {}) nodes={} live edges={}",
            v.scale_factor,
            v.shift_x,
            v.shift_y,
            v.nodes.len(),
            live_edges
        );
        for (i, n) in v.nodes.iter().enumerate() {
            println!(
                "node {i}: ({:.6},{:.6}) r={:.6} status={:?}",
                n.x, n.y, n.r, n.status
            );
        }
        let bad = v.vd_violations();
        for b in &bad {
            println!("{b}");
        }
        assert!(bad.is_empty(), "{} clearance violations", bad.len());
    }
}
