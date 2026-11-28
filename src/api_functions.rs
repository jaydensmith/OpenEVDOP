use crate::consts::*;
use crate::error::OffsetError;
use crate::real::Real;
use crate::state::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Clockwise,
    CounterClockwise,
}

/// One curve of an offset loop, in input coordinates; arcs run from `start` to `end` in `orientation`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OffsetCurve {
    Line {
        start: (f64, f64),
        end: (f64, f64),
    },
    Arc {
        start: (f64, f64),
        end: (f64, f64),
        center: (f64, f64),
        orientation: Orientation,
    },
}

impl<R: Real> Engine<R> {
    pub(crate) fn api_add_seg(&mut self, xs: f64, ys: f64, xe: f64, ye: f64) {
        let i1 = self.store_pnt(xs, ys);
        let i2 = self.store_pnt(xe, ye);
        self.store_seg(i1, i2);
    }

    /// Arcs are stored counter-clockwise from `i1` to `i2`; `ccw` records the input orientation.
    pub(crate) fn api_add_arc(
        &mut self,
        start: (f64, f64),
        end: (f64, f64),
        center: (f64, f64),
        ccw: bool,
    ) {
        let i1 = self.store_pnt(start.0, start.1);
        let i2 = self.store_pnt(end.0, end.1);
        if ccw {
            self.store_arc(i1, i2, center.0, center.1, true);
        } else {
            self.store_arc(i2, i1, center.0, center.1, false);
        }
    }

    pub(crate) fn api_compute_vd(
        &mut self,
        recover_arcs: bool,
        apx_abs: f64,
    ) -> Result<(), OffsetError> {
        loop {
            self.compute_vd(recover_arcs, apx_abs)?;
            if let Some(err) = self.failure.take() {
                return Err(err);
            }
            if !self.restart {
                return Ok(());
            }
        }
    }

    /// Offsets at distances `t, t + d, t + 2d, ...`; with neither side selected both sides are
    /// offset.
    pub(crate) fn api_compute_offset(
        &mut self,
        t: f64,
        d: f64,
        left: bool,
        right: bool,
    ) -> Result<(), OffsetError> {
        if !(ZERO < t) {
            return Ok(());
        }
        let (left, right) = if left && right {
            (false, false)
        } else {
            (left, right)
        };
        if !self.vd_computed {
            return Err(OffsetError::VdNotComputed);
        }
        self.new_offset_needed = true;
        self.compute_offset(t, d, left, right);
        Ok(())
    }

    fn unscale(&self, x: f64, y: f64) -> (f64, f64) {
        (
            (R::of(x) / self.scale_factor + self.shift_x).stored(),
            (R::of(y) / self.scale_factor + self.shift_y).stored(),
        )
    }

    /// The traced offset loops in input coordinates, every curve kept, however short, so each ends
    /// where the next starts.
    pub(crate) fn offset_curves(&self) -> Vec<Vec<OffsetCurve>> {
        self.offset_loops
            .iter()
            .map(|l| {
                let elems = &self.offset_elems[l.start..=l.end];
                (0..elems.len())
                    .map(|k| self.offset_curve(elems[k], elems[(k + 1) % elems.len()]))
                    .collect()
            })
            .collect()
    }

    fn offset_curve(&self, elem: OffsetElem, next: OffsetElem) -> OffsetCurve {
        let start = self.unscale(elem.x, elem.y);
        let end = self.unscale(next.x, next.y);
        match elem.kind {
            OffsetKind::Line => OffsetCurve::Line { start, end },
            OffsetKind::PntArc => {
                let c = self.pnts[elem.site];
                OffsetCurve::Arc {
                    start,
                    end,
                    center: self.unscale(c.x, c.y),
                    orientation: Orientation::Clockwise,
                }
            }
            OffsetKind::ApxArcInside | OffsetKind::ApxArcOutside => {
                let a = self.arcs[elem.site];
                let orientation = if elem.kind == OffsetKind::ApxArcInside {
                    Orientation::CounterClockwise
                } else {
                    Orientation::Clockwise
                };
                OffsetCurve::Arc {
                    start,
                    end,
                    center: self.unscale(a.cx, a.cy),
                    orientation,
                }
            }
        }
    }
}
