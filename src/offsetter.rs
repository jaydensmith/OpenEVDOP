use crate::api_functions::{OffsetCurve, Orientation};
use crate::error::OffsetError;
use crate::real::Real;
use crate::state::Engine;
use x87_extended::Extended;

const DEFAULT_OFFSET_RADIUS: f64 = 5.0;
const DEFAULT_APX_ARCS: f64 = 0.001;

/// The arithmetic an offset is computed in.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Arithmetic {
    /// `f64` throughout.
    #[default]
    Double,
    /// The x87's 80-bit registers at 64-bit precision; values are rounded to `f64` where they are
    /// stored.
    Extended,
}

/// Offset computation settings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OffsetParameters {
    pub left_offset: bool,
    pub right_offset: bool,
    pub offset_radius: f64,
    pub recover_arcs: bool,
    /// Maximum absolute deviation when approximating input arcs; lower values cost more time.
    pub apx_arcs: f64,
    /// Step between successive offsets (`offset_radius`, `+ depletion`, ...); 0 yields one offset.
    pub depletion: f64,
    /// Start points within this distance of the previous end point are snapped onto it.
    pub tolerance: f64,
    pub arithmetic: Arithmetic,
}

impl Default for OffsetParameters {
    fn default() -> Self {
        Self {
            left_offset: false,
            right_offset: true,
            offset_radius: DEFAULT_OFFSET_RADIUS,
            recover_arcs: true,
            apx_arcs: DEFAULT_APX_ARCS,
            depletion: 0.0,
            tolerance: 0.0,
            arithmetic: Arithmetic::Double,
        }
    }
}

/// Offset-curve computation over input line segments and circular arcs.
#[derive(Debug, Clone, Copy)]
enum InputElement {
    Line {
        start: (f64, f64),
        end: (f64, f64),
    },
    Arc {
        start: (f64, f64),
        end: (f64, f64),
        center: (f64, f64),
        ccw: bool,
    },
}

#[derive(Debug, Clone)]
pub struct Offsetter {
    pub parameters: OffsetParameters,
    last_end: Option<(f64, f64)>,
    input: Vec<InputElement>,
    loops: Vec<Vec<OffsetCurve>>,
}

impl Default for Offsetter {
    fn default() -> Self {
        Self::new()
    }
}

impl Offsetter {
    pub fn new() -> Self {
        Self {
            parameters: OffsetParameters::default(),
            last_end: None,
            input: Vec::new(),
            loops: Vec::new(),
        }
    }

    /// Clears the input geometry, the computed offsets and the end point used for snapping.
    pub fn reset_input(&mut self) {
        self.input.clear();
        self.loops.clear();
        self.last_end = None;
    }

    /// Adds a segment; returns the gap if its start was snapped onto the previous end point.
    pub fn add_line(&mut self, start: (f64, f64), end: (f64, f64)) -> Option<f64> {
        let (start, snapped) = self.snap_start(start);
        self.input.push(InputElement::Line { start, end });
        self.last_end = Some(end);
        snapped
    }

    /// Adds a circular arc from `start` to `end` around `center`;
    /// returns the gap if its start was snapped onto the previous end point.
    pub fn add_arc(
        &mut self,
        start: (f64, f64),
        end: (f64, f64),
        center: (f64, f64),
        orientation: Orientation,
    ) -> Option<f64> {
        let (start, snapped) = self.snap_start(start);
        let ccw = orientation == Orientation::CounterClockwise;
        self.input.push(InputElement::Arc {
            start,
            end,
            center,
            ccw,
        });
        self.last_end = Some(end);
        snapped
    }

    /// Builds the Voronoi diagram of the input and traces the offset loops; may be called again
    /// after changing parameters. On failure no offsets are reported and the input is kept.
    pub fn compute_offset(&mut self) -> Result<(), OffsetError> {
        let result = match self.parameters.arithmetic {
            Arithmetic::Double => self.computed::<f64>(),
            Arithmetic::Extended => self.computed::<Extended>(),
        };
        match result {
            Ok(loops) => {
                self.loops = loops;
                Ok(())
            }
            Err(error) => {
                self.loops = Vec::new();
                Err(error)
            }
        }
    }

    /// The offset loops of the input in arithmetic `R`.
    fn computed<R: Real>(&self) -> Result<Vec<Vec<OffsetCurve>>, OffsetError> {
        let p = self.parameters;
        // The engine scales, cleans and approximates its input in place, so every run starts fresh.
        let mut engine = Engine::<R>::new();
        for element in &self.input {
            match *element {
                InputElement::Line { start, end } => {
                    engine.api_add_seg(start.0, start.1, end.0, end.1)
                }
                InputElement::Arc {
                    start,
                    end,
                    center,
                    ccw,
                } => engine.api_add_arc(start, end, center, ccw),
            }
        }
        engine
            .api_compute_vd(p.recover_arcs, p.apx_arcs)
            .and_then(|()| {
                engine.api_compute_offset(
                    p.offset_radius,
                    p.depletion,
                    p.left_offset,
                    p.right_offset,
                )
            })
            .map(|()| engine.offset_curves())
    }

    pub fn num_offset_loops(&self) -> usize {
        self.loops.len()
    }

    pub fn num_offset_curves(&self, loop_index: usize) -> Option<usize> {
        self.loops.get(loop_index).map(Vec::len)
    }

    pub fn curve(&self, loop_index: usize, curve_index: usize) -> Option<OffsetCurve> {
        self.loops.get(loop_index)?.get(curve_index).copied()
    }

    fn snap_start(&self, start: (f64, f64)) -> ((f64, f64), Option<f64>) {
        let Some(prev) = self.last_end else {
            return (start, None);
        };
        let gap = match self.parameters.arithmetic {
            Arithmetic::Double => gap::<f64>(start, prev),
            Arithmetic::Extended => gap::<Extended>(start, prev),
        };
        if gap <= self.parameters.tolerance && gap != 0.0 {
            (prev, Some(gap))
        } else {
            (start, None)
        }
    }
}

/// The distance from `prev` to `start` in arithmetic `R`.
fn gap<R: Real>(start: (f64, f64), prev: (f64, f64)) -> f64 {
    let (dx, dy) = (
        (R::of(start.0) - prev.0).stored(),
        (R::of(start.1) - prev.1).stored(),
    );
    (R::of(dy) * dy + R::of(dx) * dx).sqrt().stored()
}
