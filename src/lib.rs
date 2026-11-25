// Comparisons are negated deliberately: their behaviour on NaN is part of the algorithm.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

mod api_functions;
mod approx;
mod clean_data;
mod compute;
mod consts;
mod data;
mod data_off;
mod desperate;
mod error;
mod grid;
mod heap;
mod misc;
mod numerics;
mod offset;
mod offsetter;
mod pnt_pnt_pnt;
mod prepro;
mod seg_pnt_pnt;
mod seg_seg_pnt;
mod seg_seg_seg;
mod state;
mod tree;
mod vd_basic;
#[cfg(test)]
mod vd_check;
mod vd_data;
mod vd_pnt;
mod vd_seg;

pub use api_functions::{OffsetCurve, Orientation};
pub use error::OffsetError;
pub use offsetter::{Arithmetic, OffsetParameters, Offsetter};
