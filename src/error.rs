use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffsetError {
    ArcsNotApproximated,
    TooManyRestarts,
    VdNotComputed,
    NumericalFailure,
}

impl fmt::Display for OffsetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OffsetError::ArcsNotApproximated => {
                write!(
                    f,
                    "input contains arcs but apx_arcs <= 0, so they cannot be approximated"
                )
            }
            OffsetError::TooManyRestarts => {
                write!(
                    f,
                    "Voronoi computation restarted too often due to numerical trouble"
                )
            }
            OffsetError::VdNotComputed => write!(f, "Voronoi diagram has not been computed"),
            OffsetError::NumericalFailure => {
                write!(
                    f,
                    "no Voronoi node could be computed, even by the fallback routines"
                )
            }
        }
    }
}

impl std::error::Error for OffsetError {}
