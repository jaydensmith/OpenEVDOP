//! The arithmetic an offset is computed in: `f64`, or the x87's 80-bit registers.

use std::ops::{Add, Div, Mul, Neg, Sub};
use x87_extended::{Extended, double};

/// A number an offset is computed in. Values are kept as `f64` wherever they are stored; an
/// expression is evaluated in `Self` and stored with `stored`.
pub(crate) trait Real:
    Copy
    + PartialEq
    + PartialOrd
    + PartialEq<f64>
    + PartialOrd<f64>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    fn of(value: f64) -> Self;
    /// The value as stored to a double.
    fn stored(self) -> f64;
    fn abs(self) -> Self;
    fn sqrt(self) -> Self;
    fn atan2(self, x: Self) -> Self;
    /// Truncated toward zero to an `i32`.
    fn truncate(self) -> i32;
}

impl Real for f64 {
    fn of(value: f64) -> f64 {
        value
    }

    fn stored(self) -> f64 {
        self
    }

    fn abs(self) -> f64 {
        f64::abs(self)
    }

    fn sqrt(self) -> f64 {
        f64::sqrt(self)
    }

    fn atan2(self, x: f64) -> f64 {
        double::atan2(self, x)
    }

    fn truncate(self) -> i32 {
        self as i32
    }
}

impl Real for Extended {
    fn of(value: f64) -> Extended {
        Extended::of(value)
    }

    fn stored(self) -> f64 {
        self.to_double()
    }

    fn abs(self) -> Extended {
        Extended::abs(self)
    }

    fn sqrt(self) -> Extended {
        Extended::sqrt(self)
    }

    fn atan2(self, x: Extended) -> Extended {
        self.fpatan(x)
    }

    /// The x87's integer indefinite, `i32::MIN`, for a NaN, an infinity or a value outside `i32`.
    fn truncate(self) -> i32 {
        Extended::truncate(self)
            .and_then(|whole| i32::try_from(whole).ok())
            .unwrap_or(i32::MIN)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::FAR_NODE_DIST2;
    use x87_extended::Extended;

    /// A temporary this many times `FAR_NODE_DIST2` overflows a double but not a register.
    const OVERFLOWING_FACTOR: f64 = 20.0;
    const BACK_IN_RANGE_DIVISOR: f64 = 40.0;

    fn kept_by<R: Real>() -> f64 {
        (R::of(1.0) + 2f64.powi(-60) - 1.0).stored()
    }

    fn overflow_in<R: Real>() -> f64 {
        (R::of(FAR_NODE_DIST2) * OVERFLOWING_FACTOR / BACK_IN_RANGE_DIVISOR).stored()
    }

    fn truncated_by<R: Real>() -> i32 {
        (R::of(3.0) - 2f64.powi(-62)).truncate()
    }

    #[test]
    fn extended_keeps_what_f64_loses() {
        let (double, extended) = (kept_by::<f64>(), kept_by::<Extended>());
        println!("(1 + 2^-60) - 1: f64 {double:e}, Extended {extended:e}");
        assert_eq!(double, 0.0);
        assert_eq!(extended, 2f64.powi(-60));
    }

    #[test]
    fn a_temporary_overflows_only_when_stored() {
        let (double, extended) = (overflow_in::<f64>(), overflow_in::<Extended>());
        println!("FAR_NODE_DIST2 × 20 / 40: f64 {double:e}, Extended {extended:e}");
        assert!(double.is_infinite());
        assert_eq!(extended, FAR_NODE_DIST2 / 2.0);
    }

    #[test]
    fn truncate_reads_the_register() {
        let (double, extended) = (truncated_by::<f64>(), truncated_by::<Extended>());
        println!("3 - 2^-62 truncated: f64 {double}, Extended {extended}");
        assert_eq!(double, 3);
        assert_eq!(extended, 2);
        assert_eq!(Extended::of(f64::NAN).truncate(), None);
        assert_eq!(
            <Extended as Real>::truncate(Extended::of(f64::NAN)),
            i32::MIN
        );
        assert_eq!(<Extended as Real>::truncate(Extended::of(3e9)), i32::MIN);
    }
}
