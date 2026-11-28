# OpenEVDOP

Open source Rust engine for computing 2D Euclidean Voronoi diagrams of points, straight-line segments, and circular arcs. It computes offset curves of input line segments and circular arcs.

```rust
use open_evdop::{OffsetCurve, Orientation, Offsetter};

let mut v = Offsetter::new();
v.parameters.offset_radius = 1.0;
v.add_line((0.0, 0.0), (10.0, 0.0));
v.add_arc((10.0, 0.0), (10.0, 4.0), (10.0, 2.0), Orientation::CounterClockwise);
// ...close the shape...
v.compute_offset()?;
for l in 0..v.num_offset_loops() {
    for c in 0..v.num_offset_curves(l).unwrap() {
        match v.curve(l, c).unwrap() {
            OffsetCurve::Line { start, end } => { /* ... */ }
            OffsetCurve::Arc { start, end, center, orientation } => { /* ... */ }
        }
    }
}
```

## Gotchas

- **Default parameters.** `offset_radius = 5.0`, `apx_arcs = 0.001`, `right_offset = true`, `recover_arcs = true`, `left_offset = false`, `depletion = 0`, `tolerance = 0`.
- **Selecting both sides.** Setting both `left_offset` and `right_offset` produces offsets on both sides, which is the same as setting neither. The side is relative to each input element's direction, so for a counter-clockwise closed shape, right means outside.
- **Offset orientation.** Offset loops are oriented so that arcs around input vertices always run clockwise. Offsets of input arcs run clockwise when they lie outside the input arc and counter-clockwise when they lie inside.
- **One-sided offsets need consistently oriented input.** If the input has a dangling stub or mixed directions, a one-sided offset may return no loops for the side the stub doesn't consistently bound. A two-sided offset still works.
- **Arcs require approximation.** Input arcs are approximated by segments with a maximum deviation of `apx_arcs`. With `apx_arcs <= 0` and arcs in the input, `compute_offset` returns `Err(OffsetError::ArcsNotApproximated)`.
- **Recovery is approximate.** With `recover_arcs = true` the output follows the input arcs, but only to within `apx_arcs`. With `recover_arcs = false` you get the offset of the approximating polygon: many short lines and vertex arcs.
- **Errors report no offsets.** After an `Err` from `compute_offset` there are no loops, but your input is kept, so you can fix the parameters and call it again. `compute_offset` can also be called again after changing parameters.
- **Start-point snapping.** In `add_line` and `add_arc`, a start point is moved onto the previous end point only when `0 < gap <= tolerance`. A gap of exactly 0, or one larger than `tolerance`, is left alone. The methods return `Some(gap)` when they snapped. `reset_input` clears the previous end point, so new geometry never snaps to discarded geometry.
- **Successive offsets.** `depletion` is the offset step: a non-zero value produces loops at `radius`, `radius + depletion`, `radius + 2·depletion`, … up to the largest meaningful offset. With 0 you get a single offset.
- **Intersecting input is repaired, not rejected.** Crossing or touching segments are split and the computation restarts, up to 500 times before `Err(OffsetError::TooManyRestarts)`. If no Voronoi node can be computed even by the fallback routines, the result is `Err(OffsetError::NumericalFailure)`.
- **No offset is computed for non-positive radii.** For `offset_radius <= 1e-11`, `compute_offset` returns `Ok` with zero loops.
- **Every curve is returned.** Each curve ends bit for bit where the next starts. At a tangent joint, such as where two arcs form a circle, the loop holds a vertex arc of no length: its two ends are the same point, or lie within rounding of each other in either order. Treat such an arc as having no sweep, whatever its orientation says.
- **Arithmetic.** `Arithmetic::Extended` computes like an x87 FPU at 64-bit precision, rounding to `f64` only where values are stored; it is much slower than `Double`. Both modes use correctly rounded `sin`, `cos` and `acos`, so results are identical on every platform; each call takes tens of microseconds, so input with many arcs takes longer to approximate than with the platform's own functions.

## Tests

- `tests/offset_geometry.rs` and `tests/stress.rs` check the geometry in both arithmetic modes. Every curve must end bit for bit where the next starts, and every loop must lie at exactly the offset distance from the input. The inputs include random polygons of up to 5,000 vertices.
- `tests/reference_outputs.rs` replays `tests/reference/*.script` and compares the result with the recorded output: `*.expected` for `Arithmetic::Double` and `*.extended.expected` for `Arithmetic::Extended`. Values must agree exactly.

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.
