-- | ignore

import "dist"
import "cbrng"

module R = rademacher_distribution i32 u32 i64 squares32
module N = normal_distribution f32 u32 i64 squares32

-- ==
-- entry: rademacher_is_statistically_uniform
-- compiled random input { [100]i64 1_000_000i64 }
-- output { true }
entry rademacher_is_statistically_uniform seeds n =
  map (\seed ->
         let dist = R.construct seed ()
         -- `\sigma^2 = 1`; the sum `\sigma^2 = n`; the sum `\sigma = sqrt(n)`.
         -- 0.0063% of the seeds should be outside of `4*sqrt(n)`.
         let sigma = f32.i64 n |> f32.sqrt |> f32.ceil |> i32.f32 |> (*) 4_i32
         in map (R.rand dist) (iota n) |> reduce (+) 0_i32 |> i32.abs |> (>=) sigma)
      seeds
  |> and

-- ==
-- entry: normal_is_statistically_correct
-- compiled random input { [100]i64 1_000i64 }
-- output { true }
entry normal_is_statistically_correct seeds n =
  map (\seed ->
         let dist = N.construct seed {mean = 0f32, stddev = 1f32}
         -- `\sigma^2 = 1`; the sum `\sigma^2 = n`; the sum `\sigma = sqrt(n)`.
         -- 0.0063% of the seeds should be outside of `4*sqrt(n)`.
         let sigma = f32.i64 n |> f32.sqrt |> f32.ceil |> (*) 4_f32
         in map (N.rand dist) (iota n) |> reduce (+) 0_f32 |> f32.abs |> (>=) sigma)
      seeds
  |> and
