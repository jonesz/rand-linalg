-- | ignore

import "dist"
import "cbrng"

def sleeve = 19950406i64
module R = rademacher_distribution i32 u32 i64 squares32

-- ==
-- entry: rademacher_is_statistically_uniform
-- input  { 1_000_000_000i64 }
-- output { true }
entry rademacher_is_statistically_uniform n =
  let dist = R.construct sleeve ()
  -- `\sigma^2 = 1`; the sum `\sigma^2 = n`; the sum `\sigma = sqrt(n)`.
  -- 0.0063% of the seeds should be outside of `4*sqrt(n)`.
  let sigma = f32.i64 n |> f32.sqrt |> f32.ceil |> i32.f32 |> (*) 4_i32
  in map (R.rand dist) (iota n) |> reduce (+) 0_i32 |> i32.abs |> (>=) sigma
