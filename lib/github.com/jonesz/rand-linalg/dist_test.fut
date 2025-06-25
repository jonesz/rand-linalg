-- | ignore

import "dist"
import "../../diku-dk/cpprandom/random"

module R = rademacher_distribution f32 u32 minstd_rand {def seed = 7i32}

-- ==
-- entry: rademacher_is_uniform
-- input  { 1_000_000i64 }
-- output { 500_000i64 }
entry rademacher_is_uniform n =
  map (\i -> i32.i64 i |> R.rand ()) (iota n) |> filter ((==) 1f32) |> length
