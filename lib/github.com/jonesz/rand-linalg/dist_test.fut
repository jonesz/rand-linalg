-- | ignore

import "dist"
import "cbrng"

module D = rademacher_distribution f32 u32 (squares32 {def key = 194956i64})

-- ==
-- entry: rademacher_is_uniform
-- input  { 1_000_000i64 }
-- output { 500_000i64 }
entry rademacher_is_uniform n =
  map (D.rand ()) (iota n) |> filter ((==) 1f32) |> length
