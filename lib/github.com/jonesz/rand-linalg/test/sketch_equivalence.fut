-- | ignore

import "../sketch"
import "../test_matrices"
import "../../cbrng-fut/distribution"
import "../../cbrng-fut/cbrng"
import "../../../diku-dk/linalg/linalg"

module SK = mk_sketch f32 (rademacher_distribution f32 u32 squares32)
module TM = tm f32
module L = mk_linalg f32

-- ==
-- entry: left_sketch_equivalence
-- compiled random input { i64 }
-- output { true }
entry left_sketch_equivalence seed =
  let seed = squares32.construct seed
  let A = TM.lowRankHiNoise (seed + 1) 50 100
  let oracle = L.matvecmul_row (transpose A)
  let sk_1 = SK.left.dense.sketch seed () 20 A |> flatten
  let sk_2 = SK.left.oracle.sketch seed () 20 oracle |> flatten
  in map2 (==) sk_1 sk_2 |> and

-- ==
-- entry: right_sketch_equivalence
-- compiled random input { i64 }
-- output { true }
entry right_sketch_equivalence seed =
  let seed = squares32.construct seed
  let A = TM.lowRankHiNoise (seed + 1) 50 100
  let oracle = L.matvecmul_row A
  let sk_1 = SK.right.dense.sketch seed () 20 A |> flatten
  let sk_2 = SK.right.oracle.sketch seed () 20 oracle |> flatten
  in map2 (f32.==) sk_1 sk_2 |> and
