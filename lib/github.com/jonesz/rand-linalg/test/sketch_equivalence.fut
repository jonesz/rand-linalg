-- | ignore

import "../sketch"
import "../../cbrng-fut/cbrng"
import "../../../diku-dk/linalg/linalg"

module SK = mk_gaussian_embedding f32 u32 squares32
module L = mk_linalg f32

-- ==
-- entry: left_sketch_equivalence_ld_lo_sketch
-- compiled random input { i64 [100][10]f32 }
-- output { true }
entry left_sketch_equivalence_ld_lo_sketch sketch_seed A =
  let sketch_seed = SK.dist.engine.construct sketch_seed

  let sk_1 = SK.sketch_ld sketch_seed 20 A |> flatten
  let sk_2 = SK.sketch_lo sketch_seed 20 (L.matvecmul_row (transpose A)) |> flatten
  let sk_3 =
    let O = SK.sketch sketch_seed 20 100
    in L.matmul O A |> flatten

  in map3 (\a b c -> a == b && b == c) sk_1 sk_2 sk_3 |> and

-- ==
-- entry: right_sketch_equivalence_ld_lo_sketch
-- compiled random input { i64 [10][100]f32 }
-- output { true }
entry right_sketch_equivalence_ld_lo_sketch sketch_seed A =
  let sketch_seed = SK.dist.engine.construct sketch_seed

  let sk_1 = SK.sketch_rd sketch_seed 20 A |> flatten
  let sk_2 = SK.sketch_ro sketch_seed 20 (L.matvecmul_row A) |> flatten
  let sk_3 =
    let O = SK.sketch sketch_seed 20 100
    in L.matmul A (transpose O) |> flatten

  in map3 (\a b c -> a == b && b == c) sk_1 sk_2 sk_3 |> and
