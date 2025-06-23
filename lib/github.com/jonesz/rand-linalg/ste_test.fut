-- | ignore

import "ste"
import "dist"
import "../../diku-dk/linalg/linalg"
import "../../diku-dk/cpprandom/random"

module D = rademacher_distribution f32 minstd_rand {def seed = 06041995i32}
module H = hutchinson f32
module L = mk_linalg f32

-- ==
-- entry: test_hutchinson_rademacher
-- input { [[2.0f32, 2.0f32, 2.0f32], [5f32, 0.5f32, 33f32], [2.5f32, 3.0f32, 4.0f32]] }
-- output { 6.5f32 }
entry test_hutchinson_rademacher [n] (A: [n][n]f32) =
  let matvec x = L.matvecmul_row A x
  in H.ste 100i32 D.rand matvec

-- ==
-- entry: test_hutchinson_rademacher_random
-- compiled random input { [1000][1000]f32 }
-- output { 0.001f32 }
entry test_hutchinson_rademacher_random [n] (A: [n][n]f32) =
  let tr = map (\i -> A[i][i]) (iota n) |> reduce (+) 0f32
  let matvec x = L.matvecmul_row A x
  in (f32.-) tr (H.ste 100i32 D.rand matvec) |> f32.abs
