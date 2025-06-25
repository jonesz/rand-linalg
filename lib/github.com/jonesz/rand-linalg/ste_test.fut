-- | ignore

import "ste"
import "dist"
import "../../diku-dk/linalg/linalg"
import "../../diku-dk/cpprandom/random"

module D = rademacher_distribution f32 u32 minstd_rand {def seed = 06041995i32}
module H = hutchinson f32
module L = mk_linalg f32

-- ==
-- entry: test_hutchinson_rademacher_eye
-- input { }
-- output { 0.01f32 }
entry test_hutchinson_rademacher_eye =
  let n = 1000
  let A = L.eye n
  let tr = map (\i -> A[i][i]) (iota n) |> reduce (+) 0f32
  let matvec x = L.matvecmul_row A x
  -- https://www.ethanepperly.com/index.php/2024/01/04/how-good-can-stochastic-trace-estimates-be/
  -- NOTE: This is for symmetric/positive definit matrices -- our test is this, but it's I...
  in (f32.-) (H.ste 10_000i32 D.rand matvec) tr |> f32.abs |> flip (f32./) tr

-- ==
-- entry: test_hutchinson_rademacher_diagonally_dominant
-- compiled random input { [1000]f32 }
-- output { 0.01f32 }
entry test_hutchinson_rademacher_diagonally_dominant [n] (A: [n]f32) =
  let A = L.todiag A
  let tr = map (\i -> A[i][i]) (iota n) |> reduce (+) 0f32
  let matvec x = L.matvecmul_row A x
  -- https://www.ethanepperly.com/index.php/2024/01/04/how-good-can-stochastic-trace-estimates-be/
  -- NOTE: This is for symmetric/positive definit matrices -- our test is diagonally dominant, but
  -- there's nothing off the diagonal.
  in (f32.-) (H.ste 10_000i32 D.rand matvec) tr |> f32.abs |> flip (f32./) tr

-- ==
-- entry: test_hutchinson_rademacher_random
-- compiled random input { [1000][1000]f32 }
-- output { 0.01f32 }
entry test_hutchinson_rademacher_random [n] (A: [n][n]f32) =
  let tr = map (\i -> A[i][i]) (iota n) |> reduce (+) 0f32
  let matvec x = L.matvecmul_row A x
  -- https://www.ethanepperly.com/index.php/2024/01/04/how-good-can-stochastic-trace-estimates-be/
  -- NOTE: This is for symmetric/positive definit matrices, whereas our test isn't.
  in (f32.-) (H.ste 10_000i32 D.rand matvec) tr |> f32.abs |> flip (f32./) tr
