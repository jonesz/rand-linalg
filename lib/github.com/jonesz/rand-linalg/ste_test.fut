-- | ignore

import "ste"
import "dist"
import "../../diku-dk/linalg/linalg"

module D = rademacher f32 { def seed = 06041995i64 }
module H = hutchinson f32
module L = mk_linalg f32

-- ==
-- entry: test_hutchinson_rademacher
-- input { [[1.0f32, 2.0f32], [3.0f32, 4.0f32]] }
-- output { 5.0f32 }
entry test_hutchinson_rademacher [n] (A: [n][n]f32) =
	let matvec x = L.matvecmul_row A x 
	in H.ste 100i64 D.rand matvec

-- ==
-- entry: test_hutchinson_rademacher_random
-- compiled random input { [10][10]f32 }
-- output { true }
entry test_hutchinson_rademacher_random [n] (A: [n][n]f32) =
	let tr = map (\i -> A[i][i]) (iota n) |> reduce (+) 0f32
	let matvec x = L.matvecmul_row A x 
	in (H.ste 100i64 D.rand matvec) == tr
