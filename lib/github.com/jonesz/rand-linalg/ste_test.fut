-- | ignore

import "ste"
import "dist"
import "test_matrices"
import "../../diku-dk/linalg/linalg"
import "cbrng"

def sleeve = 19950406i64
module D = rademacher_distribution f32 u32 i64 squares32
module H = hutchinson f32
module L = mk_linalg f32

def tr [n] (A: [n][n]f32) = map (\i -> A[i][i]) (iota n) |> reduce (+) 0f32

-- ==
-- entry: test_hutchinson_polyDecaySlow
-- input { [1i64, 10i64, 1000i64] [10i64, 100i64, 10000i64] 1i64 }
-- output { true }
entry test_hutchinson_polyDecaySlow R n m =
  let dist = D.construct sleeve ()
  in map2 (\R_i n_i ->
             let A = polyDecaySlow R_i n_i
             let matvec x = L.matvecmul_row A x
             in (f32.-) (H.ste m (D.rand dist) matvec) (tr A) |> f32.abs |> (f32.>=) 0.000001_f32)
          R
          n
     |> and

-- ==
-- entry: test_hutchinson_polyDecayMed
-- input { [1i64, 10i64, 1000i64] [10i64, 100i64, 10000i64] 1i64 }
-- output { true }
entry test_hutchinson_polyDecayMed R n m =
  let dist = D.construct sleeve ()
  in map2 (\R_i n_i ->
             let A = polyDecayMed R_i n_i
             let matvec x = L.matvecmul_row A x
             in (f32.-) (H.ste m (D.rand dist) matvec) (tr A) |> f32.abs |> (f32.>=) 0.000001_f32)
          R
          n
     |> and

-- ==
-- entry: test_hutchinson_polyDecayFast
-- input { [1i64, 10i64, 1000i64] [10i64, 100i64, 10000i64] 1i64 }
-- output { true }
entry test_hutchinson_polyDecayFast R n m =
  let dist = D.construct sleeve ()
  in map2 (\R_i n_i ->
             let A = polyDecayFast R_i n_i
             let matvec x = L.matvecmul_row A x
             in (f32.-) (H.ste m (D.rand dist) matvec) (tr A) |> f32.abs |> (f32.>=) 0.000001_f32)
          R
          n
     |> and

-- ==
-- entry: test_hutchinson_expDecaySlow
-- input { [1i64, 10i64, 1000i64] [10i64, 100i64, 10000i64] 1i64 }
-- output { true }
entry test_hutchinson_expDecaySlow R n m =
  let dist = D.construct sleeve ()
  in map2 (\R_i n_i ->
             let A = expDecaySlow R_i n_i
             let matvec x = L.matvecmul_row A x
             in (f32.-) (H.ste m (D.rand dist) matvec) (tr A) |> f32.abs |> (f32.>=) 0.000001_f32)
          R
          n
     |> and

-- ==
-- entry: test_hutchinson_expDecayMed
-- input { [1i64, 10i64, 1000i64] [10i64, 100i64, 10000i64] 1i64 }
-- output { true }
entry test_hutchinson_expDecayMed R n m =
  let dist = D.construct sleeve ()
  in map2 (\R_i n_i ->
             let A = expDecayMed R_i n_i
             let matvec x = L.matvecmul_row A x
             in (f32.-) (H.ste m (D.rand dist) matvec) (tr A) |> f32.abs |> (f32.>=) 0.000001_f32)
          R
          n
     |> and

-- ==
-- entry: test_hutchinson_expDecayFast
-- input { [1i64, 10i64, 1000i64] [10i64, 100i64, 10000i64] 1i64 }
-- output { true }
entry test_hutchinson_expDecayFast R n m =
  let dist = D.construct sleeve ()
  in map2 (\R_i n_i ->
             let A = expDecayFast R_i n_i
             let matvec x = L.matvecmul_row A x
             in (f32.-) (H.ste m (D.rand dist) matvec) (tr A) |> f32.abs |> (f32.>=) 0.000001_f32)
          R
          n
     |> and
