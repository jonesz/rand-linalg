-- | ignore

import "ste"
import "dist"
import "test_matrices"
import "../../diku-dk/linalg/linalg"
import "cbrng"

module mk_chebyshev_rademacher_test (R: real) (S: ste with t = R.t) = {
  module L = mk_linalg R
  module D = rademacher_distribution R u32 i64 squares32

  def ste A m seed =
    let matvec = L.matvecmul_row A
    in S.ste m (D.rand (D.construct seed ())) matvec

  def tr [n] (A: [n][n]R.t) =
    map (\i -> A[i][i]) (iota n) |> reduce (R.+) (R.i64 0)

  -- | Variance for the Rademacher test vector.
  def var [n] (A: [n][n]R.t) (m: i64) =
    let f_sq = map (\i -> (R.**) i (R.i64 2)) (flatten A) |> reduce (R.+) (R.i64 0)
    let d_sq = map (\i -> (R.**) A[i][i] (R.i64 2)) (iota n) |> reduce (R.+) (R.i64 0)
    in (R.-) f_sq d_sq |> flip (R./) (R.i64 m)

  -- | P{|X_s - tr(A)| >= t * tr(A)} <= Var[X] / s(tr A)^2 t^2.
  def chebyshev_bound A t (s: i64) =
    let bot = tr A |> flip (R.**) (R.i64 2) |> (R.*) ((R.**) t (R.i64 2)) |> (R.*) (R.i64 s)
    in var A s |> flip (R./) bot

  -- TODO: Move the `ste` code here; then place the check for the bound
  -- in the test function.	tiwd
  def test [n] (A: [n][n]R.t) m seed =
    let est_tr = ste A m seed
    let act_tr = tr A
    let bound = chebyshev_bound A (R.f32 0.0001_f32) m
    in (R.-) est_tr act_tr |> R.abs |> (R.>=) bound
}

module H = mk_chebyshev_rademacher_test f64 (hutchinson f64)
module HPP = mk_chebyshev_rademacher_test f64 (hutchplusplus f64)

module T = tm f64
module L = mk_linalg f32

-- ==
-- entry: test_matmul_w_matvec
-- random input { [100][10][10]f32 [100][10][3]f32 }
-- output { true }
entry test_matmul_w_matvec A B =
  map2 (\A_i B_i ->
          let j = L.matmul A_i B_i |> flatten
          let k = matmul_w_matvec (L.matvecmul_row A_i) B_i |> flatten
          in map2 (==) j k
             |> and)
       A
       B
  |> and

-- ~~
-- ~~ HUTCHINSON TESTS
-- ~~

-- entry: test_hutchinson_polyDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecaySlow_chebyshev seed R n m =
  H.test (T.polyDecaySlow R n) m seed

-- ==
-- entry: test_hutchinson_polyDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayMed_chebyshev seed R n m =
  H.test (T.polyDecayMed R n) m seed

-- ==
-- entry: test_hutchinson_polyDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayFast_chebyshev seed R n m =
  H.test (T.polyDecayFast R n) m seed

-- ==
-- entry: test_hutchinson_expDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecaySlow_chebyshev seed R n m =
  H.test (T.expDecaySlow R n) m seed

-- ==
-- entry: test_hutchinson_expDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayMed_chebyshev seed R n m =
  H.test (T.expDecayMed R n) m seed

-- ==
-- entry: test_hutchinson_expDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayFast_chebyshev seed R n m =
  H.test (T.expDecayFast R n) m seed

-- ==
-- entry: test_hutchinson_lowRankLowNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankLowNoise_chebyshev seed R n m =
  H.test (T.lowRankLowNoise (seed + 1) R n) m seed

-- ==
-- entry: test_hutchinson_lowRankMedNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankMedNoise_chebyshev seed R n m =
  H.test (T.lowRankMedNoise (seed + 1) R n) m seed

-- ==
-- entry: test_hutchinson_lowRankHiNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankHiNoise_chebyshev seed R n m =
  H.test (T.lowRankHiNoise (seed + 1) R n) m seed

-- ~
-- ~ HUTCHPLUSPLUS TESTS
-- ~

-- entry: test_hutchplusplus_polyDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_polyDecaySlow_chebyshev seed R n m =
  HPP.test (T.polyDecaySlow R n) m seed

-- ==
-- entry: test_hutchplusplus_polyDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_polyDecayMed_chebyshev seed R n m =
  HPP.test (T.polyDecayMed R n) m seed

-- ==
-- entry: test_hutchplusplus_polyDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_polyDecayFast_chebyshev seed R n m =
  HPP.test (T.polyDecayFast R n) m seed

-- ==
-- entry: test_hutchplusplus_expDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_expDecaySlow_chebyshev seed R n m =
  HPP.test (T.expDecaySlow R n) m seed

-- ==
-- entry: test_hutchplusplus_expDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_expDecayMed_chebyshev seed R n m =
  HPP.test (T.expDecayMed R n) m seed

-- ==
-- entry: test_hutchplusplus_expDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_expDecayFast_chebyshev seed R n m =
  HPP.test (T.expDecayFast R n) m seed

-- ==
-- entry: test_hutchplusplus_lowRankLowNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchplusplus_lowRankLowNoise_chebyshev seed R n m =
  HPP.test (T.lowRankLowNoise (seed + 1) R n) m seed

-- ==
-- entry: test_hutchplusplus_lowRankMedNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchplusplus_lowRankMedNoise_chebyshev seed R n m =
  HPP.test (T.lowRankMedNoise (seed + 1) R n) m seed

-- ==
-- entry: test_hutchplusplus_lowRankHiNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchplusplus_lowRankHiNoise_chebyshev seed R n m =
  HPP.test (T.lowRankHiNoise (seed + 1) R n) m seed
