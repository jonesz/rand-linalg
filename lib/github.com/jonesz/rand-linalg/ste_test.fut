-- | ignore

import "../../diku-dk/linalg/linalg"
import "../cbrng-fut/cbrng"
import "../cbrng-fut/distribution"
import "ste"
import "test_matrices"

module mk_chebyshev_rademacher_test (R: real) (S: ste with t = R.t) = {
  module L = mk_linalg R
  module D = rademacher_distribution R u32 i64 squares32

  def ste A s seed = S.ste s (D.rand seed) (L.matvecmul_row A)
  def tr [n] (A: [n][n]R.t) = map (\i -> A[i][i]) (iota n) |> reduce (R.+) (R.i64 0)

  -- Variance for the Rademacher test vector.
  def var [n] (A: [n][n]R.t) =
    let f_sq = map (\i -> (R.**) i (R.i64 2)) (flatten A) |> reduce (R.+) (R.i64 0)
    let d_sq = map (\i -> (R.**) A[i][i] (R.i64 2)) (iota n) |> reduce (R.+) (R.i64 0)
    in (R.-) f_sq d_sq |> (R.*) (R.i64 2)

  -- P{|X_s - tr(A)| >= t * tr(A)} <= Var[X] / s(tr A)^2 t^2.
  def chebyshev_bound A s t =
    let top = var A 
    let bot = (R.**) (tr A) (R.i64 2) |> (R.*) ((R.**) t (R.i64 2)) |> (R.*) (R.i64 s)
    in (R./) top bot 

  def test [n] (A: [n][n]R.t) s (t: f32) seed =
    let t = R.f32 t
    let m = 100i64
    let bound = chebyshev_bound A s t

    let f seed =
      let seed = squares32.construct seed
      let true_tr = tr A
      let esti_tr = ste A s seed
      let left  = (R.-) esti_tr true_tr |> R.abs
      let right = (R.*) true_tr t
      in (R.>=) left right |> i64.bool
 
    in map (i64.+ seed) (iota m) |> map f
      |> reduce (+) 0i64 |> R.i64 |> flip (R./) (R.i64 m) |> flip (R.<=) bound
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
entry test_hutchinson_polyDecaySlow_chebyshev seed R n s =
  H.test (T.polyDecaySlow R n) s 1.0_f32 seed

-- entry: test_hutchinson_polyDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayMed_chebyshev seed R n s =
  H.test (T.polyDecayMed R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_polyDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayFast_chebyshev seed R n s =
  H.test (T.polyDecayFast R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_expDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecaySlow_chebyshev seed R n s =
  H.test (T.expDecaySlow R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_expDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayMed_chebyshev seed R n s =
  H.test (T.expDecayMed R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_expDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayFast_chebyshev seed R n s =
  H.test (T.expDecayFast R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_lowRankLowNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankLowNoise_chebyshev seed R n s =
  H.test (T.lowRankLowNoise (seed + 1) R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_lowRankMedNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankMedNoise_chebyshev seed R n s =
  H.test (T.lowRankMedNoise (seed + 1) R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchinson_lowRankHiNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankHiNoise_chebyshev seed R n s =
  H.test (T.lowRankHiNoise (seed + 1) R n) s 1.0_f32 seed

-- ~
-- ~ HUTCHPLUSPLUS TESTS
-- ~

-- entry: test_hutchplusplus_polyDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_polyDecaySlow_chebyshev seed R n s =
  HPP.test (T.polyDecaySlow R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_polyDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_polyDecayMed_chebyshev seed R n s =
  HPP.test (T.polyDecayMed R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_polyDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_polyDecayFast_chebyshev seed R n s =
  HPP.test (T.polyDecayFast R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_expDecaySlow_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_expDecaySlow_chebyshev seed R n s =
  HPP.test (T.expDecaySlow R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_expDecayMed_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_expDecayMed_chebyshev seed R n s =
  HPP.test (T.expDecayMed R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_expDecayFast_chebyshev
-- compiled random input { i64 10i64 100i64 3i64 }
-- compiled random input { i64 90i64 100i64 3i64 }
-- output { true }
entry test_hutchplusplus_expDecayFast_chebyshev seed R n s =
  HPP.test (T.expDecayFast R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_lowRankLowNoise_chebyshev
-- compiled random input { i64 10i64 100i64 30i64 }
-- compiled random input { i64 90i64 100i64 300i64 }
-- output { true }
entry test_hutchplusplus_lowRankLowNoise_chebyshev seed R n s =
  HPP.test (T.lowRankLowNoise (seed * 0xA) R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_lowRankMedNoise_chebyshev
-- compiled random input { i64 10i64 100i64 30i64 }
-- compiled random input { i64 10i64 100i64 300i64 }
-- output { true }
entry test_hutchplusplus_lowRankMedNoise_chebyshev seed R n s =
  HPP.test (T.lowRankMedNoise (seed * 0xA) R n) s 1.0_f32 seed

-- ==
-- entry: test_hutchplusplus_lowRankHiNoise_chebyshev
-- compiled random input { i64 10i64 100i64 30i64 }
-- compiled random input { i64 10i64 100i64 300i64 }
-- output { true }
entry test_hutchplusplus_lowRankHiNoise_chebyshev seed R n s =
  HPP.test (T.lowRankHiNoise (seed * 0xA) R n) s 1.0_f32 seed
