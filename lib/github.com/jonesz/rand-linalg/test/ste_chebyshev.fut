-- | ignore
-- Tests against the Chebyshev bound. Based on the variance of the input
-- matrix, the estimator should be within some distance of the true trace
-- some percentage of the time.

import "../ste"
import "../test_matrices"
import "../../cbrng-fut/distribution"
import "../../cbrng-fut/cbrng"
import "../../../diku-dk/linalg/linalg"

module mk_chebyshev_rademacher (R: real) (S: ste with t = R.t) = {
  module L = mk_linalg R
  module D = rademacher_distribution R u32 squares32

  def ste A s seed = S.ste s (D.rand seed ()) (L.matvecmul_row A)
  def tr A = L.fromdiag A |> reduce (R.+) (R.i64 0)

  -- Variance for the Rademacher test vector.
  def var [n] (A: [n][n]R.t) =
    let f_sq = map (\i -> (R.**) i (R.i64 2)) (flatten A) |> reduce (R.+) (R.i64 0)
    let d_sq = tabulate n (\i -> (R.**) A[i][i] (R.i64 2)) |> reduce (R.+) (R.i64 0)
    in (R.-) f_sq d_sq |> (R.*) (R.i64 2)

  -- P{|X_s - tr(A)| >= t * tr(A)} <= Var[X] / s(tr A)^2 t^2.
  def chebyshev_bound A s t =
    let top = var A
    let bot = (R.**) (tr A) (R.i64 2) |> (R.*) ((R.**) t (R.i64 2)) |> (R.*) (R.i64 s)
    in (R./) top bot

  -- Assert that the Chebyshev bound holds for the STE.
  def test [n] (A: [n][n]R.t) s (t: f32) seed =
    let t = R.f32 t
    let m = 100i64
    let bound = chebyshev_bound A s t

    let f seed =
      let seed = squares32.construct seed
      let true_tr = tr A
      let esti_tr = ste A s seed
      let left = (R.-) esti_tr true_tr |> R.abs
      let right = (R.*) true_tr t
      in (R.>=) left right |> i64.bool

    in tabulate m (i64.+ seed) |> map f
       |> reduce (+) 0i64
       |> R.i64
       |> flip (R./) (R.i64 m)
       |> flip (R.<=) bound
}

module HC = mk_chebyshev_rademacher f32 (hutchinson f32)
module T = tm f32

-- The trace estimator should be within 10% of the actual trace.
def t = 0.10_f32

-- ==
-- entry: test_hutchinson_lowRankLowNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankLowNoise_chebyshev seed R n s =
  HC.test (T.lowRankLowNoise (seed + 1) R n) s t seed

-- ==
-- entry: test_hutchinson_lowRankMedNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankMedNoise_chebyshev seed R n s =
  HC.test (T.lowRankMedNoise (seed + 1) R n) s t seed

-- ==
-- entry: test_hutchinson_lowRankHiNoise_chebyshev
-- compiled random input { i64 10i64 100i64 10i64 }
-- compiled random input { i64 90i64 100i64 100i64 }
-- output { true }
entry test_hutchinson_lowRankHiNoise_chebyshev seed R n s =
  HC.test (T.lowRankHiNoise (seed + 1) R n) s t seed

module HPPC = mk_chebyshev_rademacher f32 (hutchplusplus f32)

-- ==
-- entry: test_hutchplusplus_lowRankLowNoise_chebyshev
-- compiled random input { i64 10i64 100i64 9i64 }
-- compiled random input { i64 90i64 100i64 99i64 }
-- output { true }
entry test_hutchplusplus_lowRankLowNoise_chebyshev seed R n s =
  HPPC.test (T.lowRankLowNoise (seed + 1) R n) s t seed

-- ==
-- entry: test_hutchplusplus_lowRankMedNoise_chebyshev
-- compiled random input { i64 10i64 100i64 9i64 }
-- compiled random input { i64 90i64 100i64 99i64 }
-- output { true }
entry test_hutchplusplus_lowRankMedNoise_chebyshev seed R n s =
  HPPC.test (T.lowRankMedNoise (seed + 1) R n) s t seed

-- ==
-- entry: test_hutchplusplus_lowRankHiNoise_chebyshev
-- compiled random input { i64 10i64 100i64 9i64 }
-- compiled random input { i64 90i64 100i64 99i64 }
-- output { true }
entry test_hutchplusplus_lowRankHiNoise_chebyshev seed R n s =
  HPPC.test (T.lowRankHiNoise (seed + 1) R n) s t seed
