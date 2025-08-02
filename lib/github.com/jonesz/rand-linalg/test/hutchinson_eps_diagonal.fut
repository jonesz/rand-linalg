-- | ignore
-- | Tests for the Hutchinson estimator, specifically those that expect the estimated `X_hat`
-- to be within `eps` of the true `tr(A)`.

import "../ste"
import "../test_matrices"
import "../../cbrng-fut/distribution"
import "../../cbrng-fut/cbrng"
import "../../../diku-dk/linalg/linalg"

module mk_eps_test (R: real) (S: ste with t = R.t) = {
  module L = mk_linalg R
  module D = rademacher_distribution R u32 i64 squares32

  -- | A constant `eps`, for which the stochastic estimator is allowed to deviate.
  def eps = R.f32 0.00001_f32

  def tr A = L.fromdiag A |> reduce (R.+) (R.i64 0)

  -- | Return the difference between an estimated trace and the actual trace.
  def tr_diff seed A samples =
    let seed = squares32.construct seed
    let exp_tr = tr A
    let est_tr = S.ste samples (D.rand seed) (L.matvecmul_row A)
    in (R.-) exp_tr est_tr |> R.abs

  def test seed A samples = tr_diff seed A samples |> (R.<) eps
}

module H = mk_eps_test f32 (hutchinson f32)
module T = tm f32

-- The following tests only require 1 sample because they're only diagonal matrices.

-- entry: test_hutchinson_polyDecaySlow_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecaySlow_eps seed R n s =
  let A = T.polyDecaySlow R n
  in H.test seed A s

-- entry: test_hutchinson_polyDecayMed_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayMed_eps seed R n s =
  let A = T.polyDecayMed R n
  in H.test seed A s

-- entry: test_hutchinson_polyDecayFast_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayFast_eps seed R n s =
  let A = T.polyDecayFast R n
  in H.test seed A s

-- entry: test_hutchinson_expDecaySlow_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecaySlow_eps seed R n s =
  let A = T.expDecaySlow R n
  in H.test seed A s

-- entry: test_hutchinson_expDecayMed_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayMed_eps seed R n s =
  let A = T.expDecayMed R n
  in H.test seed A s

-- entry: test_hutchinson_expDecayFast_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayFast_eps seed R n s =
  let A = T.expDecayFast R n
  in H.test seed A s
