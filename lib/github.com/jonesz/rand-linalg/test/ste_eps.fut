-- | ignore
-- Check whether the STE is within `eps` of the actual trace; the majority
-- of these tests are on diagonal matrices.

import "../ste"
import "../../cbrng-fut/distribution"
import "../../cbrng-fut/cbrng"
import "../../../diku-dk/linalg/linalg"
import "tro_matrices"

module mk_eps_test (R: real) (S: ste with t = R.t) = {
  module L = mk_linalg R
  module D = rademacher_distribution R u32 squares32

  -- A constant `eps`, for which the stochastic estimator is allowed to deviate.
  local def eps = R.f32 0.001_f32
  local def tr A = L.fromdiag A |> reduce (R.+) (R.i64 0)

  -- Return the difference between an estimated trace and the actual trace.
  def tr_diff seed A k =
    let seed = squares32.construct seed
    let exp_tr = tr A
    let est_tr = S.ste k (D.rand seed ()) (L.matvecmul_row A)
    in (R.-) exp_tr est_tr |> R.abs

  def test seed A k = tr_diff seed A k |> (R.>) eps
}

module T = mk_tro f32

-- The Hutchinson estimator should require a single sample.
module H = mk_eps_test f32 (hutchinson f32)

-- ==
-- entry: test_hutchinson_polyDecaySlow_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecaySlow_eps seed R n s =
  let A = T.polyDecaySlow R n
  in H.test seed A s

-- ==
-- entry: test_hutchinson_polyDecayMed_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayMed_eps seed R n s =
  let A = T.polyDecayMed R n
  in H.test seed A s

-- ==
-- entry: test_hutchinson_polyDecayFast_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_polyDecayFast_eps seed R n s =
  let A = T.polyDecayFast R n
  in H.test seed A s

-- ==
-- entry: test_hutchinson_expDecaySlow_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecaySlow_eps seed R n s =
  let A = T.expDecaySlow R n
  in H.test seed A s

-- ==
-- entry: test_hutchinson_expDecayMed_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayMed_eps seed R n s =
  let A = T.expDecayMed R n
  in H.test seed A s

-- ==
-- entry: test_hutchinson_expDecayFast_eps
-- compiled random input { i64 10i64 100i64 1i64 }
-- compiled random input { i64 90i64 100i64 1i64 }
-- output { true }
entry test_hutchinson_expDecayFast_eps seed R n s =
  let A = T.expDecayFast R n
  in H.test seed A s
