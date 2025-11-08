-- | Randomized SVD.
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"
import "rangefinder"
import "qr/qr"

-- | Randomized SVD.
module type rsvd = {
  -- | The underlying scalar type.
  type t
  -- | The underlying distribution this sketches from.
  module dist: cbrng_distribution

  -- | Compute the rsvd for a matrix A with `m >> n`.
  -- the argument `l` is the target rank `k` + `p` where `p` could be {k, 5, 10}, etc.
  val rsvd [m] [n] : dist.engine.k -> [m][n]t -> (l: i64) -> ([m][l]t, [l][l]t, [l][n]t)
}

-- | Compute the SVD via one-sided Jacobi iterations.
module mk_one_sided_jacobi_slow (R: real) : {
  type t
  val svd [l] : [l][l]t -> ([l][l]t, [l][l]t, [l][l]t)
} with t = R.t = {
	type t = R.t
	module LA = mk_linalg R

	-- Compute a 2x2 Schur Decomposition (8.5.2).
	-- https://math.ecnu.edu.cn/~jypan/Teaching/books/2013%20Matrix%20Computations%204th.pdf
	let schurr_decomp A_k (i: i64, j: i64) : (t, t) =
		let col_i = A_k[0:, i]
		let col_j = A_k[0:, j]

		-- form b_ij, b_jj, b_ii
		let (b_ij, b_jj, b_ii) =
			let tmp = map2 (LA.dotprod) [col_i, col_j, col_i] [col_j, col_j, col_i]
			in (tmp[0], tmp[1], tmp[2])

		-- TODO: Stopping criteria should be based on this condition, not number of sweeps.
		-- This is the convergence formula out of the `utk.edu` Dongarra PDF.
		let cond =
			let eps = R.f32 1e-4f32
			in (R.>=) (R.abs b_ij) ((R.*) b_jj b_ii |> R.sqrt |> (R.*) eps)

		-- TODO: Drop this conditional entirely? If they are orthogonal, do we just get the identity matrix?
		in if cond
		then
			let tau = (R.-) b_jj b_ii |> flip (R./) ((R.*) (R.i64 2) b_ij)
			-- TODO: This is `sgn` implementation which deviates from Golub.
			let t   = (R.*) tau tau |> (R.+) (R.i64 1) |> R.sqrt |> (R.+) (R.abs tau) |> (R./) (R.sgn tau)
			let cs  = (R.*) t t |> (R.+) (R.i64 1) |> R.sqrt |> (R./) (R.i64 1)
			let sn  = (R.*) cs t
			in (cs, sn)
		else
			(R.i64 1, R.i64 0)

	-- A chess round-robin scheduling algorithm; fix the head player, rotate, and then
	-- pair each player at the ends.
	let round_robin [l] (xs: [l]i64) : [l]i64 = [0i64] ++ rotate 1 (tail xs) :> [l]i64
	let pairs [l] (xs: [l]i64) : [l/2i64](i64, i64) =
		map (\i ->
			-- `0` should pair with `l`, `1` should pair with `l-1`, so forth...
			let a = drop i xs |> head
			let b = reverse xs |> drop i |> head
			in (a, b)
		) (iota (l / 2i64))

	-- Compute a single parallel sweep over the entirety of schedule.
	let parallel_sweep [l] (A_k: [l][l]t) (V_k: [l][l]t) =
		let xs = iota l
		let (A_k, V_k, _) = 
			loop (A_k, V_k, xs) = (A_k, V_k, xs) for _cycle < l do
				let col_pairs = pairs xs
				let rotations = map (schurr_decomp A_k) col_pairs

				let J =
					loop J = (LA.eye l) for ((cs, sn), (p, q)) in (zip rotations col_pairs) do 
						J with [p, p] = cs
					  	  with [p, q] = sn
					  	  with [q, p] = sn |> R.neg
					  	  with [q, q] = cs

				let (A_k, V_k) =
					let tmp = map (flip (LA.matmul) J) [A_k, V_k]
					in (tmp[0], tmp[1])

				in (A_k, V_k, round_robin xs)

		in (A_k, V_k)

    -- | Compute (U, S, V_T).
	def svd [l] (A: [l][l]t) : ([l][l]t, [l][l]t, [l][l]t) =
		let (A_hat, V) =
			-- TODO: See the note about convergence above. 10~ sweeps should be enough.
			loop (A_k, V_k) = (A, LA.eye l) for _iter < 10 do
				parallel_sweep A_k V_k

        -- The singular values are the euclidean norms of each column.
		let S = transpose A_hat |> map (LA.vecnorm) |> LA.todiag

        -- Compute U; (A_hat = US) == (A_hat * S_inv = U)
		let U =
			let S_inv = LA.fromdiag S |> map (\x -> (R./) (R.i64 1) x) |> LA.todiag
			in LA.matmul A_hat S_inv

		in (U, S, (transpose V))
}

module mk_rsvd (R: real) (T: rangefinder with t = R.t) : rsvd with t = R.t = {
  module dist = T.dist
  type t = R.t

  module LA = mk_linalg R
  module TQ = mk_householder_thin_qr R
  module SVD = mk_one_sided_jacobi_slow R

  -- Retun (U, S, V^T); Assuming that m >> n.
  def rsvd [m] [n] seed (A: [m][n]t) (l: i64) =

    -- Compute `Q` with the rangefinder.
    let Q: [m][l]t = T.rangefinder seed l A
    -- Form `C = Q*A`, the approximation; note: l <= n.
    let C: [l][n]t = LA.matmul (transpose Q) A

    -- Following the randomized rangefinder problem `l << n`...
    -- See [Dongarra](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf);
    -- "if m << n, it is more efficient to perform an LQ factorization of A".
    -- LQ can be achieved through our `thin_qr` and a transpose.
    let (Q_qr, R_qr) : ([n][l]t, [l][l]t) = transpose C |> TQ.qr ()
    let (Q_qr, L_qr) : ([l][n]t, [l][l]t) = (transpose Q_qr, transpose R_qr)

    let (U_hat, S, V_T) = SVD.svd L_qr

    -- Revert the LQ transformation.
    let V = LA.matmul V_T Q_qr
    -- Form the SVD..
    let U = LA.matmul Q U_hat

    in (U, S, V)
}

import "../cbrng-fut/cbrng"
import "sketch"

-- | RSVD with Gaussian embeddings and the default rangefinder.
module mk_rsvd_default (R: real) : rsvd with t = R.t = {
  local module EMB = mk_gaussian_embedding R u32 squares32
  local module RSVD = mk_rsvd R (mk_rangefinder_dense_default R EMB)

  type t = R.t
  module dist = RSVD.dist

  def rsvd = RSVD.rsvd
}

