-- | Randomized SVD.
import "rangefinder"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"
import "qr/qr"

local module one_sided_jacobi_serial(R: real) = {
	type t = R.t

	module TQ = mk_householder_thin_qr R
	module LA = mk_linalg R

	def serial_schedule n =
		-- TODO: This filter operation sucks; we know the number of indices per sweep too.
		map (\j -> map (\k -> (j, k)) (iota n)) (iota n) |> flatten |> filter (\(j, k) -> j > k)

	def sweep [l] (U: [l][l]t) (S: [l][l]t) (V: [l][l]t) : ([l][l]t, [l][l]t, [l][l]t) =
		let svd_2v2_rot a_pp a_pq a_qp a_qq =
			-- [2x2 Symmetric Schur Decomposition](https://web.stanford.edu/class/cme335/lecture7.pdf)
			-- We aren't symmetric, so the denominator is a + d.
			let (cr, sr) =
				let top = a_qq - a_pp
				let bot = a_pq + a_qp
				in if (top == (R.i64 0) && bot == (R.i64 0))
					then ((R.i64 1), (R.i64 0))
					else
						let tau = (R./) top bot
						let t = ???
						let c = 1 / sqrt(1 + t^2) 
						let s = c * t
			
		foldl (\(U_i, S_i, V_i) (i, j) ->
			let a_ii = S[i][i]
			let a_ij = S[i][j]
			let a_ji = S[j][i]
			let a_jj = S[j][j]
			in ???)
		(U, S, V) (serial_schedule l)

	-- Return (U, S, V^T)
	def svd_econ [l][n] (A: [l][n]t): ([l][l]t, [l][l]t, [l][n]t) =

		-- Following the randomized rangefinder problem `l << n`; we need a square
		-- matrix before the two-sided jacobi iteration: compute the LQ factorization.
		-- See [Dongarra](https://icl.utk.edu/files/publications/2018/icl-utk-1341-2018.pdf);
		-- "if m << n, it is more efficient to perform an LQ factorization of A".
		-- LQ can be achieved through our `thin_qr` and a transpose.
		--
		-- There's a chance where stability doesn't matter and A * A^T is reasonable.

		let (Q, R) = transpose A |> TQ.qr ()     -- `thin_qr([n][l]t)``
		let (Q, L) = (transpose Q, transpose R)

		-- At this point, we compute the SVD on the matrix `L`, which is `[l][l]t` (square)
		-- and suitable for a two-sided jacobi.

		-- TODO: We should sweep based on some evaluation of the off-diagonal elements
		-- and their convergence to zero... instead we'll sweep ten times.
		let (U, S, V)
			= loop (U, S, V) = (LA.eye l, L, LA.eye l) for _sweep_num < 10 do
				sweep U S V

		-- TODO: Sort the singular values?
		-- TODO: Change all negative singular values to positive.
		-- TODO: Determine whether this construction of V is correct.
		let V = LA.matmul (transpose V) Q
		in (U, S, V)

}

module type svd = {
	type t
	type seed
	val svd [m][n] : seed -> [m][n]t -> i64 -> ()
}

module randomized_svd (R: real) (T: rangefinder with t = R.t) : svd with t = R.t = {
	type t = R.t
	type seed = T.sk.dist.engine.k

	module L = mk_linalg R

	-- Assuming that m >> n.
	def svd [m][n] seed (A: [m][n]t) (l: i64) =
		-- Form `Q = qr_econ(AQ)`. 
		let Q: [m][l]t = T.rangefinder seed l A
		-- Form `C = Q*A`, the approximation; note: l <= n.
		let C: [l][n]t = L.matmul (transpose Q) A

		let (U_hat, S, V): ([l][l]t, [l][m]t, [m][m]t) = ???
		let U = L.matmul Q U_hat
		-- in (U, S, V)
		in ???
}
