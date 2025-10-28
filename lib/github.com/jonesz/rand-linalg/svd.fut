-- | Randomized SVD.
import "rangefinder"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"

module type svd = {
	type t
	type seed 

	val svd [L][N] : seed -> (A: [L][N]t) -> (k: i64) -> ([L][k]t, [k][k]t, [N][k]t)
}

module randomized_svd (R: real) (T: rangefinder with t = R.t) : svd with t = R.t = {
	type t = R.t
	type seed = T.sk.dist.engine.k

	module L = mk_linalg R

	def svd [m][n] seed (A: [m][n]t) (k: i64) =
		-- Form `Q = qr_econ(AQ)`.
		let Q = T.rangefinder seed k A
		-- Form `C = Q*A`, the approximation.
		let C = L.matmul (transpose Q) A

		-- This should be a "small dense SVD".
		-- let (U_hat, S, V) = svd_compact C
		-- let U = L.matmul Q U_hat
		-- in (U, S, V)
		in ???
}
