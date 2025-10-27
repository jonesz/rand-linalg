-- | Randomized SVD.
import "sketch"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/linalg"
import "../../diku-dk/linalg/qr"

local module type svd = {
	type t
	val svd [m][n] : (A: [m][n]t) -> ([m][m]t, [m][n]t, [n][n]t)
}

-- https://www.irisa.fr/sage/bernard/publis/SVD-Chapter06.pdf
-- https://web.stanford.edu/class/cme335/lecture7.pdf
local module two_sided_jacobi_serial(R: real): svd with t = R.t = {
	module L = mk_linalg R
	type t = R.t

	-- Compute the rotation matrix for `J(p, q, theta)`.
	let jacobi_rotation_matrix (p: i64) (q: i64) theta (n: i64) =
		let J = L.eye n

		let c = R.cos theta
		let s = R.sin theta

		let J = J with [p][q] = s
		let J = J with [p][p] = c          -- [[ c, s]
		let J = J with [q][p] = (R.neg s)  --  [-s, c]]
		let J = J with [q][q] = c

		in J

	def serial_schedule n =
		-- TODO: This filter operation sucks; we know the number of indices per sweep too.
		map (\j -> map (\k -> (j, k)) (iota n)) (iota n) |> flatten |> filter (\(j, k) -> j > k)

	def svd [m][n] (A: [m][n]t): ([m][m]t, [m][n]t, [n][n]t) =
		let U = L.eye n
		let V = L.eye n

		let t_k a_k =
			-- sign(a_k) / (|a_k| + sqrt(1 + a_k^2))
			(R.sgn a_k) / (R.abs a_k |> (R.+) (R.sqrt ((R.+) (R.i64 1) (R.* a_k a_k))))

		let S = L.matmul (transpose J_l) S
		let S = L.matmul S J_r

		let U = L.matmul U J_l
		let V = L.matmul V J_r

		-- TODO: Post-processing.
		in (U, S, V)
}

-- module randomized_svd (R: real) = {
-- 	let t = R.t
-- 
-- 	-- TODO: This requires a block size; it's perhaps better to implement MGS?
-- 	module QR = mk_block_householder R
-- 
-- 	def svd =
-- 		let Q = randomized_rangefinder f d B
-- 		let C = L.matmul (transpose Q) B
-- 		let (U, S, V) = svd C
-- 		let B_hat = L.matmul Q U |> flip (L.matmul) S |> flip (L.matmul) V
-- 		in ???
-- }
