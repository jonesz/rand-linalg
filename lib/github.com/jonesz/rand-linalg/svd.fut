-- | Randomized SVD.
import "sketch"
import "../cbrng-fut/distribution"
import "../../diku-dk/linalg/qr"

local module type svd = {
	type t
	val svd [m][n] : (A: [m][n]t) -> ([m][m]t, [m][n]t, [n][n]t)
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
