-- | Solving the randomized rangefinder problem.
import "sketch"

-- | A randomized rangefinder that computes `Q` for the rangefinder problem.
module mk_rangefinder (D: sketch_right_dense) (QR: { val qr [m][n] : [m][n]D.t -> ([m][n]D.t, [n][n]D.t) }): {
	module sk: sketch_right_dense
	val rangefinder [m][n] : sk.dist.engine.k -> (l: i64) -> [m][n]D.t -> [m][l]D.t
} = {
	module sk = D
	def rangefinder [m][n] seed l (B: [m][n]D.t) =
		let Y: [m][l]D.t = D.sketch seed l B
		let (q, _) = QR.qr Y
		in q
}

-- | A randomized rangefinder that computes `Q` for the rangefinder problem, with a subspace iteration before.
module mk_rangefinder_subspace (D: sketch_right_dense) (QR: { val qr [n][m] : [m][n]D.t -> ([m][n]D.t, [n][n]D.t) }): {
	module sk: sketch_right_dense
	val rangefinder [m][n] : sk.dist.engine.k -> (l: i64) -> (q: i64) -> [m][n]D.t -> [m][l]D.t
} = {
	module sk = D
	def rangefinder seed l q B = ???
}
