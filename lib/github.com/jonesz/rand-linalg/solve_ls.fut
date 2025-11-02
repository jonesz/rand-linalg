-- | Solving linear systems via sketching.
import "sketch"

-- | Solve a linear system via Sketch-and-Solve.
module mk_SAS (R: real) (S: sketch_left_dense with t = R.t) = {
	type t = R.t
	module dist = S.dist

	def SAS [m][n] seed (d: i64) (A: [m][n]t) (b: [m]t) (solver: [d][n]t -> [d]t -> [n]t) =
		let SA = S.sketch seed d A 
		let Sb = S.sketch seed d (transpose [b]) |> flatten :> [d]t
		in solver SA Sb
}
