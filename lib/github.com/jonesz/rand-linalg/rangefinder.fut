-- | Solving the randomized rangefinder problem.
import "sketch"
import "qr/qr"
import "../cbrng-fut/distribution"

module type rangefinder = {
  type t 
  module dist: cbrng_distribution
  val rangefinder [m] [n] : dist.engine.k -> (l : i64) -> [m][n]t -> [m][l]t
}

module type rangefinder_oracle = {
  type t 
  module dist: cbrng_distribution
  val rangefinder [m] [n] : dist.engine.k -> (l : i64) -> ([n]t -> [m]t) -> [m][l]t
}

-- | A randomized rangefinder that computes `Q` for the rangefinder problem.
-- TODO: Parameterize the Rangefinder with the econ QR.
-- module mk_rangefinder (D: sketch_right_dense) (QR: {val qr [m] [n] : [m][n]D.t -> ([m][n]D.t, [n][n]D.t)})
module mk_rangefinder (R: real) (D: sketch_right_dense with t = R.t)
  : rangefinder with t = R.t = {
  type t = R.t
  module dist = D.dist

  local module QR = mk_householder_thin_qr R

  def rangefinder [m] [n] seed l (B: [m][n]t) =
    let Y: [m][l]t = D.sketch seed l B
    let (q, _) = QR.qr () Y
    in q
}
