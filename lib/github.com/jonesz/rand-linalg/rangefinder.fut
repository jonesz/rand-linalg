-- | Solving the randomized rangefinder problem.
import "sketch"
import "qr/qr"
import "../cbrng-fut/distribution"

module type rangefinder = {
  -- | The underlying scalar type.
  type t 
  -- | The underlying distribution the rangefinder samples a test matrix \Omega from.
  module dist: cbrng_distribution

  val rangefinder [m] [n] : dist.engine.k -> (l : i64) -> [m][n]t -> [m][l]t
}

module type rangefinder_oracle = {
  -- | The underlying scalar type.
  type t 
  -- | The underlying distribution the rangefinder samples a test matrix \Omega from.
  module dist: cbrng_distribution

  val rangefinder [m] [n] : dist.engine.k -> (l : i64) -> ([n]t -> [m]t) -> [m][l]t
}

-- | A randomized rangefinder that computes `Q` for the rangefinder problem.
module mk_rangefinder
  (R: real) (D: sketch_right_dense with t = R.t)
  (QR: { val qr [m] [n] : [m][n]R.t -> ([m][n]D.t, [n][n]R.t) })
  : rangefinder with t = R.t = {

  type t = R.t
  module dist = D.dist

  def rangefinder [m] [n] seed l (B: [m][n]t) =
    let Y: [m][l]t = D.sketch seed l B
    let (q, _) = QR.qr Y
    in q
}

-- | A randomized rangefinder that computes `Q` for the rangefinder problem via a matrix-mul oracle.
module mk_rangefinder_oracle
  (R: real) (D: sketch_right_oracle with t = R.t)
  (QR: { val qr [m] [n] : [m][n]R.t -> ([m][n]D.t, [n][n]R.t) })
  : rangefinder_oracle with t = R.t = {

  type t = D.t
  module dist = D.dist

  def rangefinder [m] seed l oracle =
    let Y: [m][l]t = D.sketch seed l oracle
    let (q, _) = QR.qr Y
    in q
}

-- | The default randomized rangefinder, utilizing the `householder_thin_qr`.
module mk_rangefinder_dense_default (R: real) (D: sketch_right_dense with t = R.t)
  = mk_rangefinder R D
    {
      module QR = mk_householder_thin_qr R
      def qr = QR.qr ()
    }
