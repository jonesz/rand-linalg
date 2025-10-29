-- | Solving the randomized rangefinder problem.

import "sketch"

module type rangefinder = {
  type t 
  module sk: sketch_right_dense
  val rangefinder [m] [n] : sk.dist.engine.k -> (l : i64) -> [m][n]t -> [m][l]t
}

module type rangefinder_oracle = {
  type t 
  module sk: sketch_right_oracle
  val rangefinder [m] [n] : sk.dist.engine.k -> (l : i64) -> ([n]t -> [m]t) -> [m][l]t
}

-- | A randomized rangefinder that computes `Q` for the rangefinder problem.
module mk_rangefinder (D: sketch_right_dense) (QR: {val qr [m] [n] : [m][n]D.t -> ([m][n]D.t, [n][n]D.t)})
  : rangefinder with t = D.t = {
  type t = D.t
  module sk = D

  def rangefinder [m] [n] seed l (B: [m][n]D.t) =
    let Y: [m][l]D.t = D.sketch seed l B
    let (q, _) = QR.qr Y
    in q
}

-- | A randomized rangefinder that computes `Q` for the rangefinder problem via a matrix-free oracle.
module mk_rangefinder_oracle (D: sketch_right_oracle) (QR: {val qr [m] [n] : [m][n]D.t -> ([m][n]D.t, [n][n]D.t)})
  : rangefinder_oracle with t = D.t = {
  type t = D.t
  module sk = D

  def rangefinder [m] [n] seed l (oracle: [n]t -> [m]t) =
    let Y: [m][l]D.t = D.sketch seed l oracle
    let (q, _) = QR.qr Y
    in q
}

-- -- | A randomized rangefinder that computes `Q` for the rangefinder problem, with a subspace iteration before.
-- module mk_rangefinder_subspace (D: sketch_right_dense) (QR: {val qr [n] [m] : [m][n]D.t -> ([m][n]D.t, [n][n]D.t)})
--   : {
--       module sk: sketch_right_dense
--       val rangefinder [m] [n] : sk.dist.engine.k -> (l: i64) -> (q: i64) -> [m][n]D.t -> [m][l]D.t
--     } = {
--   module sk = D
-- 
--   def rangefinder seed l q B = ???
-- }
