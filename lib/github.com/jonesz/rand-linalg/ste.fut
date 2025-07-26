-- | Routines for computing a stochastic trace estimation.
import "../../diku-dk/linalg/linalg"
import "../../diku-dk/linalg/qr"

-- | Given a primitive `w -> Aw`, compute `AX = B`.
def matmul_w_matvec [n][m] 't (matvec: [n]t -> [n]t) (X: [n][m]t) =
    map (matvec) (transpose X) |> transpose

module type ste = {
  -- | The underlying scalar type of the matrix.
  type t

  -- | Estimate a trace.
  val ste [n] : (samples: i64) -> (cbrng_rand: i64 -> t) -> (matvec: [n]t -> [n]t) -> t
}

-- | The Monte Carlo Girard-Hutchinson trace estimator.
module hutchinson (R: real)
  : ste
    with t = R.t = {
  type t = R.t

  -- | Compute a single sample `X_i = w_i*(Aw_i)`.
  def sample cbrng_rand matvec i =
    -- TODO: There should be a conjugate transpose here.
    let dotprod a b =
      map2 (R.*) a b |> reduce (R.+) (R.i64 0)

    -- Compute the i'th random vector.
    let r_vec [n] i =
      map (\j -> i * n + j |> cbrng_rand) (iota n)
    
    let w = r_vec i
    -- `X = w*(Aw)`.
    in matvec w |> dotprod w
    
  -- | Form the Monte Carlo trace estimate `X_hat_s = 1/s * sum(X_i)`.
  def mcte [s] (samples: [s]t) =
    reduce (R.+) (R.i64 0) samples |> flip (R./) (R.i64 s)

  -- | Estimate the trace.
  def ste samples cbrng_rand matvec =
    map (sample cbrng_rand matvec) (iota samples) |> mcte
}

-- | The adaptive Hutch++ algorithm presented in https://arxiv.org/pdf/2010.09649
module hutchplusplus (R: real)
  : ste
    with t = R.t = {
  type t = R.t

  module L = mk_linalg R
  module QR = mk_gram_schmidt R

  -- Generate a random vector for a specific iteration.
  local def r_vec 't [n] (rand: i64 -> t) (iter: i64) : [n]t =
    map (\i -> (iter * n) + i |> rand) (iota n)

  local def tr A =
    L.fromdiag A |> reduce (R.+) (R.i64 0)

  def ste [n] (m: i64) (rand: i64 -> t) (matvec: [n]t -> [n]t) =
    let mot = m / 3
    -- "m_over_three"
    let z: [n + n][mot]t = map (r_vec rand) (iota (n + n))
    let S: [n][mot]t = take n z :> [n][mot]t
    let G: [n][mot]t = drop n z :> [n][mot]t

    -- A: [n][n]t; S: [n][mot]t; AS = [n][mot]t
    let AS: [n][mot]t = matmul_w_matvec matvec S
    let (Q: [n][n]t, _) = QR.qr AS

    -- tr(Q^T AQ)
    let left = map (matvec) Q |> L.matmul (transpose Q) |> tr
    -- (3/m) tr(G^T(I-QQ^T)A(I-QQ^T)G)
    let right =
      let I_QQT = L.matmul Q (transpose Q) |> L.matsub (L.eye n)
      in L.matmul I_QQT G                  -- I_QQ^T G
        |> matmul_w_matvec matvec          -- A (I-QQ^T) G
        |> L.matmul I_QQT                  -- (I-QQ^T) A (I-QQ^T) G
        |> L.matmul (transpose G)          -- G^T (I-QQ^T) A (I-QQ^T) G
        |> tr
        |> (R.*) (3_f32 / (f32.i64 m) |> R.f32)
    in (R.+) left right
}
