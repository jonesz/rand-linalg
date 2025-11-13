-- | ignore

import "../qr/qr"
import "../../../diku-dk/linalg/linalg"

module H = mk_householder_thin_qr f32
module L = mk_linalg f32

local def tol = 1e-4f32

-- Frobenius Norm: ||A||_F = sqrt(tr(A^TA))
local def f_norm A =
	L.matmul (transpose A) A |> L.fromdiag |> reduce (f32.+) 0f32 |> f32.sqrt
	
-- ==
-- entry: householder_thin_qr_reconstruction
-- compiled random input { [100][1000][10]f32 }
-- output { true }
entry householder_thin_qr_reconstruction As =
	-- Check `||A - QR||_F / ||A||_F < tol`.
	let f A =
		let (Q, R) = H.qr () A
		let A_F    = f_norm A
		let A_QR_F = L.matmul Q R |> L.matsub A |> f_norm
		let tmp = A_QR_F / A_F
		in tol > tmp

	in map f As |> and

-- ==
-- entry: householder_thin_qr_orthogonality
-- compiled random input { [100][1000][10]f32 }
-- output { true }
entry householder_thin_qr_orthogonality [m][n][p] (As: [m][n][p]f32) =
	-- Test `||I - Q^TQ||_F < tol`.
	let f A =
		let (Q, R) = H.qr () A
		let norm = L.matmul (transpose Q) Q |> L.matsub (L.eye p) |> f_norm
		in tol > norm
	in map f As |> and
