-- | ignore

import "tsqr"
import "../test_matrices"
import "../../../diku-dk/linalg/linalg"

module QR_ECON = qr_econ f32
module TM = tm f32
module L = mk_linalg f32

-- | All elements squared, summed, then take the square root.
local def f_norm A =
	let A = flatten A
	in map2 (*) A A |> reduce (+) 0 |> f32.sqrt
	
-- ==
-- entry: test_econ_qr
-- compiled random input { [100][1000000][100]f32 }
-- output { true }
entry test_econ_qr rm =
	let f A =
		let (q, r) = QR_ECON.qr A
		let A_norm = f_norm A
		let A_minus_QR_norm = L.matmul q r |> L.matsub A |> f_norm
		in A_minus_QR_norm / A_norm |> (<=) 0.001_f32
	in map (f) rm
