-- | ignore

import "tsqr"
import "../test_matrices"
import "../../../diku-dk/linalg/linalg"

module QR_ECON = qr_econ f32
module TM = tm f32
module L = mk_linalg f32

def frobenius_norm A =
	map2 (f32.*) (flatten A) (flatten A) |> reduce (f32.+) 0f32 |> f32.sqrt
	
-- ==
-- entry: test_econ_qr
-- compiled random input { [100][1000000][100]f32 }
-- output { true }
entry test_econ_qr rm =
	map (\rm_i -> let (q, r) = QR_ECON.qr rm_i
		let f_norm = frobenius_norm rm_i
		let qr_norm = L.matmul q r |> L.matsub rm_i |> frobenius_norm 
		let ratio = qr_norm / f_norm 
		in ratio <= 0.001_f32) rm |> and
