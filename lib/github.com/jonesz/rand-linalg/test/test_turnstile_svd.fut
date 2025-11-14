-- | ignore
import "../svd"
import "../../../diku-dk/linalg/linalg"

module LA = mk_linalg f32
module TSVD = mk_sketchy_svd_default f32 

-- ==
-- entry: test_turnstile_orthogonality_U
-- compiled random input { i64 [500][500]f32 }
-- output { true }
entry test_turnstile_orthogonality_U [l] seed (A: [l][l]f32) =
	let k = 20i64
	-- let r = 20i64
	let (U, _, _) =
		let (X, Y, Z, s) = TSVD.initialize l l k seed
		let (X, Y, Z) =
			loop (X_i, Y_i, Z_i) = (X, Y, Z) for i < l do
				TSVD.linear_update_col s X_i Y_i Z_i i A[0:, i]
		in TSVD.sketchy_svd s X Y Z -- r

	let U_T_U = LA.matmul (transpose U) U
	let I_l = LA.eye k
 	in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten U_T_U) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_turnstile_orthogonality_V
-- compiled random input { i64 [500][500]f32 }
-- output { true }
entry test_turnstile_orthogonality_V [l] seed (A: [l][l]f32) =
	let k = 20i64
	-- let r = 20i64
	let (_, _, V_T) =
		let (X, Y, Z, s) = TSVD.initialize l l k seed
		let (X, Y, Z) =
			loop (X_i, Y_i, Z_i) = (X, Y, Z) for i < l do
				TSVD.linear_update_col s X_i Y_i Z_i i A[0:, i]
		in TSVD.sketchy_svd s X Y Z -- r

	let V_T_V = LA.matmul V_T (transpose V_T)
	let I_l = LA.eye k
 	in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten V_T_V) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_turnstile_reconstruction
-- compiled random input { i64 [500][500]f32 }
-- output { true }
entry test_turnstile_reconstruction [l] seed (A: [l][l]f32) =
	let k = 20i64
	-- let r = 20i64
	let (U, S, V_T) =
		let (X, Y, Z, s) = TSVD.initialize l l k seed
		let (X, Y, Z) =
			loop (X_i, Y_i, Z_i) = (X, Y, Z) for i < l do
				TSVD.linear_update_col s X_i Y_i Z_i i A[0:, i]
		in TSVD.sketchy_svd s X Y Z -- r

	let A_reconstructed = LA.matmul (LA.matmul U S) V_T
	in 
 		let f_norm X = flatten X |> map (\x -> x * x) |> f32.sum |> f32.sqrt -- Frobenius Norm.
  		let top = map2 (\a b -> (f32.-) a b |> f32.abs) (flatten A) (flatten A_reconstructed) |> unflatten |> f_norm
  		let bot = f_norm A
  		in (f32./) top bot |> (f32.>=) 1.0_f32
