-- | ignore

import "../svd"
import "tro_matrices"
import "../../../diku-dk/linalg/linalg"

module LA = mk_linalg f32
module OSJS = mk_one_sided_jacobi f32
module RSVD = mk_rsvd_default f32

-- ==
-- entry: test_orthogonality_U_one_sided_jacobi
-- compiled random input { [10][10]f32 }
-- output { true }
entry test_orthogonality_U_one_sided_jacobi [l] (A: [l][l]f32) =
 let (U, _, _) = OSJS.svd A
 let U_T_U = LA.matmul (transpose U) U
 let I_l = LA.eye l
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten U_T_U) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_orthogonality_V_one_sided_jacobi
-- compiled random input { [10][10]f32 }
-- output { true }
entry test_orthogonality_V_one_sided_jacobi [l] (A: [l][l]f32) =
 let (_, _, V_T) = OSJS.svd A
 let V_T_V = LA.matmul V_T (transpose V_T)
 let I_l = LA.eye l
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten V_T_V) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_reconstruction_one_sided_jacobi
-- compiled random input { [10][10]f32 }
-- output { true }
entry test_reconstruction_one_sided_jacobi A =
 let (U, S, V_T) = OSJS.svd A
 let A_reconstructed = LA.matmul (LA.matmul U S) V_T
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten A) (flatten A_reconstructed) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_orthogonality_U_rsvd
-- compiled random input { i64 [500][100]f32 }
-- output { true } 
entry test_orthogonality_U_rsvd [m] [n] seed (A: [m][n]f32) =
 let seed = RSVD.dist.engine.construct seed
 let (U, _, _) = RSVD.rsvd seed A 10i64
 let U_T_U = LA.matmul (transpose U) U
 let I_l = LA.eye 10i64
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten U_T_U) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_orthogonality_V_rsvd
-- compiled random input { i64 [500][100]f32 }
-- output { true }
entry test_orthogonality_V_rsvd [m] [n] seed (A: [m][n]f32) = 
 let seed = RSVD.dist.engine.construct seed
 let (_, _, V_T) = RSVD.rsvd seed A 10i64
 let V_T_V = LA.matmul V_T (transpose V_T)
 let I_l = LA.eye 10i64
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten V_T_V) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

-- ==
-- entry: test_reconstruction_rsvd
-- compiled random input { i64 [500][100]f32 }
-- output { true } 
entry test_reconstruction_rsvd [m] [n] seed (A: [m][n]f32) =
 let f_norm X = flatten X |> map (\x -> x * x) |> f32.sum |> f32.sqrt -- Frobenius Norm.
 let seed = RSVD.dist.engine.construct seed
 let (U, S, V_T) = RSVD.rsvd seed A 10i64
 let A_reconstructed = LA.matmul (LA.matmul U S) V_T
 in
  let top = map2 (\a b -> (f32.-) a b |> f32.abs) (flatten A) (flatten A_reconstructed) |> unflatten |> f_norm
  let bot = f_norm A
  in (f32./) top bot |> (f32.>=) 1.0_f32
 
