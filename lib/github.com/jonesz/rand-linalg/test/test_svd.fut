-- | ignore

import "../svd"
import "tro_matrices"
import "../../../diku-dk/linalg/linalg"

module RSVD = mk_rsvd_default f32
module TM = mk_tro f32
module LA = mk_linalg f32

module OSJS = mk_one_sided_jacobi_serial f32

-- ==
-- entry: test_orthogonality_U_one_sided_serial_jacobi
-- compiled random input { [100][100]f32 }
-- output { 1.0f32 }
entry test_orthogonality_U_one_sided_serial_jacobi [l] (A: [l][l]f32) =
 let (U, _, _) = OSJS.svd A
 let U_T_U = LA.matmul (transpose U) U
 let I_l = LA.eye l
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten U_T_U) (flatten I_l) |> f32.sum

-- ==
-- entry: test_orthogonality_V_one_sided_serial_jacobi
-- compiled random input { [100][100]f32 }
-- output { true }
entry test_orthogonality_V_one_sided_serial_jacobi [l] (A: [l][l]f32) =
 let (_, _, V_T) = OSJS.svd A
 let V_T_V = LA.matmul V_T (transpose V_T)
 let I_l = LA.eye l
 in map2 (\a b -> (f32.-) a b |> f32.abs) (flatten V_T_V) (flatten I_l) |> f32.sum |> (f32.>=) 0.01f32

