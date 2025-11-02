-- | ignore

import "../svd"
import "../rangefinder"
import "../sketch"
import "../../cbrng-fut/cbrng"
import "tro_matrices"
import "../../../diku-dk/linalg/linalg"

-- module G = mk_gaussian_embedding f32 u32 squares32
-- module R = mk_rangefinder f32 G.dense.right
-- module SVD = randomized_svd f32 R
-- module TM = mk_tro f32
-- module L = mk_linalg f32

-- TODO: Figure out how to test this. I've run through the REPL and it works fine.

-- -- ==
-- -- entry: test_linear_system_via_svd
-- -- compiled random input { i64 }
-- -- output { true }
-- entry test_linear_system_via_svd seed =
--  	let svd_seed = SVD.dist.engine.construct seed
-- 	let mat = TM.lowRankHiNoise (seed + 1) 50 100
-- 	let (U, S, V_T) = SVD.svd_econ svd_seed mat 30
-- 	let mat_reconstructed = L.matmul (L.matmul U S) V_T
-- 	in map2 (\a b -> a - b |> f32.abs) (flatten mat) (flatten mat_reconstructed) |> f32.sum |> (f32.>=) 0.01f32
