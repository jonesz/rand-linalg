-- | ignore
import "../svd"
import "../rangefinder"
import "../sketch"
import "../../cbrng-fut/cbrng"
import "tro_matrices"

module G = mk_gaussian_embedding f32 u32 squares32
module R = mk_rangefinder f32 G.dense.right
module SVD = randomized_svd f32 R

-- -- ==
-- -- entry: test_linear_system_via_svd
-- -- compiled random input { i64 }
-- -- output { true }
-- entry test_linear_system_via_svd seed =
-- 	let 
