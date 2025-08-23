-- | ignore

import "../embedding"
import "../test_matrices"
import "../../cbrng-fut/cbrng"

module T = tm f32
module G = mk_gaussian_embedding f32 u32 squares32

-- ==
-- entry: test_embedding_gaussian
entry test_embedding_gaussian = 
	let m = T.lowRankHiNoise 0x1337 50i64 100i64
	in G.dense.embed (squares32.construct 0x1338) 30i64 m
