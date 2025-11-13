-- | ignore

import "../distribution"
import "../../cbrng-fut/cbrng"

module SSD = mk_random_sparse_signs f32 u32 squares32

-- TODO: I don't think the underlying `squares32` distrbution is uniform under
-- `[E.min, E.max]`.

-- ==
-- entry: test_sparse_sign_distribution
-- random input { i64 1000000i64 }
-- output { [0.75f32, 0.125f32, 0.125f32] }
entry test_sparse_sign_distribution seed sz =
	let cnt z =
		\xs -> map (f32.== z) xs |> map (i64.bool) |> i64.sum |> f32.i64

	let xs = tabulate sz (SSD.rand (squares32.construct seed) 0.25f32)
	let (z, n, o) =
		let tmp = map (\i -> cnt i xs) [0f32, -1f32, 1f32]
		in (tmp[0], tmp[1], tmp[2])

	let sz_f32 = f32.i64 sz
	in [z / sz_f32, n / sz_f32, o / sz_f32]
