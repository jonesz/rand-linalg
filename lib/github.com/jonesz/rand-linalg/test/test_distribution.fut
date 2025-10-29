--| ignore
import "../distribution"
import "../../cbrng-fut/cbrng"

module SSD = mk_random_sparse_signs f32 u32 squares32

-- ==
-- entry: test_sparse_sign_distribution
-- random input { i64 1000000i64 }
-- output { [0.75f32, 0.125f32, 0.125f32] }
entry test_sparse_sign_distribution seed sz =
	let xs = tabulate sz (SSD.rand (squares32.construct seed) 0.25f32)
	let zero_cnt = map (\x_i -> x_i == 0f32)  xs |> map (i64.bool) |> i64.sum |> f32.i64
	let none_cnt = map (\x_i -> x_i == 1f32)  xs |> map (i64.bool) |> i64.sum |> f32.i64
	let pone_cnt = map (\x_i -> x_i == -1f32) xs |> map (i64.bool) |> i64.sum |> f32.i64
	in [zero_cnt / (f32.i64 sz), none_cnt / (f32.i64 sz), pone_cnt / (f32.i64 sz)]
