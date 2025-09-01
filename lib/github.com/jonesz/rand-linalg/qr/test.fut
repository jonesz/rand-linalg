def f k =
	loop (acc, B) = (1, []) for i < k do
		let acc = acc
		let B = map (\z -> map (\x -> z + k) (iota acc)) (iota acc)
		in (acc + 1, B)
