.PHONY: test clean bench

test: test_qr.fut qr.fut
	futhark test $<

bench: bench_qr.fut qr.fut
	futhark bench $<

clean:
	$(RM) *.actual *.expected qr.c qr test_qr.c test_qr bench_qr.c bench_qr
