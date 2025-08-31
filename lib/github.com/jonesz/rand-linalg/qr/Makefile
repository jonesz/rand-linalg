.PHONY: test clean bench

test: test_tsqr.fut tsqr.fut
	futhark test $<

bench: bench_qr_econ.fut tsqr.fut
	futhark bench $<

clean:
	$(RM) *.actual *.expected tsqr.c tsqr test_tsqr.c test_tsqr bench_qr_econ.c bench_qr_econ
