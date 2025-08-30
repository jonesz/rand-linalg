test: test_tsqr.fut tsqr.fut
	futhark test $<

.PHONY: clean

clean:
	$(RM) *.actual *.expected tsqr.c tsqr test_tsqr.c test_tsqr bench_qr_econ.c bench_qr_econ
