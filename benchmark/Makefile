all: bench_qr bench_svd

.PHONY: bench_qr bench_svd clean

bench_qr:
	futhark bench --backend=multicore qr.fut

bench_svd:
	futhark bench --backend=multicore svd.fut

clean:
	$(RM) *.c *.expected *.actual
	$(RM) qr svd
