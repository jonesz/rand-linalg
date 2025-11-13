FUTHARK?=futhark

TESTS=test_*.fut

.PHONY: test
test: $(TESTS)
	$(FUTHARK) test $(TESTS)

.PHONY: clean
clean:
	$(RM) test_distribution test_qr test_sketch_equivalence test_ste_chebyshev test_ste_eps test_svd test_turnstile_svd *.c *.actual *.expected
