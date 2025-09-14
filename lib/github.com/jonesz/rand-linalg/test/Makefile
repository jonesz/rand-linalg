test: test_sketch_equivalence

test_sketch_equivalence: sketch_equivalence.fut
	futhark test $^

.PHONY: clean 

clean: rm_sketch_equivalence

rm_sketch_equivalence:
	$(RM) sketch_equivalence sketch_equivalence.c
