SRC= \
	lib/github.com/jonesz/rand-linalg/ste.fut \
	lib/github.com/jonesz/rand-linalg/test_matrices.fut

TEST = \
	lib/github.com/jonesz/rand-linalg/test/hutchinson_eps_diagonal.fut

all: test

test: $(TEST) $(SRC)
	futhark test --backend=multicore $<

.PHONY: doc clean

doc:
	futhark doc -o doc/ lib/github.com/jonesz/rand-linalg

clean:
	$(RM) -rf doc/
	$(RM) lib/github.com/jonesz/rand-linalg/test/*.c
	$(RM) lib/github.com/jonesz/rand-linalg/test/*.expected
	$(RM) lib/github.com/jonesz/rand-linalg/test/*.actual
	$(RM) lib/github.com/jonesz/rand-linalg/test/hutchinson_eps_diagonal

