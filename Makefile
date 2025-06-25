all: test

test: lib/github.com/jonesz/rand-linalg/ste_test.fut lib/github.com/jonesz/rand-linalg/dist_test.fut
	futhark test $^

.PHONY: doc clean

doc:
	futhark doc -o doc/ lib/github.com/jonesz/rand-linalg

clean:
	rm -rf doc/ lib/github.com/jonesz/rand-linalg/*.c lib/github.com/jonesz/rand-linalg/*.expected \
	lib/github.com/jonesz/rand-linalg/*.actual \
	lib/github.com/jonesz/rand-linalg/ste_test \
	lib/github.com/jonesz/rand-linalg/dist_test

