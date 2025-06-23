all: doc

.PHONY: doc clean

doc:
	futhark doc -o doc/ lib/github.com/jonesz/rand-linalg

clean:
	rm -rf doc/
