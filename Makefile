all:
	./setup.py build_ext --inplace

.PHONY: clean
clean:
	rm -rf *.so *.pyc *.c build 