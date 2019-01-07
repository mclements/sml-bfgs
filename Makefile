
test-smlnj:
	sml @SMLquiet -m test-smlnj.cm

c:
	gcc `pkg-config --cflags libRmath` vmmin.c -lRmath
	./a.out


