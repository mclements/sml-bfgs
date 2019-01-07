
test-smlnj:
	sml @SMLquiet -m test-smlnj.cm

c:
	gcc `pkg-config --cflags libRmath` vmmin.c `pkg-config --libs libRmath`
	./a.out


