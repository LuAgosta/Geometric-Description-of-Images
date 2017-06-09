GCC_FLAGS = -g -std=c99 -O3 -D_GNU_SOURCE
LINK_FLAGS = -ltiff -ljpeg -lpng -lm 

all: flst flstRGB

fsaddles.o: fsaddles.c
	gcc $(GCC_FLAGS) -c $^


flst.o: flst.c
	gcc $(GCC_FLAGS) -c $^

flst_reconstruct.o: flst_reconstruct.c
	gcc $(GCC_FLAGS) -c $^

flst_pixels.o: flst_pixels.c
	gcc $(GCC_FLAGS) -c $^

flst_boundary.o: flst_boundary.c
	gcc $(GCC_FLAGS) -c $^



flst_bilinear.o: flst_bilinear.c
	gcc $(GCC_FLAGS) -c $^

flstb_dual.o: flstb_dual.c
	gcc $(GCC_FLAGS) -c $^

flstb_dualchain.o: flstb_dualchain.c
	gcc $(GCC_FLAGS) -c $^

flstb_boundary.o: flstb_boundary.c
	gcc $(GCC_FLAGS) -c $^

flstb_quantize.o: flstb_quantize.c
	gcc $(GCC_FLAGS) -c $^




iio.o: iio.c 
	gcc $(GCC_FLAGS) -c $^

mw_error.o: mw_error.c 
	gcc $(GCC_FLAGS) -c $^

shape.o: shape.c 
	gcc $(GCC_FLAGS) -c $^

fsignal.o: fsignal.c 
	gcc $(GCC_FLAGS) -c $^

fimage.o: fimage.c
	gcc $(GCC_FLAGS) -c $^

list.o: list.c
	gcc $(GCC_FLAGS) -c $^

main_flst.o: main_flst.c
	gcc $(GCC_FLAGS) -c $^

flstRGB.o: flstRGB.c
	gcc $(GCC_FLAGS) -c $^



flst: flst_boundary.o flst_pixels.o flst.o flst_reconstruct.o  fsaddles.o \
      flstb_dualchain.o flstb_boundary.o flstb_dual.o flst_bilinear.o \
      flstb_quantize.o \
      fsignal.o mw_error.o shape.o fimage.o main_flst.o iio.o list.o 
	gcc $(GCC_FLAGS)  -o flst $^ $(LINK_FLAGS)


flstRGB: flst_boundary.o flst_pixels.o flst.o flst_reconstruct.o  fsaddles.o \
      flstb_dualchain.o flstb_boundary.o flstb_dual.o flst_bilinear.o \
      flstb_quantize.o  \
      fsignal.o mw_error.o shape.o fimage.o flstRGB.o iio.o list.o 
	gcc $(GCC_FLAGS)  -o flstRGB $^ $(LINK_FLAGS)



clean:
	rm -f *.o flst flstRGB *.nd *~ *# 

.PHONY: clean 

