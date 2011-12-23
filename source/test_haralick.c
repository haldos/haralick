
//		+-----------------------------------+
//		| Haralick's Edge Detector          |
//		| Implemented by Haldo Spontón		|
//		| Last modified: 22/dic/2011		|
//		+-----------------------------------+

#include "iio.c"
#include "2dconvolution.c"
#include <time.h>


int main(int argc, char *argv[]) {
	if (argc != 4) {
		printf("Usage: %s input_image_1 rhozero output\n", argv[0]);
	} else {

		// Execution time:
		double start = (double)clock();
	
		// Parameters
		float rhozero = atof(argv[2]);
	
		// Load input image (using iio)
		int w, h, pixeldim;
		float *im_orig = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
		fprintf(stderr, "Input image loaded:\t %dx%d image with %d channel(s).\n", w, h, pixeldim);

		// Grayscale conversion (if necessary)
		double *im = malloc(w*h*sizeof(double));
		if (im == NULL){
			fprintf(stderr, "Out of memory...\n");
			exit(EXIT_FAILURE);
		}
		int z;
		int zmax = w*h;
		if (pixeldim==3){
			for(z=0;z<zmax;z++){
				im[z] =  (double)(6968*im_orig[3*z] + 23434*im_orig[3*z + 1] + 2366*im_orig[3*z + 2])/32768;
			}
			fprintf(stderr, "images converted to grayscale\n");
		} else {
			for(z=0;z<zmax;z++){
				im[z] = (double)im_orig[z];
			}
			fprintf(stderr, "images are already in grayscale\n");
		}

		// Haralick's masks for computing k1 to k10
		double masks[10][25] = { {-13,   2,   7,  2, -13,   2,  17,  22,  17,   2,  7,  22,  27, 22,  7,  2,  17, 22, 17,  2, -13,   2,  7,  2, -13},
								 { 31,  -5, -17, -5,  31, -44, -62, -68, -62, -44,  0,   0,   0,  0,  0, 44,  62, 68, 62, 44, -31,   5, 17,  5, -31},
								 { 31, -44,   0, 44, -31,  -5, -62,   0,  62,   5, -17, -68,  0, 68, 17, -5, -62,  0, 62,  5,  31, -44,  0, 44, -31},
 								 {  2,   2,   2,  2,   2,  -1,  -1,  -1,  -1,  -1,  -2,  -2, -2, -2, -2, -1,  -1, -1, -1, -1,   2,   2,  2,  2,   2},
								 {  4,   2,   0, -2,  -4,   2,   1,   0,  -1,  -2,   0,   0,  0,  0,  0, -2,  -1,  0,  1,  2,  -4,  -2,  0,  2,   4},
								 {  2,  -1,  -2, -1,   2,   2,  -1,  -2,  -1,   2,   2,  -1, -2, -1,  2,  2,  -1, -2, -1,  2,   2,  -1, -2, -1,   2},
								 { -1,  -1,  -1, -1,  -1,   2,   2,   2,   2,   2,   0,   0,  0,  0,  0, -2,  -2, -2, -2, -2,   1,   1,  1,  1,   1},
								 { -4,  -2,   0,  2,   4,   2,   1,   0,  -1,  -2,   4,   2,  0, -2, -4,  2,   1,  0, -1, -2,  -4,  -2,  0,  2,   4},
								 { -4,   2,   4,  2,  -4,  -2,   1,   2,   1,  -2,   0,   0,  0,  0,  0,  2,  -1, -2, -1,  2,   4,  -2, -4, -2,   4},
								 { -1,   2,   0, -2,   1,  -1,   2,   0,  -2,   1,  -1,   2,  0, -2,  1, -1,   2,  0, -2,  1,  -1,   2,  0, -2,   1} };
		double weights[10] = {1/175, 1/420, 1/420, 1/70, 1/100, 1/70, 1/60, 1/140, 1/140, 1/60};

		// Initialise edge image
		float *edges = malloc(w*h*sizeof(float));

		// Zero-padding
		int wx = (w+8);
		int hx = (h+8);
		double *aux = calloc(wx*hx,sizeof(double));
		int i,j,fila,col;
		int imax = wx*hx;
		for(i=0;i<imax;i++){
			fila = (int)(i/wx);
			col = i-(wx*fila);	
			if ( (fila>=4)&&(col>=4)&&(fila<h+4)&&(col<w+4) ) {
				aux[i] = im[(col-4)+(w*(fila-4))];
			}
		}
		
		// Haralick's algorithm
		// TODO
		int i_zp;
		double k[10];
		for(fila=0;fila<h;fila++){
			for(col=0;col<w;col++){
				i = col + w*fila;				// original image & edges images index
				i_zp = (col+4) + wx*(fila+4);	// zero padded image index
				
				// k1 to k10 (note: k1 is not necessary)
			}
		}
		

		// Save output image
		iio_save_image_float_vec(argv[3], edges, w, h, 1);
		fprintf(stderr, "Output Image saved in %s:\t %dx%d image with %d channel(s).\n", argv[3], w, h, pixeldim);
	
		// Free memory
		free(im_orig);
		free(im);

		fprintf(stderr, "haralick's edge detector computation done.\n");

		// Execution time:
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		fprintf(stderr, "execution time: %1.3f s.\n", exectime);		

		return 0;
	
	} // else (argc)
}
