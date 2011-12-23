
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
		double weights[10] = {175, 420, 420, 70, 100, 70, 60, 140, 140, 60};
		int i,j;
		for(i=0;i<10;i++){
			for(j=0;j<25;j++){
				masks[i][j] /= weights[i];
			}
		}

		// Initialise edge image
		float *edges = calloc(w*h,sizeof(float));

		// Zero-padding
		int wx = (w+8);
		int hx = (h+8);
		double *aux = calloc(wx*hx,sizeof(double));
		int fila,col;
		int imax = wx*hx;
		for(i=0;i<imax;i++){
			fila = (int)(i/wx);
			col = i-(wx*fila);	
			if ( (fila>=4)&&(col>=4)&&(fila<h+4)&&(col<w+4) ) {
				aux[i] = im[(col-4)+(w*(fila-4))];
			}
		}
		
		// Haralick's algorithm
		int i_zp, u, v, num_edges;
		num_edges = 0;
		double k[10];
		int *offsets = get_neighbors_offset(wx, 5);
		double acum;
		double C2, C3, denom, sintheta, costheta;
		for(fila=0;fila<h;fila++){
			for(col=0;col<w;col++){
				i = col + w*fila;				// original image & edges images index
				i_zp = (col+4) + wx*(fila+4);	// zero padded image index
				
				// need neighborhood
				double *neighborhood = get_neighborhood(aux, i_zp, 5, offsets);

				// k1 to k10 (note: k1 (u=0) is not necessary)
				for(u=1;u<10;u++){
					acum = 0;
					for(v=0;v<25;v++){
						acum += neighborhood[v]*masks[u][v];
					}
					k[u] = acum;
				}
			
				// compute C2 and C3
				denom = k[1]*k[1] + k[2]*k[2];
				C2 = ( k[1]*k[1]*k[3] + k[1]*k[2]*k[4] + k[2]*k[2]*k[5] ) / denom;
				C3 = ( k[1]*k[1]*k[1]*k[6] + k[1]*k[1]*k[2]*k[7] + k[1]*k[2]*k[2]*k[8] + k[2]*k[2]*k[2]*k[9] ) / ( denom*sqrt(denom) );

				denom = fabs(C2 / (3*C3));	// reuse of denom...
				// edge pixel conditions
				if ( (C3<0) && (denom<rhozero) ) {
					edges[i] = 255;
					num_edges += 1;
				}

				// free neighborhood
				free_neighborhood(neighborhood);
			}
		}
		
		// Result
		fprintf(stderr, "%d edge points found...\n", num_edges);

		// Save output image
		iio_save_image_float_vec(argv[3], edges, w, h, 1);
		fprintf(stderr, "Output Image saved in %s:\t %dx%d image with %d channel(s).\n", argv[3], w, h, pixeldim);
	
		// Free memory
		free(im_orig);
		free(im);
		free(aux);
		free(edges);

		fprintf(stderr, "haralick's edge detector computation done.\n");

		// Execution time:
		double finish = (double)clock();
		double exectime = (finish - start)/CLOCKS_PER_SEC;
		fprintf(stderr, "\texecution time: %1.3f s.\n", exectime);		

		return 0;
	
	} // else (argc)
}
