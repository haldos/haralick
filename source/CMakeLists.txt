SET(LIBS iio 2dconvolution)
ADD_LIBRARY(iio SHARED iio.c iio.h)
ADD_LIBRARY(2dconv SHARED 2dconvolution.c 2dconvolution.h)
