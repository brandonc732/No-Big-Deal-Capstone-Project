/*
     author: John Hughes
      email: jph264@psu.edu
       date: August 18, 2008
description: This header contains code to process a TIFF file and populate an
             instance of class Pixels with the data from the file.

Copyright (C) 2008 John Hughes

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef TIFF_H_
#define TIFF_H_

#include <tiffio.h>  // This header is part of libtiff.
#include "types.h"
//#include <cmath>

int w;

// This function's second parameter is the width of a pixel, in nanometers.

Pixels * processTIFF(TIFF * image, size_t pixelWidth)
{
	/*
	We place some restrictions on a TIFF file.

	1) It must use 16 bits per sample.
	2) It must define its width.
	3) We assume that the image is square.
	*/

    uint8 bps;
    if ((TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &bps) == 0) || (bps != 16))
    {
        std::cerr << "Your TIFF image must have 16 bits per sample.\n";
        return NULL;
    }
    uint32 width;
    if(TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &width) == 0)
    {
        std::cerr << "Your TIFF image does not define its width.\n";
        return NULL;
    }
    const tsize_t STRIP_SIZE = TIFFStripSize(image);
    const size_t STRIPS = TIFFNumberOfStrips(image),
                 BUFFER_SIZE = TIFFNumberOfStrips(image) * STRIP_SIZE;
    char * buffer;
    if ((buffer = new char[BUFFER_SIZE]) == NULL)
    {
        std::cerr << "I could not allocate enough memory to store the uncompressed image.\n";
        return NULL;
    }
    size_t imageOffset = 0;
    int result;
    for (size_t j = 0; j < STRIPS; j++)
    {
        if ((result = TIFFReadEncodedStrip(image, j, buffer + imageOffset, STRIP_SIZE)) == -1)
        {
            std::cerr << "Strip number " << j << "could not be read.\n";
            return NULL;
        }
        imageOffset += result;
    }
    size_t n = width * width;
    Pixels * pixels = new Pixels(n);
    unsigned short * ptr = (unsigned short *)buffer;        //*************************
    size_t k = 0;
    double * x = pixels->getX(),
           * y = pixels->getY(),
           * z = pixels->getZ();
    for (size_t i = 0; i < width; i++)
        for (size_t j = 0; j < width; j++)
        {
            z[k] = ptr[k];
            x[k] = j * w + w / 2.0;
            y[k++] = i * w + w / 2.0;
            //std::cout << x[k - 1] << " " << y[k - 1] << " " << z[k - 1] << "\n";
        }
    delete [] buffer;
    return pixels;
}

#endif /*TIFF_H_*/
