/*
     author: John Hughes
      email: jph264@psu.edu
       date: August 18, 2008
description: This header implements a couple of types. The first is class
             Pixels, which is used throughout the application. The second is
             struct Function, which is used by the numerical differentiation
             routines.

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

#ifndef TYPES_H_
#define TYPES_H_

#include <vector>

class Pixels
{
    private:

        size_t _n;     // This is the number of observations (pixels).
        double * _x,   // (_x, _y) marks the center of a pixel, and _z is the
               * _y,   // pixel's intensity.
               * _z;

    public:

        Pixels(size_t n) : _n(n),
                           _x(new double[n]),
                           _y(new double[n]),
                           _z(new double[n]) {}
        ~Pixels()
        {
            delete [] _x;
            delete [] _y;
            delete [] _z;
        }
        size_t getN() const { return _n; }
        double * getX() const { return _x; }
        double * getY() const { return _y; }
        double * getZ() const { return _z; }
};

/*
A struct of type Function is used to specify a function to be differentiated.
The first field is a pointer to said function, and the second field supplies
any auxiliary data to the function.
*/

struct Function
{
	double (*f)(const std::vector<double>& params, void * aux);
	void * aux;
};

#endif /* TYPES_H_ */
