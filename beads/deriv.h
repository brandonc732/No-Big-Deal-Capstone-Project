/*
     author: John Hughes
      email: jph264@psu.edu
       date: August 18, 2008
description: This header implements numerical second-order partial
             differentiation, which we need to compute the observed Fisher
             information matrix.

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

#ifndef DERIV_H_
#define DERIV_H_

#include <gsl/gsl_math.h>
#include "types.h"

size_t beads;  // This global is accessed throughout the application.

// Function rate evaluates the Gaussian function centered at (x0, y0) and
// having standard deviation S and scaling parameter A.

double rate(double x, double y, double x0, double y0, double A, double S)
{
    return A * exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (S * S));
}

/*
This function evaluates the log-likelihood function at the point params.
The pixel data is passed to logL as auxiliary data. See our paper for more
information.
*/

double logL(const std::vector<double>& params, void * aux)
{
	Pixels * pixels = (Pixels *)aux;
	size_t n = pixels->getN();
	double * x = pixels->getX(),
		   * y = pixels->getY(),
		   * z = pixels->getZ(),
		   * x0 = new double[beads],
		   * y0 = new double[beads],
		   * A = new double[beads];
	size_t i,
		   j;
	for (i = 0, j = 0; i < 3 * beads; i += 3, j++)
	{
		x0[j] = params[i];
		y0[j] = params[i + 1];
		A[j] = params[i + 2];
	}
	double S = params[i],
		   B = params[i + 1],
	       theta = params[i + 2];
	double value = 0;
	for (size_t i = 0; i < n; i++)
	{
		double fi = B;
		for (size_t j = 0; j < beads; j++)
			fi += rate(x[i], y[i], x0[j], y0[j], A[j], S);
		value += log(fi + theta) + (z[i] - fi) * (z[i] - fi) / (fi + theta);
	}
	delete [] x0;
	delete [] y0;
	delete [] A;
	return -0.5 * (n * log(2 * M_PI) + value);
}

/*
The following function computes the second-order partial derivative of function
F with respect to parameters i and j. The central difference formula is used.
The last two parameters are an initial step size and a convergence criterion.
*/ 

double differentiate(Function * F, const std::vector<double>& params, size_t i, size_t j, double step, double tolerance)
{
	double (*f)(const std::vector<double>& params, void * aux) = F->f;
	size_t p = params.size();
	void * aux = F->aux;
	std::vector<double> args(p);
	double h = step,
		   last,
		   current,
		   first,
		   second,
		   third,
		   fourth;
	for (size_t k = 0; k < p; k++)
	     args[k] = params[k];
	if (i == j)
	{
        args[i] = params[i] + h;
        first = (*f)(args, aux);
        args[i] = params[i];
        second = (*f)(args, aux);
        args[i] = params[i] - h;
        third = (*f)(args, aux);
	    current = (first - 2 * second + third) / (h * h);
	}
	else
	{
        args[i] = params[i] + h;
        args[j] = params[j] + h;
        first = (*f)(args, aux);
        args[j] = params[j] - h;
        second = (*f)(args, aux);
        args[i] = params[i] - h;
        fourth = (*f)(args, aux);
        args[j] = params[j] + h;
        third = (*f)(args, aux);
        current = (first - second - third + fourth) / (4 * h * h);
	}
	do
	{
		last = current;
		h /= 2;
		if (i == j)
	    {
		    args[i] = params[i] + h;
		    first = (*f)(args, aux);
		    args[i] = params[i];
		    second = (*f)(args, aux);
		    args[i] = params[i] - h;
		    third = (*f)(args, aux);
			current = (first - 2 * second + third) / (h * h);
		}
		else
		{
			args[i] = params[i] + h;
			args[j] = params[j] + h;
			first = (*f)(args, aux);
			args[j] = params[j] - h;
			second = (*f)(args, aux);
			args[i] = params[i] - h;
			fourth = (*f)(args, aux);
			args[j] = params[j] + h;
			third = (*f)(args, aux);
			current = (first - second - third + fourth) / (4 * h * h);
		}
	}
	while (fabs(current - last) > tolerance);
	return current;
}

#endif /* DERIV_H_ */
