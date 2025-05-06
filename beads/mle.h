/*
     author: John Hughes
      email: jph264@psu.edu
       date: August 18, 2008
description: This header implements parameter estimation and related tasks.
             It makes heavy use of the GNU Scientific Library.

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

#ifndef MLE_H_
#define MLE_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "deriv.h"

//std::ofstream OF("simasse.txt");

#define FIT(i) gsl_vector_get(solver->x, i)
#define ERR(i) sqrt(gsl_matrix_get(V, i, i))
#define VAR(i) gsl_matrix_get(V, i, i)
#define COV(i, j) gsl_matrix_get(V, i, j)
#define GET(i) gsl_vector_get(params, i)

// We use these values as initial estimates of S, B, A_j, and theta.

const double S0 = 200,
             B0 = 200,
             A0 = 10000,
             theta0 = 10;

/*
The following four functions implement OLS estimation of beta \ {theta}.
See the GSL documentation for more information about the GSL least-squares
framework.
*/

int OLS_f(const gsl_vector * params, void * data, gsl_vector * f)
{
    Pixels * pixels = (Pixels *)data;
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
        x0[j] = GET(i);
        y0[j] = GET(i + 1);
        A[j] = GET(i + 2);
    }
    double S = GET(i),
           B = GET(i + 1);
    for (i = 0; i < n; i++)
    {
        double Zi = B;
        for (j = 0; j < beads; j++)
            Zi += rate(x[i], y[i], x0[j], y0[j], A[j], S);
        gsl_vector_set(f, i, Zi - z[i]);
    }
    delete [] x0;
    delete [] y0;
    delete [] A;
    return GSL_SUCCESS;
}

int OLS_df(const gsl_vector * params, void * data, gsl_matrix * J)
{
    Pixels * pixels = (Pixels *)data;
    size_t n = pixels->getN();
    double * x = pixels->getX(),
           * y = pixels->getY(),
           * x0 = new double[beads],
           * y0 = new double[beads],
           * A = new double[beads];
    size_t i,
           j;
    for (i = 0, j = 0; i < 3 * beads; i += 3, j++)
    {
        x0[j] = GET(i);
        y0[j] = GET(i + 1);
        A[j] = GET(i + 2);
    }
    double S = GET(i),
           B = GET(i + 1);
    double * rates = new double[beads];
    for (i = 0; i < n; i++)
    {
        double f = B;
        double Snumerator = 0;
        size_t j;
        for (j = 0; j < beads; j++)
        {
            rates[j] = rate(x[i], y[i], x0[j], y0[j], A[j], S);
            f += rates[j];
            Snumerator += ((x[i] - x0[j]) * (x[i] - x0[j]) + (y[i] - y0[j]) * (y[i] - y0[j])) * rates[j];
        }
        size_t k;
        for (j = 0, k = 0; j < 3 * beads; j += 3, k++)
        {
            gsl_matrix_set(J, i, j, 2 * (x[i] - x0[k]) / (S * S) * rates[k]);
            gsl_matrix_set(J, i, j + 1, 2 * (y[i] - y0[k]) / (S * S) * rates[k]);
            gsl_matrix_set(J, i, j + 2, rates[k] / A[k]);
        }
        gsl_matrix_set(J, i, j, 2 * Snumerator / (S * S * S));
        gsl_matrix_set(J, i, j + 1, 1);
    }
    delete [] x0;
    delete [] y0;
    delete [] A;
    delete [] rates;
    return GSL_SUCCESS;
}

int OLS_fdf(const gsl_vector * params, void * data, gsl_vector * f, gsl_matrix * J)
{
    OLS_f(params, data, f);
    OLS_df(params, data, J);
    return GSL_SUCCESS;
}

/*
The following function is the driver for our OLS fitting. The function returns
the approximate information criterion described in our paper.
*/

double OLS_fit(void * data, size_t p, std::vector<double>& beta)
{
    Pixels * pixels = (Pixels *)data;
    size_t n = pixels->getN();
    const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver * solver = gsl_multifit_fdfsolver_alloc(T, n, p);
    gsl_multifit_function_fdf f;
    f.f = &OLS_f;
    f.df = &OLS_df;
    f.fdf = &OLS_fdf;
    f.n = n;
    f.p = p;
    f.params = data;
    double * init = new double[p];
    size_t j;
    for (j = 0; j < p; j++)
    	init[j] = beta[j];
    gsl_vector_view params = gsl_vector_view_array(init, p);
    gsl_multifit_fdfsolver_set(solver, &f, &params.vector);
    int status;
    size_t iteration = 0;
    do
    {
        iteration++;
        gsl_multifit_fdfsolver_iterate(solver);
        status = gsl_multifit_test_delta(solver->dx, solver->x, 0.0001, 0.0001);
    }
    while (status == GSL_CONTINUE && iteration < 500);
    for (j = 0; j < p; j++)
        beta[j] = FIT(j);
    double RSS = gsl_blas_dnrm2(solver->f);
    RSS *= RSS;
    double IC = p * sqrt(n) + n * log(RSS / n);
    gsl_multifit_fdfsolver_free(solver);
    delete [] init;
    return IC;
}

/*
The next two functions implement our maximum likelihood estimation. We maximize
the likelihood by minimizing the negative log-likelihood using the Nelder-Mead
simplex algorithm. See the GSL documentation for more information regarding the
GSL minimization framework.
*/

double MLE_f(const gsl_vector * params, void * data)
{
    Pixels * pixels = (Pixels *)data;
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
        x0[j] = GET(i);
        y0[j] = GET(i + 1);
        A[j] = GET(i + 2);
    }
    double S = GET(i),
           B = GET(i + 1),
           theta = GET(i + 2);
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
    return value;
}

// This function is the driver for our MLE. The function returns the value of
// the log-likelihood at the (perhaps local) maximum.

double estimate(void * data, size_t p, std::vector<double>& beta)
{
    const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer * solver = gsl_multimin_fminimizer_alloc(T, p);
    gsl_multimin_function f;
    f.n = p;
    f.f = &MLE_f;
    f.params = data;
    gsl_vector * init = gsl_vector_alloc(p);
    size_t i,
           j;
    for (j = 0; j < p; j++)
        gsl_vector_set(init, j, beta[j]);
    gsl_vector * ss = gsl_vector_alloc(p);
    gsl_vector_set_all(ss, 1);
    gsl_multimin_fminimizer_set(solver, &f, init, ss);
    int status;
    size_t iteration = 0;
    double size;
    do
    {
        iteration++;
        status = gsl_multimin_fminimizer_iterate(solver);
        if (status)
        {
        	std::cout << "error\n";
        	break;
        }
        size = gsl_multimin_fminimizer_size(solver);
        status = gsl_multimin_test_size(size, 0.1);
    }
    while (status == GSL_CONTINUE && iteration < 2500);
    //std::cout << "iterations = " << iteration << std::endl;
    for (i = 0, j = 0; i < 3 * beads; i += 3, j++)
	{
		beta[i] = FIT(i);
		beta[i + 1] = FIT(i + 1);
		beta[i + 2] = FIT(i + 2);
	}
	beta[i] = fabs(FIT(i));
	beta[i + 1] = FIT(i + 1);
	beta[i + 2] = FIT(i + 2);
    gsl_vector_free(init);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(solver);
    size_t n = ((Pixels *)data)->getN();
    double logL = -0.5 * (n * log(2 * M_PI) + solver->fval);
    return logL;
}

// MLE_cov computes and returns the approximate covariance matrix for the
// estimators.

gsl_matrix * MLE_cov(Pixels * pixels, const std::vector<double>& beta)
{
    Function F;
    F.f = &logL;
    F.aux = pixels;
    double df;
	size_t p = beta.size();
    gsl_matrix * J = gsl_matrix_alloc(p, p);

    // First compute the observed Fisher information matrix.

    for (size_t i = 0; i < p; i++)
        for (size_t j = i; j < p; j++)
    	{
        	df = differentiate(&F, beta, i, j, 0.1, 1e-3);
        	gsl_matrix_set(J, i, j, -df);
        	gsl_matrix_set(J, j, i, -df);
        }
    gsl_permutation * perm = gsl_permutation_alloc(p);
    int signum;

    // Then LU decompose the matrix.

    gsl_linalg_LU_decomp(J, perm, &signum);
    gsl_matrix * V = gsl_matrix_alloc(p, p);

    // Finally, invert the matrix.

    gsl_linalg_LU_invert(J, perm, V);
    gsl_matrix_free(J);
    gsl_permutation_free(perm);
    return V;
}

// This function prints the parameter estimates and their approximate standard
// errors to the supplied output stream.

void printResults(const std::vector<double>& beta, Pixels * pixels, std::ostream& out)
{
	gsl_matrix * V = MLE_cov(pixels, beta);
	size_t i,
	       j;

	/*for (i = 0; i < p; i++)
	{
		for (j = 0; j < p; j++)
			out << gsl_matrix_get(cov, i, j) << " ";
		out << std::endl;
	}*/

	for (i = 0, j = 0; i < 3 * beads; i += 3, j++)
	{
		out << "    x" << j << " = " << beta[i] << " +/- " << ERR(i) << "\n"
			<< "    y" << j << " = " << beta[i + 1] << " +/- " << ERR(i + 1) << "\n"
			<< "    A" << j << " = " << beta[i + 2] << " +/- " << ERR(i + 2) << "\n";
	}
	out << "     S = " << beta[i] << " +/- " << ERR(i) << "\n"
		<< "     B = " << beta[i + 1] << " +/- " << ERR(i + 1) << "\n"
		<< " theta = " << beta[i + 2] << " +/- " << ERR(i + 2) << std::endl;

	/*for (i = 0; i < 3 * beads; i += 3)
		out << beta[i + 2] << " " << VAR(i + 2) << " ";
	out << beta[i] << " " << VAR(i) << " " << beta[i + 1] << " " << VAR(i + 1) << " " << beta[i + 2] << " " << VAR(i + 2) << std::endl;
    */
	/*for (i = 0, j = 0; i < 3 * beads; i += 3, j++)
	    OF << beta[i] << " " << beta[i + 1] << " " << beta[i + 2] << " ";
    OF << beta[i] << " " << beta[i + 1] << " " << beta[i + 2] << std::endl;*/

    gsl_matrix_free(V);
}

/*
This function takes a two-dimensional matrix of pixel values and returns the
location of the center of the pixel with the largest intensity. The location
is given in nanometer offsets from the image's top left corner.
*/

void findMax(size_t& x0, size_t& y0, double ** map, size_t dim, size_t width)
{
    double maxz = 0;
    size_t i;
    int j,
        maxi = 0,
        maxj = 0;
    for (i = 0; i < dim; i++)
        for (j = 0; j < int(dim); j++)
        {
            if (map[i][j] > maxz)
            {
                maxz = map[i][j];
                maxi = i;
                maxj = j;
            }
        }
    x0 = static_cast<size_t>(maxj * width + width / 2.0);
    y0 = static_cast<size_t>(maxi * width + width / 2.0);
    maxi -= 3;
    maxj -= 3;
    for (i = 0; i < 7; i++)
    {
        for (j = 0; j < 7; j++)
        	if (maxi >= 0 && maxi < int(dim) && maxj + j > 0 && maxj + j < int(dim))
                map[maxi][maxj + j] = 0;
        maxi++;
    }
}

/*
findBeads builds a two-dimensional array of pixel values and uses repeated calls
to findMax to fill beta with initial estimates of bead locations. The global
variable beads determines how many calls to findMax should be made.
*/

void findBeads(Pixels * pixels, size_t width, std::vector<double>& beta)
{
    size_t dim = static_cast<size_t>(sqrt(pixels->getN()));
    double ** map = new double *[dim];
    size_t j;
    for (j = 0; j < dim; j++)
        map[j] = new double[dim];
    size_t i,
           k = 0;
    double * z = pixels->getZ();
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            map[i][j] = z[k++];      //*******************
    size_t x0,
           y0;
    for (j = 0, k = 0; j < beads; j++, k += 3)
    {
        findMax(x0, y0, map, dim, width);
        beta[k] = x0;
        beta[k + 1] = y0;
        beta[k + 2] = A0;
    }
    for (j = 0; j < dim; j++)
        delete [] map[j];
    delete [] map;
}

#endif /*MLE_H_*/
