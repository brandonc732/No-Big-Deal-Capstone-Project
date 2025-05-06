/*
     author: John Hughes
      email: jph264@psu.edu
       date: August 18, 2008
description: This is the application's entry point.

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

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "tiff.h"
#include "mle.h"

#include <chrono>

using namespace std;

// This string shows how the program must be invoked.

const string USAGE = "usage: beads width file [dumpfile]\n";

void dumpToFile(Pixels * pixels, const vector<double>& beta, const string& filename);

int main(int argc, char * argv[])
{
    if (argc < 3)
    {
        cout << USAGE;
        return -1;
    }
    w = atoi(argv[1]);  // How wide is a pixel?
    if (w <= 0)
    {
        cout << "You must supply a valid pixel width, in nanometers.\n\n" << USAGE;
        return -1;
    }
    bool dump = false;
    if (argc > 3)
        dump = true;  // Should we dump the pixel data, fitted values, and
                      // residuals to a file?

    TIFFSetWarningHandler(NULL);
    TIFF * image = TIFFOpen(argv[2], "r");  // Try to open the image file.
    if (image == NULL)
    {
        cout << argv[2] << " could not be opened.\n";
        return -1;
    }
    Pixels * pixels = processTIFF(image, w);  // Extract the pixel data.
    TIFFClose(image);
    if (pixels == NULL)
    {
        cout << "Exiting due to undetermined error.\n";
        return -1;
    }
    beads = 0;
    size_t p = 3;
    vector<vector<double> > beta(1);    // Store vectors of estimates here.
    beta[0].resize(3);
    findBeads(pixels, w, beta[0]);
    beta[0][0] = S0;
    beta[0][1] = B0;
    beta[0][2] = theta0;

    /*
    Now carry out the preliminary phase of our algorithm. This estimates the
    number of beads and stores estimates of all parameters except theta. See the
    paper for details.
    */

    double IC = OLS_fit((void *)pixels, p - 1, beta[0]),
           oldIC;
    //for (size_t j = 0; j < beta[beads].size(); j++)
    //    cout << beta[beads][j] << endl;
    do
    {
    	//cout << "beads = " << beads << endl
    	//     << "   IC = " <<  IC << endl;

        cout << "trying " << beads + 1 << "beads:" << endl;

        auto start = std::chrono::high_resolution_clock::now();

    	beads++;
    	beta.resize(beads + 1);
    	oldIC = IC;
    	p += 3;
    	beta[beads].resize(p);
    	findBeads(pixels, w, beta[beads]);
    	beta[beads][p - 3] = S0;
    	beta[beads][p - 2] = B0;
    	beta[beads][p - 1] = theta0;
    	IC = OLS_fit((void *)pixels, p - 1, beta[beads]);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> duration = end - start;

        cout << "solved in " << duration.count() << " seconds, IC = " << IC << endl;

    	//for (size_t j = 0; j < beta[beads].size(); j++)
    	    //cout << beta[beads][j] << endl;
    }
    while (IC < oldIC);
    
    cout << "final beads = " << beads - 1 << endl
         << "   IC = " <<  IC << endl;
    
    //cout << "beads = " << beads << endl
         //<< "   IC = " <<  IC << endl;

    /*
    Now carry out the final stage of our algorithm. This refines the previous
    parameter estimates, estimates theta, and uses likelihood ratio tests to
    determine the number of beads.
    */

    /*double logLr = estimate((void *)pixels, p, beta[beads]),
           logLf,
           G2;
    for (size_t j = 0; j < beta[beads].size(); j++)
        cout << beta[beads][j] << endl;
    cout << "    p = " << p << endl
         << "beads = " << beads << endl
         << "logLf = " <<  logLf << endl
         << "logLr = " << logLr << endl
         << "  G^2 = " << G2 << endl;
    do
    {
    	beads--;
    	p -= 3;
    	logLf = logLr;
    	logLr = estimate((void *)pixels, p, beta[beads]);
    	for (size_t j = 0; j < beta[beads].size(); j++)
    		cout << beta[beads][j] << endl;
    	G2 = -2 * (logLr - logLf);
    	cout << "    p = " << p << endl
    	     << "beads = " << beads << endl
    	     << "logLf = " <<  logLf << endl
    	     << "logLr = " << logLr << endl
    	     << "  G^2 = " << G2 << endl;
    	if (beads == 0)
    	    break;
    }
    while (G2 < 7.8);
    beads++;
    p += 3;*/
    beads--;
    p -= 3;

    cout << "printing results" << beads << endl;

    printResults(beta[beads], pixels, cout);  // Send the results to STDOUT.

    cout << "dump=" << dump << endl;
    
	if (dump)
		dumpToFile(pixels, beta[beads], argv[3]);  // Dump, if applicable.
    delete pixels;
    return 0;
}

/*
This function dumps the pixel data, the fitted values, the raw residuals, and
the standardized residuals to a file. The file can be read directly into an R
data frame. Our R package, FIONAdiag, can use this data to provide various
diagnostics of the fit.
*/

void dumpToFile(Pixels * pixels, const vector<double>& beta, const string& filename)
{
	ofstream OF(filename.c_str());
	if (! OF.is_open())
		return;
	OF << "x y z fitted raw_res std_res" << endl;
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
	    x0[j] = beta[i];
	    y0[j] = beta[i + 1];
	    A[j] = beta[i + 2];
	}
	double S = beta[i],
	       B = beta[i + 1],
	       theta = beta[i + 2];
	for (size_t i = 0; i < n; i++)
	{
	    double fit = B;
	    for (size_t j = 0; j < beads; j++)
	        fit += rate(x[i], y[i], x0[j], y0[j], A[j], S);
	    double res = z[i] - fit;
	    OF << x[i] << " " << y[i] << " " << z[i] << " " << fit << " " << res << " " << res / sqrt(fit + theta) << endl;
	}
	OF.close();
	delete [] x0;
	delete [] y0;
	delete [] A;
}

