/* --------------------------------------------------------------------------------
 * 
 *                              THIS FILE IS PART OF                               
 * 
 * AnaMorph: A Framework for Geometric Modelling, Consistency Analysis and Surface
 * Mesh Generation of Anatomically Reconstructed Neuron Morphologies
 * 
 * Web site: http://www.anamorph.net
 * 
 * Copyright (c) 2013-2014, Konstantin Mörschel.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 * 
 *    This product includes software developed by Konstantin Mörschel for
 *    AnaMorph (http://www.anamorph.net).
 * 
 * 4. Neither the name "AnaMorph" nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without
 *    specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS OF ANAMORPH ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS OF ANAMORPH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * -------------------------------------------------------------------------------- */

#include "Polynomial.hh"
#include "BivariatePolynomial.hh"

/* print coefficients: specializations for R = F = {float, double} */
template<>
void
Polynomial<float, float>::printCoeff() const
{
    uint32_t        i;
    static char     s[1024];
    static char     tmp[256];

    sprintf(s, "[ ");
    for (i = 0; i < degree + 1; i++) {
        sprintf(tmp, "%+20.13E ", coeff[i]);
        strncat(s, tmp, 256);
    }
    strncat(s, "]\n", 2);
    debugl(0, "%s", s);
}

template<>
void
Polynomial<double, double>::printCoeff() const
{
    uint32_t        i;
    static char     s[1024];
    static char     tmp[256];

    sprintf(s, "[ ");
    for (i = 0; i < degree + 1; i++) {
        sprintf(tmp, "%+20.13E ", coeff[i]);
        strncat(s, tmp, 256);
    }
    strncat(s, "]\n", 2);
    debugl(0, "%s", s);
    printf("%s", s);
}

/* write plot file suitable for gnuplot: specializations for float and double */
template<>
void
Polynomial<float, float>::writePlotFile(
    float           t0,
    float           t1,
    uint32_t        ticks,
    std::string     filename) const
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout) {
        fprintf(stderr, "Polynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str() );
        return;;
    }

    uint32_t    i;
    float       t;

    for (i = 0; i <= ticks; i++) {
        t = t0 + ((float)i / (float)ticks)*(t1 - t0);
        fprintf(fout, "%12.5E %12.5E\n", t, this->eval(t));
    }
    fclose(fout);
}

template<>
void
Polynomial<double, double>::writePlotFile(
    double          t0,
    double          t1,
    uint32_t        ticks,
    std::string     filename) const
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout) {
        fprintf(stderr, "Polynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str() );
        return;;
    }

    uint32_t    i;
    double      t;

    for (i = 0; i <= ticks; i++) {
        t = t0 + ((double)i / (double)ticks)*(t1 - t0);
        fprintf(fout, "%12.5E %12.5E\n", t, this->eval(t));
    }
    fclose(fout);
}

/* print coefficients: specializations for R = F = {float, double} */
template <>
void
BivariatePolynomial<float, float>::printCoeff() const
{
    uint32_t i, j;
    for (i = 0; i < m + 1; i++) {
        printf("( ");
        for (j = 0; j < n + 1; j++) {
            printf("%+20.13E ", coeff(i, j));
        }
        printf(")\n");
    }
}

template <>
void
BivariatePolynomial<double, double>::printCoeff() const
{
    uint32_t i, j;
    for (i = 0; i < m + 1; i++) {
        printf("( ");
        for (j = 0; j < n + 1; j++) {
            printf("%+20.13E ", coeff(i, j));
        }
        printf(")\n");
    }
}

/* write plot file suitable for gnuplot: specializations for float and double */
template <>
void
BivariatePolynomial<float, float>::writePlotFile(uint32_t ticks, std::string filename) const
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout) {
        fprintf(stderr, "BivariatePolynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str() );
        return;;
    }

    uint32_t    i, j;
    float       x, y;

    for (i = 0; i <= ticks; i++) {
        x = (float)i / (float)ticks;
        for (j = 0; j <= ticks; j++) {
            y = (float)j / (float)ticks;
            fprintf(fout, "%12.5E %12.5E %12.5E\n", x, y, this->eval(x, y) );
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}

template <>
void
BivariatePolynomial<double, double>::writePlotFile(uint32_t ticks, std::string filename) const
{
    FILE *fout = fopen(filename.c_str(), "w");
    if (!fout) {
        fprintf(stderr, "BivariatePolynomial::writePlotFile(): can't open file \"%s\" for writing.\n", filename.c_str() );
        return;;
    }

    uint32_t    i, j;
    double      x, y;

    for (i = 0; i <= ticks; i++) {
        x = (double)i / (double)ticks;
        for (j = 0; j <= ticks; j++) {
            y = (double)j / (double)ticks;
            fprintf(fout, "%12.5E %12.5E %12.5E\n", x, y, this->eval(x, y) );
        }
        fprintf(fout, "\n");
    }
    fclose(fout);
}
