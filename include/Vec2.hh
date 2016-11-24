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

#ifndef VEC2_H
#define VEC2_H

class Vec2 {
    private:
        double comps[2];

    public:

        Vec2() {
            this->comps[0] = 0.0;
            this->comps[1] = 0.0;
        } 

        Vec2(double x, double y) {
            this->comps[0] = x;
            this->comps[1] = y;
        }

        /* copy constructor */
        Vec2(const Vec2 &x) {
            this->comps[0]  = x.comps[0];
            this->comps[1]  = x.comps[1];
        }

        ~Vec2() {}

        Vec2
        operator+(const Vec2 &v) const 
        {
            Vec2 s;
            s.comps[0] = this->comps[0] + v.comps[0];
            s.comps[1] = this->comps[1] + v.comps[1];
            return s;
        }

        Vec2
        operator-(const Vec2 &v) const 
        {
            Vec2 s;
            s.comps[0] = this->comps[0] - v.comps[0];
            s.comps[1] = this->comps[1] - v.comps[1];
            return s;
        }

        Vec2 &
        operator=(const Vec2 &a) {
            this->comps[0] = a.comps[0];
            this->comps[1] = a.comps[1];
            return (*this);
        }

        Vec2 &
        operator*=(const double &x)
        {
            this->comps[0] *= x;
            this->comps[1] *= x;

            return (*this);
        }

        Vec2 &
        operator/=(const double &x)
        {
            this->comps[0] /= x;
            this->comps[1] /= x;

            return (*this);
        }

        Vec2 &
        operator+=(const Vec2 &a)
        {
            this->comps[0] += a.comps[0];
            this->comps[1] += a.comps[1];

            return (*this);
        }

        Vec2 &
        operator-=(const Vec2 &a)
        {
            this->comps[0] -= a.comps[0];
            this->comps[1] -= a.comps[1];

            return (*this);
        }

        bool
        operator==(const Vec2 &a) 
        {
            /* CAUTION: we use '==' for double precision floating point numbers here, in most cases, this is MOST
             * UNWISE, since there's no well defined equality for those numbers. usually, we'd have to use (a - b) <
             * machine_eps. here, we use that to compare with null vector or custom set unit vectors... */
            return (this->comps[0] == a.comps[0] &&
                    this->comps[1] == a.comps[1]);
        }

        bool
        operator!=(const Vec2 &a)
        {
            return (this->comps[0] != a.comps[0] ||
                    this->comps[1] != a.comps[1]);
        }

        double &
        operator[](const uint32_t i)
        {
            if (i >= 2) {
                throw ("Vec2::operator[]: can't access member with index >= 3 in Vec2.\n");
            }
            return this->comps[i];
        }


        /* length for p2-norm (eucl) */
        double
        len2(void) const 
        {
            return sqrt(comps[0]*comps[0] + comps[1]*comps[1]);
        }

        /* len squares */
        double
        len2squared(void) const 
        {
            return (comps[0]*comps[0] + comps[1]*comps[1]);
        }

        void
        normalize(void) 
        {
            double len = this->len2();
            this->comps[0] /= len;
            this->comps[1] /= len;
        }

        void
        scale(double alpha) 
        {
            this->comps[0] *= alpha;
            this->comps[1] *= alpha;
        }

        /* scalar product, inner product */
        double
        operator*(const Vec2 &b) const
        {
            double scalprod = 0.0;

            for (unsigned i = 0; i < 2; i++) {
                scalprod += this->comps[i] * b.comps[i];
            }
            return scalprod;
        }

        /* multiply with c e |R */
        Vec2
        operator*(const double &c) const 
        {
            Vec2 result;

            for (unsigned i = 0; i < 2; i++) {
                result.comps[i] = c * this->comps[i];
            }
            return result;
        }

        /* "cross product" for 2d vectors */
        double
        cross(const Vec2 &b) const
        {
            return (comps[0]*b.comps[1] - comps[1]*b.comps[0]);
        }

        void
        set(double x, double y)
        {
            this->comps[0] = x;
            this->comps[1] = y;
        }

        void
        print()
        {
            printf("(%f, %f)\n", this->comps[0], this->comps[1]);
        }

};

#endif
