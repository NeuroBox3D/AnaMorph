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

#include "Vec3.hh"

template <typename R>
Vec3<R>::Vec3() : Vector<R>(3)
{
    this->fill(0.0);
}

template <typename R>
Vec3<R>::Vec3(const R& v)
: Vector<R>(3, v)
{}

template <typename R>
Vec3<R>::Vec3(
    R const &x,
    R const &y,
    R const &z) : Vector<R>(3)
{
    this->operator()(0) = x;
    this->operator()(1) = y;
    this->operator()(2) = z;
}

/* copy constructor / assignment operator */
template <typename R>
Vec3<R>::Vec3(Vec3 const &v)
: Vector<R>(v) {}
//{
//   Vector<R>::operator=(v);
//}

template <typename R>
Vec3<R> &
Vec3<R>::operator=(Vec3 const &v) {
    Vector<R>::operator=(v);
    return (*this);
}

/* compatibility size-based constructor and resize() that accepts size 3 only. useful if Vec3 is used as a
 * template argument that requires this type of constructor */
template <typename R>
Vec3<R>::Vec3(uint32_t n) : Vector<R>(3)
{
    if (n != 3) {
        throw("Vec3<R>::resize(size_t): given size != 3.");
    }
}

template <typename R>
Vec3<R>::~Vec3()
{
}

template <typename R>
void
Vec3<R>::resize(uint32_t size) 
{
    if (size != 3) {
        throw("Vec3<R>::resize(size_t): given size != 3.");
    }
}

/* component access */
template <typename R>
R &
Vec3<R>::operator()(uint32_t i)
{
    if (i >= 3) {
        throw ("Vec3<R>::operator[]: can't access member with index >= 3 in Vec3.\n");
    }
    return this->comps[i];
}

template <typename R>
R const &
Vec3<R>::operator()(uint32_t i) const
{
    if (i >= 3) {
        throw ("Vec3<R>::operator[]: can't access member with index >= 3 in Vec3.\n");
    }
    return this->comps[i];
}

template <typename R>
R &
Vec3<R>::operator[](uint32_t i)
{
    if (i >= 3) {
        throw ("Vec3<R>::operator[]: can't access member with index >= 3 in Vec3.\n");
    }
    return this->comps[i];
}

template <typename R>
R const &
Vec3<R>::operator[](uint32_t i) const
{
    if (i >= 3) {
        throw ("Vec3<R>::operator[]: can't access member with index >= 3 in Vec3.\n");
    }
    return this->comps[i];
}

/* arithmetic */
template <typename R>
Vec3<R>
Vec3<R>::operator+(Vec3 const &v) const 
{
    return (Vec3(*this) += v);
}

template <typename R>
Vec3<R> &
Vec3<R>::operator+=(Vec3 const &v)
{
    Vector<R>::operator+=(v);
    return (*this);
}

template <typename R>
Vec3<R>
Vec3<R>::operator-(Vec3 const &v) const 
{
    return (Vec3(*this) -= v);
}

template <typename R>
Vec3<R> &
Vec3<R>::operator-=(Vec3 const &v)
{
    Vector<R>::operator-=(v);
    return (*this);
}

template <typename R>
Vec3<R>
Vec3<R>::operator*(R const &x) const
{
    return (Vec3(*this) *= x);
}

template <typename R>
Vec3<R> &
Vec3<R>::operator*=(R const &x)
{
    Vector<R>::operator*=(x);
    return (*this);
}

template <typename R>
Vec3<R>
Vec3<R>::operator/(R const &x) const
{
    return (Vec3(*this) /= x);
}

template <typename R>
Vec3<R> &
Vec3<R>::operator/=(R const &x)
{
    Vector<R>::operator/=(x);
    return (*this);
}

/* scalar product, inner product */
template <typename R>
R
Vec3<R>::operator*(Vec3 const &v) const
{
    R scalprod = 0.0;

    for (uint32_t i = 0; i < 3; i++) {
        scalprod += this->operator[](i) * v(i);
    }
    return scalprod;
}

/* cross product */
template <typename R>
Vec3<R>
Vec3<R>::cross(Vec3 const &v) const
{
    Vec3 result;

    /*
    result.comps[0] = this->comps[1]*b.comps[2] - this->comps[2]*b.comps[1];
    result.comps[1] = this->comps[2]*b.comps[0] - this->comps[0]*b.comps[2];
    result.comps[2] = this->comps[0]*b.comps[1] - this->comps[1]*b.comps[0];
    */
    result[0] = this->operator()(1)*v(2) - this->operator()(2)*v(1);
    result[1] = this->operator()(2)*v(0) - this->operator()(0)*v(2);
    result[2] = this->operator()(0)*v(1) - this->operator()(1)*v(0);

    return result;
}

template <typename R>
bool
Vec3<R>::operator==(Vec3 const &v) const
{
    /* CAUTION: we use '==' for type R (usually R = float, double, long double) here, in most cases, this is MOST
     * UNWISE, since there's no well defined equality for those numbers. usually, we'd have to use (a - b) <
     * machine_eps. here, we use that to compare with null vector or custom set vectors which are supposed to match
     * bitwise with respect to IEEE 754 */
    return (this->operator()(0) == v(0) &&
            this->operator()(1) == v(1) &&
            this->operator()(2) == v(2)
        );
}

template <typename R>
bool
Vec3<R>::operator!=(Vec3 const &v)
{
    return !((*this) == v);
}

/* tuple-like comparison operators: a < b iff a[0] < b[0] && .. && a[2] < b[2].
 *
 * a > b iff !(a == b) && !(a < b) */
template <typename R>
bool 
Vec3<R>::operator<(Vec3 const &v) const {
    return (this->operator()(0) < v(0) &&
            this->operator()(1) < v(1) &&
            this->operator()(2) < v(2)
        );
}

template <typename R>
bool
Vec3<R>::operator>(Vec3 const &v) const{
    return !((*this) == v || (*this) < v);
}

template <typename R>
bool
Vec3<R>::operator>=(Vec3 const &v) const{
    return !((*this) < v);
}

template <typename R>
bool
Vec3<R>::operator<=(Vec3 const &v) const{
    return !((*this) > v);
}

/* length for p2-norm (eucl) */
template <typename R>
R
Vec3<R>::len2(void) const 
{
    return std::sqrt(this->len2squared());
}

/* squared length for p2-norm */
template <typename R>
R
Vec3<R>::len2squared(void) const 
{
    return (
            this->operator()(0) * this->operator()(0) +
            this->operator()(1) * this->operator()(1) +
            this->operator()(2) * this->operator()(2)
        );
}

template <typename R>
Vec3<R> &
Vec3<R>::normalize(void) 
{
    R len = this->len2();
    this->comps[0] /= len;
    this->comps[1] /= len;
    this->comps[2] /= len;

    return (*this);
}

template <typename R>
void
Vec3<R>::print() const
{
    printf("(%f, %f, %f)\n", this->comps[0], this->comps[1], this->comps[2]);
}

template <typename R>
void
Vec3<R>::print_debugl(uint32_t level) const
{
    debugTabInc();
    debugl(level, "(%f, %f, %f)\n", this->comps[0], this->comps[1], this->comps[2]);
    debugTabDec();
}


