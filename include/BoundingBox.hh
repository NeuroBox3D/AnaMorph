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

#ifndef BOUNDING_BOX_HH
#define BOUNDING_BOX_HH

#include "aux.hh"

template <typename R>
class BoundingBox {
    private:
        Vec3<R>             coord_min;
        Vec3<R>             coord_max;

    public:
        BoundingBox()
        {
            using Aux::Numbers::inf;
            using Aux::VecMat::onesVec3;

            this->coord_min = onesVec3<R>() * inf<R>();
            this->coord_max = onesVec3<R>() * (-inf<R>());
        }

        BoundingBox(const Vec3<R>& v1, const Vec3<R>& v2)
        {
            using Aux::VecMat::minVec3;
            using Aux::VecMat::maxVec3;

            minVec3<R>(coord_min, v1, v2);
            maxVec3<R>(coord_max, v1, v2);
        }

        BoundingBox(BoundingBox const &x)
        {
            this->coord_min = x.coord_min;
            this->coord_max = x.coord_max;
        }
    
        BoundingBox &
        operator=(BoundingBox const &x)
        {
            this->coord_min = x.coord_min;
            this->coord_max = x.coord_max;

            return (*this);
        }

       ~BoundingBox()
       {
       }

        bool
        operator<(BoundingBox const &x) const
        {
            return (
                this->coord_min[0] > x.coord_min[0] &&
                this->coord_min[1] > x.coord_min[1] &&
                this->coord_min[2] > x.coord_min[2] &&
                this->coord_max[0] < x.coord_max[0] &&
                this->coord_max[1] < x.coord_max[1] &&
                this->coord_max[2] < x.coord_max[2]);
        }

        bool
        operator<=(BoundingBox const &x) const
        {
            return (
                this->coord_min[0] >= x.coord_min[0] &&
                this->coord_min[1] >= x.coord_min[1] &&
                this->coord_min[2] >= x.coord_min[2] &&
                this->coord_max[0] <= x.coord_max[0] &&
                this->coord_max[1] <= x.coord_max[1] &&
                this->coord_max[2] <= x.coord_max[2]);
        }

        bool
        operator>(BoundingBox const &x) const
        {
            return (x < (*this));
        }

        bool
        operator>=(BoundingBox const &x) const
        {
            return (x <= (*this));
        }

        bool
        operator&&(BoundingBox const &x) const {
            return !(
                (this->coord_min[0] > x.coord_max[0] || this->coord_max[0] < x.coord_min[0]) ||
                (this->coord_min[1] > x.coord_max[1] || this->coord_max[1] < x.coord_min[1]) ||
                (this->coord_min[2] > x.coord_max[2] || this->coord_max[2] < x.coord_min[2])
            );
        }

        BoundingBox &
        extend(
            R       percentage,
            Vec3<R> offset_lb = { 0, 0, 0})
        {
            using Aux::VecMat::maxVec3;
            using Aux::VecMat::fabsVec3;

            Vec3<R> offset;
            maxVec3<R>(offset, offset_lb, fabsVec3<R>(this->coord_max - this->coord_min) * percentage);

            this->coord_min    -= offset;
            this->coord_max    += offset;

            return (*this);
        }

        BoundingBox &
        update(std::list<Vec3<R>> const &l)
        {
            for (auto &v : l) {
                Aux::VecMat::minVec3<R>(coord_min, coord_min, v);
                Aux::VecMat::maxVec3<R>(coord_max, coord_max, v);
            }

            return (*this);
        }

        BoundingBox &
        update(BoundingBox const &x) 
        {
            Aux::VecMat::minVec3<R>(coord_min, coord_min, x.min());
            Aux::VecMat::maxVec3<R>(coord_max, coord_max, x.max());
            return (*this);
        }

        Vec3<R>
        min() const
        {
            return (this->coord_min);
        }
        
        Vec3<R>
        max() const
        {
            return (this->coord_max);
        }
};

#endif
