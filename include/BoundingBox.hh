/*
 * This file is part of
 *
 * AnaMorph: a framework for geometric modelling, consistency analysis and surface
 * mesh generation of anatomically reconstructed neuron morphologies.
 * 
 * Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group
 * Author: Konstantin Mörschel
 * 
 * AnaMorph is free software: Redistribution and use in source and binary forms,
 * with or without modification, are permitted under the terms of the
 * GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 *
 * (3) Neither the name "AnaMorph" nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * (4) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Mörschel K, Breit M, Queisser G. Generating neuron geometries for detailed
 *   three-dimensional simulations using AnaMorph. Neuroinformatics (2017)"
 * "Grein S, Stepniewski M, Reiter S, Knodel MM, Queisser G.
 *   1D-3D hybrid modelling – from multi-compartment models to full resolution
 *   models in space and time. Frontiers in Neuroinformatics 8, 68 (2014)"
 * "Breit M, Stepniewski M, Grein S, Gottmann P, Reinhardt L, Queisser G.
 *   Anatomically detailed and large-scale simulations studying synapse loss
 *   and synchrony using NeuroBox. Frontiers in Neuroanatomy 10 (2016)"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

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
        : coord_min(Aux::Numbers::inf<R>()), coord_max(-Aux::Numbers::inf<R>())
        {}

        BoundingBox(const Vec3<R>& v1, const Vec3<R>& v2)
        {
            using Aux::VecMat::minVec3;
            using Aux::VecMat::maxVec3;

            minVec3<R>(coord_min, v1, v2);
            maxVec3<R>(coord_max, v1, v2);
        }

        BoundingBox(BoundingBox const &x)
        : coord_min(x.coord_min), coord_max(x.coord_max)
        {}
    
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
