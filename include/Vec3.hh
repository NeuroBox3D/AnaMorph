/*
 * This file is part of
 *
 * AnaMorph: a framework for geometric modelling, consistency analysis and surface
 * mesh generation of anatomically reconstructed neuron morphologies.
 * 
 * Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group
 * Author: Markus Breit
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

#ifndef VEC3_H
#define VEC3_H

#include "StaticVector.hh"
#include <ostream>


template <typename R>
class Vec3
: public StaticVector<3, R>
{
    public:
        typedef StaticVector<3, R> base_type;

        // constructors
        Vec3();
        Vec3(R v);
        Vec3(R x, R y, R z);
        Vec3(const Vec3& v);

        // destructor
        ~Vec3();

        // assignment operator
        Vec3& operator=(const Vec3& v);

        // resizing (only for compatibility)
        void resize(uint32_t size);

        // arithmetic
        Vec3 operator+(const Vec3& v) const;
        Vec3& operator+=(const Vec3& v);
        Vec3 operator-(const Vec3& v) const;
        Vec3& operator-=(const Vec3& v);

        Vec3 operator*(R x) const;
        Vec3& operator*=(R x);
        Vec3 operator/(R x) const;
        Vec3& operator/=(R x);

        // scalar product, cross product
        R operator*(const Vec3& v) const;
        Vec3 cross(const Vec3& v) const;

        // relational operators
        bool operator==(const Vec3& v) const;
        bool operator!=(const Vec3& v);

        // tuple-like comparison operators
        bool operator<(const Vec3& v) const;
        bool operator>(const Vec3& v) const;
        bool operator>=(const Vec3& v) const;
        bool operator<=(const Vec3& v) const;

        // norm
        R len2(void) const;
        R len2squared(void) const;

        Vec3 &normalize();

        // output
        void print() const;
        void print_debugl(uint32_t level) const;

    private:
        using StaticVector<3, R>::v;
};

template <typename R>
std::ostream& operator<<(std::ostream& stream, const Vec3<R>& v);

#endif
