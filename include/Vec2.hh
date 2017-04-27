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

#ifndef VEC2_H
#define VEC2_H


#include "StaticVector.hh"

class Vec2
: public StaticVector<2, double>
{
    public:
        typedef StaticVector<2, double> base_type;

        // constructors
        Vec2();
        explicit Vec2(double v);
        Vec2(double x, double y);
        Vec2(const Vec2& v);

        // destructor
        ~Vec2();

        // assignment operator
        Vec2& operator=(const Vec2& v);

        // resizing (only for compatibility)
        void resize(uint32_t size);

        // arithmetic
        Vec2 operator+(const Vec2& v) const;
        Vec2& operator+=(const Vec2& v);
        Vec2 operator-(const Vec2& v) const;
        Vec2& operator-=(const Vec2& v);

        Vec2 operator*(double x) const;
        Vec2& operator*=(double x);
        Vec2 operator/(double x) const;
        Vec2& operator/=(double x);

        // scalar product, cross product
        double operator*(const Vec2& v) const;
        double cross(const Vec2& v) const;

        // relational operators
        bool operator==(const Vec2& v) const;
        bool operator!=(const Vec2& v);

        // tuple-like comparison operators
        bool operator<(const Vec2& v) const;
        bool operator>(const Vec2& v) const;
        bool operator>=(const Vec2& v) const;
        bool operator<=(const Vec2& v) const;

        // norm
        double len2(void) const;
        double len2squared(void) const;

        Vec2 &normalize();

        // output
        void print() const;
        void print_debugl(uint32_t level) const;

    private:
        using StaticVector<2, double>::v;
};


#endif
