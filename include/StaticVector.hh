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

#ifndef STATIC_VECTOR_HH
#define STATIC_VECTOR_HH

#include <cstdint>	// uint32_t

template<uint32_t N, typename T = double>
class StaticVector
{
    public:
        typedef StaticVector<N, T> this_type;

        StaticVector();
        explicit StaticVector(const T& x);
        StaticVector(const this_type& v);

        ~StaticVector();

        this_type& operator=(const this_type& v);

        uint32_t size() const;

        void assign(const T& x);

        T& operator()(uint32_t i);
        T operator()(uint32_t i) const;

        T& operator[](uint32_t i);
        T operator[](uint32_t i) const;

        this_type operator+(const this_type& v) const;
        this_type& operator+=(const this_type& v);

        this_type operator-(const this_type& v) const;
        this_type& operator-=(const this_type& v);

        this_type operator*(const T& alpha) const;
        this_type& operator*=(const T& alpha);

    protected:
        T v[N];
};

#include "../tsrc/StaticVector_impl.hh"

#endif // STATIC_VECTOR_HH
