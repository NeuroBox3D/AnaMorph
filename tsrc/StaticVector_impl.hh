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

#include <assert.h>
#include "StaticVector.hh"


template<uint32_t N, typename T>
StaticVector<N, T>::StaticVector()
{}

template<uint32_t N, typename T>
StaticVector<N, T>::StaticVector(const T& x)
{
    assign(x);
}

template<uint32_t N, typename T>
StaticVector<N, T>::StaticVector(const this_type& _v)
{
    for (uint32_t i = 0; i < N; ++i)
        v[i] = _v.v[i];
}

template<uint32_t N, typename T>
StaticVector<N, T>::~StaticVector()
{}

template<uint32_t N, typename T>
StaticVector<N, T>&
StaticVector<N, T>::operator=(const this_type& _v)
{
    for (uint32_t i = 0; i < N; ++i)
        v[i] = _v.v[i];
    return *this;
}

template<uint32_t N, typename T>
uint32_t
StaticVector<N, T>::size() const
{
    return N;
}

template<uint32_t N, typename T>
void
StaticVector<N, T>::assign(const T& x)
{
    for (uint32_t i = 0; i < N; ++i)
        v[i] = x;
}

template<uint32_t N, typename T>
T&
StaticVector<N, T>::operator()(uint32_t i)
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
T
StaticVector<N, T>::operator()(uint32_t i) const
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
T &
StaticVector<N, T>::operator[](uint32_t i)
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
T
StaticVector<N, T>::operator[](uint32_t i) const
{
    assert(i < N && "Invalid component requested.");
    return v[i];
}

template<uint32_t N, typename T>
StaticVector<N, T>
StaticVector<N, T>::operator+(const this_type& _v) const
{
    StaticVector r;
    for (uint32_t i = 0; i < N; ++i)
        r[i] = v[i] + _v[i];
    return r;
}

template<uint32_t N, typename T>
StaticVector<N, T> &
StaticVector<N, T>::operator+=(const this_type& _v)
{
    for (uint32_t i = 0; i < N; ++i)
        v[i] += _v[i];
    return *this;
}

template<uint32_t N, typename T>
StaticVector<N, T>
StaticVector<N, T>::operator-(const this_type& _v) const
{
    StaticVector r;
    for (uint32_t i = 0; i < N; ++i)
        r[i] = v[i] - _v[i];
    return r;
}

template<uint32_t N, typename T>
StaticVector<N, T> &
StaticVector<N, T>::operator-=(const this_type& _v)
{
    for (uint32_t i = 0; i < N; ++i)
        v[i] -= _v[i];
    return *this;
}

template<uint32_t N, typename T>
StaticVector<N, T>
StaticVector<N, T>::operator*(const T& alpha) const
{
    StaticVector r;
    for (uint32_t i = 0; i < N; ++i)
        r[i] = alpha * v[i];
    return r;
}

template<uint32_t N, typename T>
StaticVector<N, T> &
StaticVector<N, T>::operator*=(const T& alpha)
{
    for (uint32_t i = 0; i < N; ++i)
        v[i] *= alpha;
    return *this;
}
