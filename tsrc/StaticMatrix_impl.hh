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
#include "StaticMatrix.hh"


template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>::StaticMatrix()
{}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>::~StaticMatrix()
{}

template<uint32_t M, uint32_t N, typename T>
uint32_t
StaticMatrix<M, N, T>::numRows() const
{
    return M;
}

template<uint32_t M, uint32_t N, typename T>
uint32_t
StaticMatrix<M, N, T>::numCols() const
{
    return N;
}

template<uint32_t M, uint32_t N, typename T>
void StaticMatrix<M, N, T>::fill(const T& x)
{
    for (size_t i = 0; i < M; ++i)
        m[i].assign(x);
}

template<uint32_t M, uint32_t N, typename T>
T&
StaticMatrix<M, N, T>::operator()(uint32_t i, uint32_t j)
{
    assert(i < M && j < N && "Invalid component requested.");
    return m[i][j];
}

template<uint32_t M, uint32_t N, typename T>
T
StaticMatrix<M, N, T>::operator()(uint32_t i, uint32_t j) const
{
    assert(i < M && j < N && "Invalid component requested.");
    return m[i][j];
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>
StaticMatrix<M, N, T>::operator+(const this_type& _m) const
{
    this_type r;
    for (size_t i = 0; i < M; ++i)
        r.m[i] = m[i] + _m.m[i];
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T> &
StaticMatrix<M, N, T>::operator+=(const this_type& _m)
{
    for (size_t i = 0; i < M; ++i)
        m[i] += _m.m[i];
    return *this;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>
StaticMatrix<M, N, T>::operator-(const this_type& _m) const
{
    this_type r;
    for (size_t i = 0; i < M; ++i)
        r.m[i] = m[i] - _m.m[i];
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T> &
StaticMatrix<M, N, T>::operator-=(const this_type& _m)
{
    for (size_t i = 0; i < M; ++i)
        m[i] -= _m.m[i];
    return *this;
}

template<uint32_t M, uint32_t N, typename T>
template <uint32_t L>
StaticMatrix<M, L, T>
StaticMatrix<M, N, T>::operator*(const StaticMatrix<N, L, T>& _m) const
{
    StaticMatrix<M, L, T> r;
    for (size_t i = 0; i < M; ++i)
    {
        for (size_t j = 0; j < L; ++j)
        {
            T& t = r(i,j);
            for (size_t k = 0; k < N; ++j)
                t += m[i][k] * _m(k,j);
        }
    }
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T>
StaticMatrix<M, N, T>::operator*(const T& alpha) const
{
    this_type r;
    for (size_t i = 0; i < M; ++i)
        r[i] = alpha * m[i];
    return r;
}

template<uint32_t M, uint32_t N, typename T>
StaticMatrix<M, N, T> &
StaticMatrix<M, N, T>::operator*=(const T& alpha)
{
    for (size_t i = 0; i < M; ++i)
        m[i] *= alpha;
    return *this;
}

template<uint32_t M, uint32_t N, typename T>
const StaticVector<N, T>&
StaticMatrix<M, N, T>::getRow(uint32_t i) const
{
    return m[i];
}

template<uint32_t M, uint32_t N, typename T>
StaticVector<N, T>&
StaticMatrix<M, N, T>::getRow(uint32_t i)
{
    return m[i];
}

template<uint32_t M, uint32_t N, typename T>
StaticVector<M, T>
StaticMatrix<M, N, T>::getCol(uint32_t j) const
{
    StaticVector<M, T> r;
    for (size_t i = 0; i < M; ++i)
        r[i] = m[i][j];
    return r;
}

