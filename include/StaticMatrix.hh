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

#ifndef STATIC_MATRIX_HH
#define STATIC_MATRIX_HH

#include "StaticVector.hh"

template<uint32_t M, uint32_t N, typename T = double>
class StaticMatrix
{
    public:
        typedef StaticMatrix<M, N, T> this_type;
        typedef StaticVector<N, T> row_type;
        typedef StaticVector<M, T> col_type;

        StaticMatrix();
        ~StaticMatrix();

        uint32_t numRows() const;
        uint32_t numCols() const;

        void fill(const T& x);

        T& operator()(uint32_t i, uint32_t j);
        T operator()(uint32_t i, uint32_t j) const;

        this_type operator+(const this_type& m) const;
        this_type& operator+=(const this_type& m);

        this_type operator-(const this_type& m) const;
        this_type& operator-=(const this_type& m);

        template <uint32_t L>
        StaticMatrix<M, L, T> operator*(const StaticMatrix<N, L, T>& m) const;

        this_type operator*(const T& alpha) const;
        this_type& operator*=(const T& alpha);

        const StaticVector<N, T>& getRow(uint32_t i) const;
        StaticVector<N, T>& getRow(uint32_t i);
        StaticVector<M, T> getCol(uint32_t j) const;

    protected:
        StaticVector<N,T> m[M];
};

#include "../tsrc/StaticMatrix_impl.hh"

#endif // STATIC_MATRIX_HH
