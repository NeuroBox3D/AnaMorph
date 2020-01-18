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

#include "../include/Vec3.hh"
#include "debug.hh"

template <typename R>
Vec3<R>::Vec3()
{}

template <typename R>
Vec3<R>::Vec3(R _v)
{
    v[0] = _v;
    v[1] = _v;
    v[2] = _v;
}

template <typename R>
Vec3<R>::Vec3(R x, R y, R z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}


template <typename R>
Vec3<R>::Vec3(const Vec3<R>& _v)
: StaticVector<3, R>()
{
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];
}


template <typename R>
Vec3<R>::~Vec3()
{}


template <typename R>
Vec3<R>& Vec3<R>::operator=(const Vec3<R>& _v)
{
    v[0] = _v[0];
    v[1] = _v[1];
    v[2] = _v[2];

    return *this;
}


template <typename R>
void Vec3<R>::resize(uint32_t size)
{
    if (size != 3)
        throw("Vec3<R>::resize(size_t): given size != 3.");
}


template <typename R>
Vec3<R> Vec3<R>::operator+(const Vec3<R>& _v) const
{
    Vec3<R> r;
    r[0] = v[0] + _v[0];
    r[1] = v[1] + _v[1];
    r[2] = v[2] + _v[2];
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator+=(const Vec3<R>& _v)
{
    v[0] += _v[0];
    v[1] += _v[1];
    v[2] += _v[2];
    return *this;
}

template <typename R>
Vec3<R> Vec3<R>::operator-(const Vec3<R>& _v) const
{
    Vec3<R> r;
    r[0] = v[0] - _v[0];
    r[1] = v[1] - _v[1];
    r[2] = v[2] - _v[2];
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator-=(const Vec3<R>& _v)
{
    v[0] -= _v[0];
    v[1] -= _v[1];
    v[2] -= _v[2];
    return *this;
}


template <typename R>
Vec3<R> Vec3<R>::operator*(R x) const
{
    Vec3<R> r;
    r[0] = v[0] * x;
    r[1] = v[1] * x;
    r[2] = v[2] * x;
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator*=(R x)
{
    v[0] *= x;
    v[1] *= x;
    v[2] *= x;
    return *this;
}

template <typename R>
Vec3<R> Vec3<R>::operator/(R x) const
{
    Vec3<R> r;
    r[0] = v[0] / x;
    r[1] = v[1] / x;
    r[2] = v[2] / x;
    return r;
}

template <typename R>
Vec3<R>& Vec3<R>::operator/=(R x)
{
    v[0] /= x;
    v[1] /= x;
    v[2] /= x;
    return *this;
}


template <typename R>
R Vec3<R>::operator*(const Vec3<R>& _v) const
{
    return v[0]*_v[0] + v[1]*_v[1] + v[2]*_v[2];
}

template <typename R>
Vec3<R> Vec3<R>::cross(const Vec3<R>& _v) const
{
    Vec3<R> r;
    r[0] = v[1]*_v[2] - v[2]*_v[1];
    r[1] = v[2]*_v[0] - v[0]*_v[2];
    r[2] = v[0]*_v[1] - v[1]*_v[0];
    return r;
}


template <typename R>
bool Vec3<R>::operator==(const Vec3<R>& _v) const
{
    return v[0] == _v[0] && v[1] == _v[1] && v[2] == _v[2];
}

template <typename R>
bool Vec3<R>::operator!=(const Vec3<R>& _v)
{
    return v[0] != _v[0] || v[1] != _v[1] || v[2] != _v[2];
}


template <typename R>
bool Vec3<R>::operator<(const Vec3<R>& _v) const
{
    return v[0] < _v[0] && v[1] < _v[1] && v[2] < _v[2];
}

template <typename R>
bool Vec3<R>::operator>(const Vec3<R>& _v) const
{
    return v[0] > _v[0] && v[1] > _v[1] && v[2] > _v[2];
}

template <typename R>
bool Vec3<R>::operator>=(const Vec3<R>& _v) const
{
    return v[0] >= _v[0] && v[1] >= _v[1] && v[2] >= _v[2];
}

template <typename R>
bool Vec3<R>::operator<=(const Vec3<R>& _v) const
{
    return v[0] <= _v[0] && v[1] <= _v[1] && v[2] <= _v[2];
}


template <typename R>
R Vec3<R>::len2(void) const
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

template <typename R>
R Vec3<R>::len2squared(void) const
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


template <typename R>
Vec3<R>& Vec3<R>::normalize()
{
    R n = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= n;
    v[1] /= n;
    v[2] /= n;
    return *this;
}


template <typename R>
void Vec3<R>::print() const
{
    printf("(%f, %f, %f)\n", v[0], v[1], v[2]);
}

template <typename R>
void Vec3<R>::print_debugl(uint32_t level) const
{
    debugTabInc();
    debugl(level, "(%f, %f, %f)\n", v[0], v[1], v[2]);
    debugTabDec();
}



template <typename R>
std::ostream& operator<<(std::ostream& stream, const Vec3<R>& v)
{
	if (!v.size()) return stream << "()";

	std::size_t sz = v.size() - 1;
	stream << "(";
	for (std::size_t i = 0; i < sz; ++i)
		stream << v[i] << " ";
	stream << v[sz] << ")";

	return stream;
}


// make Vec3<double> usable
template class Vec3<double>;
template std::ostream& operator<<(std::ostream& stream, const Vec3<double>& v);
