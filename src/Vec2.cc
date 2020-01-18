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

#include "../include/Vec2.hh"
#include "debug.hh"


Vec2::Vec2()
{}

Vec2::Vec2(double _v)
{
    v[0] = _v;
    v[1] = _v;
}

Vec2::Vec2(double x, double y)
{
    v[0] = x;
    v[1] = y;
}


Vec2::Vec2(const Vec2& _v)
: StaticVector<2, double>()
{
    v[0] = _v[0];
    v[1] = _v[1];
}


Vec2::~Vec2()
{}


Vec2& Vec2::operator=(const Vec2& _v)
{
    v[0] = _v[0];
    v[1] = _v[1];

    return *this;
}


void Vec2::resize(uint32_t size)
{
    if (size != 2)
        throw("Vec2::resize(size_t): given size != 2.");
}


Vec2 Vec2::operator+(const Vec2& _v) const
{
    Vec2 r;
    r[0] = v[0] + _v[0];
    r[1] = v[1] + _v[1];
    return r;
}

Vec2& Vec2::operator+=(const Vec2& _v)
{
    v[0] += _v[0];
    v[1] += _v[1];
    return *this;
}

Vec2 Vec2::operator-(const Vec2& _v) const
{
    Vec2 r;
    r[0] = v[0] - _v[0];
    r[1] = v[1] - _v[1];
    return r;
}

Vec2& Vec2::operator-=(const Vec2& _v)
{
    v[0] -= _v[0];
    v[1] -= _v[1];
    return *this;
}


Vec2 Vec2::operator*(double x) const
{
    Vec2 r;
    r[0] = v[0] * x;
    r[1] = v[1] * x;
    return r;
}

Vec2& Vec2::operator*=(double x)
{
    v[0] *= x;
    v[1] *= x;
    return *this;
}

Vec2 Vec2::operator/(double x) const
{
    Vec2 r;
    r[0] = v[0] / x;
    r[1] = v[1] / x;
    return r;
}

Vec2& Vec2::operator/=(double x)
{
    v[0] /= x;
    v[1] /= x;
    return *this;
}


double Vec2::operator*(const Vec2& _v) const
{
    return v[0]*_v[0] + v[1]*_v[1];
}

double Vec2::cross(const Vec2& _v) const
{
    return v[0]*_v[1] - v[1]*_v[0];
}


bool Vec2::operator==(const Vec2& _v) const
{
    return v[0] == _v[0] && v[1] == _v[1];
}

bool Vec2::operator!=(const Vec2& _v)
{
    return v[0] != _v[0] || v[1] != _v[1];
}


bool Vec2::operator<(const Vec2& _v) const
{
    return v[0] < _v[0] && v[1] < _v[1];
}

bool Vec2::operator>(const Vec2& _v) const
{
    return v[0] > _v[0] && v[1] > _v[1];
}

bool Vec2::operator>=(const Vec2& _v) const
{
    return v[0] >= _v[0] && v[1] >= _v[1];
}

bool Vec2::operator<=(const Vec2& _v) const
{
    return v[0] <= _v[0] && v[1] <= _v[1];
}


double Vec2::len2(void) const
{
    return sqrt(v[0]*v[0] + v[1]*v[1]);
}

double Vec2::len2squared(void) const
{
    return v[0]*v[0] + v[1]*v[1];
}


Vec2& Vec2::normalize()
{
    double n = sqrt(v[0]*v[0] + v[1]*v[1]);
    v[0] /= n;
    v[1] /= n;
    return *this;
}


void Vec2::print() const
{
    printf("(%f, %f)\n", v[0], v[1]);
}

void Vec2::print_debugl(uint32_t level) const
{
    debugTabInc();
    debugl(level, "(%f, %f)\n", v[0], v[1]);
    debugTabDec();
}


