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

#ifndef UTUPLE_H
#define UTUPLE_H

template <typename T>
class uTuple {
    public:
        uTuple(std::vector<T> items)
        {
            this->n     = items.size();
            this->items = items; 
            std::sort(this->items.begin(), this->items.end());
        }

        size_t
        size() const
        {
            return this->items.size();
        }

        void
        assign(std::vector<T> items) {
            this->items = items;
            std::sort(this->items.begin(), this->items.end());
        }

        bool
        operator<(const uTuple<T> &x) const
        {
            this->checkSize(x);
            return (this->items < x.items);
        }

        bool
        operator==(const uTuple<T> &x) const
        {
            this->checkSize(x);
            return (this->items == x.items);
        }

        bool
        operator<=(const uTuple<T> &x) const
        {
            this->checkSize(x);
            return (this->items < x.items || this->items == x.items);
        }

        bool
        operator!=(const uTuple<T> &x) const
        {
            this->checkSize(x);
            return (!this->operator==(x));
        }

        bool
        operator >(const uTuple<T> &x) const
        {
            this->checkSize(x);
            return (!this->operator<=(x));
        }

        bool
        operator>=(const uTuple<T> &x) const
        {
            this->checkSize(x);
            return (!this->operator<(x));
        }

        T &
        operator[](uint32_t i) const
        {
            this->checkRange(i);
            return this->items[i];
        }

    protected:
        std::vector<T>  items;
        uint32_t        n;

        inline void
        checkSize(const uTuple<T> &x) const
        {
            if (this->size() != x.size()) {
                throw("uTuple: size mismatch.");
            }
        }

        inline void
        checkRange(uint32_t i) const
        {
            if (i >= this->n) {
                throw("uTuple: out of range..");
            }
        }
};

template <typename T>
class uPair : public uTuple<T> {
    public:
        uPair(T a, T b) : uTuple<T>({a, b})
        {
        }

        void
        assign(T a, T b) {
            uTuple<T>::assign({a, b});
        }

        T
        low() const
        {
            return (this->items[0]);
        }

        T
        high() const
        {
            return (this->items[1]);
        }

        T
        fst() const
        {
            return (this->low());
        }

        T
        snd()
        {
            return (this->high());
        }
};

template <typename T>
class uTriple : public uTuple<T> {
    public:
        uTriple(T a, T b, T c) : uTuple<T>({a, b, c})
        {
        }

        void
        assign(T a, T b, T c) {
            uTuple<T>::assign({a, b, c});
        }

        T
        low() const
        {
            return (this->items[0]);
        }

        T
        mid() const
        {
            return (this->items[1]);
        }

        T
        high() const
        {
            return (this->items[2]);
        }

        T
        fst() const
        {
            return (this->low());
        }

        T
        snd()
        {
            return (this->mid());
        }

        T
        thd()
        {
            return (this->high());
        }
};

#endif
