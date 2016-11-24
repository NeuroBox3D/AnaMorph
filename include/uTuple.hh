/* --------------------------------------------------------------------------------
 * 
 *                              THIS FILE IS PART OF                               
 * 
 * AnaMorph: A Framework for Geometric Modelling, Consistency Analysis and Surface
 * Mesh Generation of Anatomically Reconstructed Neuron Morphologies
 * 
 * Web site: http://www.anamorph.net
 * 
 * Copyright (c) 2013-2014, Konstantin Mörschel.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 * 
 *    This product includes software developed by Konstantin Mörschel for
 *    AnaMorph (http://www.anamorph.net).
 * 
 * 4. Neither the name "AnaMorph" nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without
 *    specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS OF ANAMORPH ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS OF ANAMORPH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * -------------------------------------------------------------------------------- */

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
