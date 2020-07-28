// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#ifndef NONSYMMETRICCONICOPTIMIZATION_MONOMIALSCLASS_H
#define NONSYMMETRICCONICOPTIMIZATION_MONOMIALSCLASS_H

#include<vector>
#include<iostream>

typedef std::vector<int> Monomial;
typedef std::vector<Monomial> Monomials;

class MonomialsClass {
public:
    MonomialsClass(int n, int d) : _n(n), _d(d) {
        _monomials = constructAllMonomials(_d, _n);
    };

    static Monomial prod(Monomial m1, Monomial m2) {
        assert(m1.size() == m2.size());
        Monomial m(m1.size());
        for (unsigned i = 0; i < m.size(); ++i) {
            m[i] = m1[i] + m2[i];
        }
        return m;
    }

    Monomials getMonomials() {
        return _monomials;
    }

    void print_human_readable(Monomial m) const {
        assert(m.size() == _n);
        assert(_n >= 1);
        int last_non_zero = _n - 1;
        while (last_non_zero >= 0 && m[last_non_zero] == 0) {
            last_non_zero--;
        }
        if (last_non_zero < 0) {
            std::cout << "1";
            return;
        }
        for (unsigned i = 0; i < _n; ++i) {
            if (m[i] > 0) {
                std::cout << "x_" << i + 1;
                if (m[i] > 1) {
                    std::cout << "^" << m[i];
                }
                if (i != last_non_zero) {
                    std::cout << " + ";
                }
            }
        }
    }

    void print_human_readable() {
        for (auto monomial : _monomials) {
            print_human_readable(monomial);
            std::cout << std::endl;
        }
    }

    void print() {
        for (auto monomial : _monomials) {
            for (auto i : monomial) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
    }

private:
    Monomials constructAllMonomials(unsigned const max_degree, unsigned const num_variables);

    unsigned _n;
    unsigned _d;
    Monomials _monomials;
};


#endif //NONSYMMETRICCONICOPTIMIZATION_MONOMIALSCLASS_H
