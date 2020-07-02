// VolEsti (volume computation and sampling library)
//
// Copyright (c) 2020 Bento Natura
//
// Licensed under GNU LGPL.3, see LICENCE file

#include "MonomialsClass.h"

std::vector<std::vector<int> > MonomialsClass::constructAllMonomials(unsigned const max_degree,
                                                                     unsigned const num_variables) {

    std::vector<std::vector<int>> monomial_vec;
    if (1 == num_variables) {
        for (unsigned i = 0; i <= max_degree; i++) {
            monomial_vec.push_back({(int) i});
        }
        return monomial_vec;
    }

    for (unsigned i = 0; i <= max_degree; i++) {
        std::vector<std::vector<int>> all_options_with_this_var_at_degree_i = constructAllMonomials(max_degree - i,
                                                                                                    num_variables - 1);
        for (unsigned j = 0; j < all_options_with_this_var_at_degree_i.size(); j++) {
            all_options_with_this_var_at_degree_i.at(j).insert(all_options_with_this_var_at_degree_i.at(j).begin(), i);
        }
        monomial_vec.insert(monomial_vec.end(), all_options_with_this_var_at_degree_i.begin(),
                            all_options_with_this_var_at_degree_i.end());
    }
    return monomial_vec;
}
