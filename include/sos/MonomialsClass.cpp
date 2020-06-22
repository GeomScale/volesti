#include "MonomialsClass.h"

std::vector<std::vector<int> > MonomialsClass::constructAllMonomials(int const max_degree, int const num_variables) {

    std::vector<std::vector<int>> to_return;
    if (1 == num_variables) {
        for (int i = 0; i <= max_degree; i++) {
            to_return.push_back({i});
        }
        return to_return;
    }

    for (int i = 0; i <= max_degree; i++) {
        std::vector<std::vector<int>> all_options_with_this_var_at_degree_i = constructAllMonomials(max_degree - i, num_variables - 1);
        for (int j = 0; j < all_options_with_this_var_at_degree_i.size(); j++) {
            all_options_with_this_var_at_degree_i.at(j).insert(all_options_with_this_var_at_degree_i.at(j).begin(), i);
        }
        to_return.insert(to_return.end(), all_options_with_this_var_at_degree_i.begin(),
                         all_options_with_this_var_at_degree_i.end());

    }
    return to_return;
}
