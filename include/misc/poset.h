#include <iostream>
#include <vector>


class Poset {
public:
    typedef std::pair<unsigned int, unsigned int> RT; // relation type
    typedef std::vector<RT> RV;  // relations vector

private:
    unsigned int n;  // elements will be from 0 to n-1
    RV order_relations;     // pairs of form a <= b

public:
    Poset(unsigned int _n, RV& _order_relations) : 
        n(_n), order_relations(verify(_order_relations, n))
    {
    }


    // verify if the relations are valid
    static RV verify(const RV& relations, unsigned int n)
    { 
        for (int i = 0; i < relations.size(); i++) {
            if (relations[i].first < 0 || relations[i].first >= n)
                throw "invalid elements in order relations";

            if (relations[i].second < 0 || relations[i].second >= n)
                throw "invalid elements in order relations";
        }

        // TODO: Check if corresponding DAG is actually acyclic
        
        return relations;
    }


    void print()
    {
        std::cout << "Number of elements: " << n << std::endl;
        for(int i=0; i<n; ++i) {
            std::cout << order_relations[i].first << " <= " << order_relations[i].second << std::endl;
        }
    }


    unsigned int num_elem()
    {
        return n;
    }


    unsigned int num_relations()
    {
        return order_relations.size();
    }


    RT get_relation(unsigned int idx)
    {
        return order_relations[idx];
    }

    template <typename PT_arr, typename NT>
    bool is_in(const PT_arr& pt_coeffs, NT tol=NT(0))
    {
        for (int i = 0; i < order_relations.size(); i++) {
            unsigned int a = order_relations[i].first;
            unsigned int b = order_relations[i].second;
            
            if (! (pt_coeffs(a) - pt_coeffs(b) <= tol) ) {
                return false;
            }
        }

        return true;
    }
};
