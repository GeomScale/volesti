// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis
// Copyright (c) 2020-2021 Vaibhav Thakkar

// Contributed and/or modified by Vaibhav Thakkar, as part of Google Summer of Code 2021 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef POSET_H
#define POSET_H

#include <iostream>
#include <vector>
#include <queue>

/// A class to represent a partial order set aka Poset
class Poset {
public:
    typedef std::pair<unsigned int, unsigned int> RT; // relation type
    typedef std::vector<RT> RV;  // relations vector

private:
    unsigned int n;  // elements will be from 0 to n-1
    RV order_relations;     // pairs of form a <= b

    static void sorted_list(const unsigned int &n, const RV &relations, std::vector<unsigned int> &res)
    {
        std::vector<std::vector<unsigned int> > adjList(n);
        std::vector<unsigned int> indegree(n, 0);

        for(auto x: relations) {
            adjList[x.first].push_back(x.second);
            indegree[x.second] += 1;
        }

        std::queue<unsigned int> q;
        for(unsigned int i=0; i<n; ++i) {
            if(indegree[i] == 0)
                q.push(i);
        }

        while(!q.empty()) {
            unsigned int curr = q.front();
            res.push_back(curr);
            q.pop();

            for(unsigned int i=0; i<adjList[curr].size(); ++i) {
                unsigned int adj_idx = adjList[curr][i];
                indegree[adj_idx] -= 1;
                if(indegree[adj_idx] == 0)
                    q.push(adj_idx);
            }
        }
    }

public:
    Poset() {}

    Poset(unsigned int _n, RV& _order_relations) :
        n(_n), order_relations(verify(_order_relations, n))
    {}


    // verify if the relations are valid
    static RV verify(const RV& relations, unsigned int n)
    {
        for (int i = 0; i < relations.size(); i++) {
            if (relations[i].first < 0 || relations[i].first >= n)
                throw "invalid elements in order relations";

            if (relations[i].second < 0 || relations[i].second >= n)
                throw "invalid elements in order relations";
        }

        std::vector<unsigned int> order;
        sorted_list(n, relations, order);
        
        if(order.size() < n) { // TODO: accept cycles in the poset
            throw "corresponding DAG is not acyclic";
        }

        return relations;
    }


    void print()
    {
        std::cout << "Number of elements: " << num_elem() << std::endl;

        unsigned int count = num_relations();
        for(int i=0; i<count; ++i) {
            std::cout <<"A" << order_relations[i].first << " <= A" << order_relations[i].second << std::endl;
        }
    }


    unsigned int num_elem() const
    {
        return n;
    }


    unsigned int num_relations() const
    {
        return order_relations.size();
    }


    RT get_relation(unsigned int idx) const
    {
        return order_relations[idx];
    }

    template <typename NT>
    bool is_in(const std::vector<NT>& pt_coeffs, NT tol=NT(0)) const
    {
        for (int i = 0; i < order_relations.size(); i++) {
            unsigned int a = order_relations[i].first;
            unsigned int b = order_relations[i].second;

            if (! (pt_coeffs[a] - pt_coeffs[b] <= tol) ) {
                return false;
            }
        }

        return true;
    }


    std::vector<unsigned int> topologically_sorted_list() const
    {
        std::vector<unsigned int> res;
        sorted_list(n, order_relations, res);
        return res;
    }
};

#endif