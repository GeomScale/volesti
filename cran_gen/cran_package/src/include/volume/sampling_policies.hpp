// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef SAMPLING_POLICIES_HPP
#define SAMPLING_POLICIES_HPP

struct PushBackWalkPolicy
{
    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p) const
    {
        randPoints.push_back(p);
    }
};

template <typename BallPoly>
struct CountingWalkPolicy
{
    CountingWalkPolicy(unsigned int const& nump_PBSmall, BallPoly const& PBSmall)
            :   _nump_PBSmall(nump_PBSmall)
            ,   _PBSmall(PBSmall)
    {}

    template <typename PointList, typename Point>
    void apply(PointList &randPoints,
               Point &p)
    {
        if (_PBSmall.second().is_in(p) == -1)//is in
        {
            randPoints.push_back(p);
            ++_nump_PBSmall;
        }
    }

    unsigned int get_nump_PBSmall() const
    {
        return _nump_PBSmall;
    }

private :
    unsigned int _nump_PBSmall;
    BallPoly _PBSmall;
};

#endif // SAMPLING_POLICIES_HPP
