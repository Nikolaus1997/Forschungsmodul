#include "basis_storage.h"

BasisStorage::BasisStorage(int N):
    PP_N(N),
    weights_(N+2,0.0,0.0),
    nodes_(N+2,0.0,0.0)
    {
}

const Variable &BasisStorage::weights() const
{
    return weights_;
}

const Variable &BasisStorage::nodes() const
{
    return nodes_;
}

double BasisStorage::weights(int i) const
{
    return weights_(i);
}

double &BasisStorage::weights(int i)
{
    return weights_(i);
}

double BasisStorage::nodes(int i) const
{
    return nodes_(i);
}

double &BasisStorage::nodes(int i)
{
    return nodes_(i);
}
