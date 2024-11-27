#include "basis_storage.h"

BasisStorage::BasisStorage(int N):
    PP_N_(N),
    weights_({PP_N_+2}),
    nodes_({PP_N_+2})
    {
}

const Array1D &BasisStorage::weights() const
{
    return weights_;
}

const Array1D &BasisStorage::nodes() const
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
