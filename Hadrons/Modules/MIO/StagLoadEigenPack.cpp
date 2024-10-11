#include <Hadrons/Modules/MIO/StagLoadEigenPack.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class HADRONS_NAMESPACE::MIO::TStagLoadEigenPack<FermionEigenPack<STAGIMPL>, GIMPL, STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class HADRONS_NAMESPACE::MIO::TStagLoadEigenPack<FermionEigenPack<STAGIMPL, STAGIMPLF>, GIMPL, STAGIMPL>;
#endif

