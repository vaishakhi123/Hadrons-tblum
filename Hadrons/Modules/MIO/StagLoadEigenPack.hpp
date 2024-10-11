#ifndef Hadrons_MIO_StagLoadEigenPack_hpp_
#define Hadrons_MIO_StagLoadEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Solver.hpp> //needed for FMat


BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         StagLoadEigenPack                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)


class StagLoadEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagLoadEigenPackPar,
                                    std::string, filestem,
                                    bool, multiFile,
                                    bool, redBlack,
                                    unsigned int, size,
                                    unsigned int, Ls,
                                    std::string, gaugeXform,
                                    std::string, action,
                                    double, mass);
};

template <typename Pack, typename GImpl, typename FImpl>
class TStagLoadEigenPack: public Module<StagLoadEigenPackPar>
{
public:
    typedef FermionOperator<FImpl> FMat;
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef BaseEigenPack<Field>   BasePack;

public:
    GAUGE_TYPE_ALIASES(GImpl, );
    typedef typename GImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TStagLoadEigenPack(const std::string name);
    // destructor
    virtual ~TStagLoadEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagLoadFermionEigenPack, ARG(TStagLoadEigenPack<FermionEigenPack<STAGIMPL>, GIMPL, STAGIMPL>), MIO);

/******************************************************************************
 *                    TStagLoadEigenPack implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
TStagLoadEigenPack<Pack, GImpl, FImpl>::TStagLoadEigenPack(const std::string name)
: Module<StagLoadEigenPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
std::vector<std::string> TStagLoadEigenPack<Pack, GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().gaugeXform.empty())
    {
        in = {par().gaugeXform};
    }
    in.push_back(par().action);
    
    return in;
}

template <typename Pack, typename GImpl, typename FImpl>
std::vector<std::string> TStagLoadEigenPack<Pack, GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
void TStagLoadEigenPack<Pack, GImpl, FImpl>::setup(void)
{
    GridBase  *grid, *gridIo = nullptr, *gridRb = nullptr;

    grid   = getGrid<Field>(par().Ls);
    gridRb = getGrid<Field>(par().redBlack, par().Ls);
    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = getGrid<FieldIo>(par().redBlack, par().Ls);
    }
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().size, gridRb, gridIo);
    if (!par().gaugeXform.empty())
    {
        envTmp(GaugeMat,    "tmpXform", par().Ls, grid);
        if (par().redBlack)
        {
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, gridRb);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl, typename FImpl>
void TStagLoadEigenPack<Pack, GImpl, FImpl>::execute(void)
{
    
    auto &action= envGet(FMat, par().action);
    auto &epack = envGetDerived(BasePack, Pack, getName());

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    epack.eval.resize(par().size);
    
    ComplexD minusI(0, -1.0);
    ComplexD cc;
    RealD eval;
    double mass = par().mass;
    Field temp(epack.evec[0].Grid());
    
    // make even evecs from Odd
//    if(par().redBlack==Even){
//        for (unsigned int i = 0; i < par().size; i++)
//        {
//            eval=sqrt(epack.eval[i]-mass*mass);
//            epack.evec[i].Checkerboard() = Odd;
//            action.Meooe(epack.evec[i], temp);
//            cc = minusI/eval;
//            epack.evec[i] = cc * temp; // now it's even!
//            epack.evec[i].Checkerboard() = Even;
//        }
//    } else {
        for (unsigned int i = 0; i < par().size; i++)
            epack.evec[i].Checkerboard() = Odd;
    //}
        
    if (!par().gaugeXform.empty())
    {

        LOG(Message) << "Applying gauge transformation to eigenvectors " << getName()
                     << " using " << par().gaugeXform << std::endl;
        auto &xform = envGet(GaugeMat, par().gaugeXform);
        envGetTmp(GaugeMat,    tmpXform);
        envGetTmp(GaugeMat, tmpXformOdd);

        if (par().Ls > 1)
        {
            LOG(Message) << "Creating 5d GaugeMat from " << par().gaugeXform << std::endl;
            startTimer("5-d gauge transform creation");
            for (unsigned int j = 0; j < par().Ls; j++)
            {
                InsertSlice(xform, tmpXform, j, 0);
            }
            stopTimer("5-d gauge transform creation");
        }

        pickCheckerboard(Odd, tmpXformOdd, tmpXform);
        startTimer("Transform application");
        for (unsigned int i = 0; i < par().size; i++)
        {
            LOG(Message) << "Applying gauge transformation to eigenvector i = " << i << "/" << par().size << std::endl;
            epack.evec[i].Checkerboard() = Odd;
            epack.evec[i] = tmpXformOdd * epack.evec[i];
        }
        stopTimer("Transform application");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_StagLoadEigenPack_hpp_
