/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                  std::string, noise,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver,
                                  std::string, output,
                                  double, mass,
                                  bool,        multiFile);
};

template <typename FImpl, typename Pack>
class TA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(A2AVectors, 
    ARG(TA2AVectors<FIMPL, BaseFermionEigenPack<FIMPL>>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, 
    ARG(TA2AVectors<ZFIMPL, BaseFermionEigenPack<ZFIMPL>>), MSolver);

/******************************************************************************
 *                       TA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectors<FImpl, Pack>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;

    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver + sub_string);
    in.push_back(par().noise);

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver + sub_string);
    int         Ls          = env().getObjectLs(par().action);

    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::execute(void)
{
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver + sub_string);
    auto        &noise     = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);

    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using eigenpack '" << par().eigenPack << "' ("
                     << Nl_ << " low modes) and noise '"
                     << par().noise << "' (" << noise.fermSize() 
                     << " noise vectors)" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using noise '" << par().noise << "' (" << noise.fermSize() 
                     << " noise vectors)" << std::endl;
    }
    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
        auto &epack  = envGet(Pack, par().eigenPack);

        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[il], f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[il], f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("W low mode");
    }

    // High modes
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
}

// Staggered A2A vectors
//
class StagSparseA2AVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AVectorsPar,
                                    std::string, action,
                                    std::string, gauge,
                                    std::string, eigenPack,
                                    std::string, output,
                                    int, inc,
                                    int, tinc,
                                    double, mass,
                                    bool, multiFile);
};

//
// Sparsened A2A vectors for staggered conserved current
//
template <typename FImpl, typename Pack>
class TStagSparseA2AVectors : public Module<StagSparseA2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsLowStaggered<FImpl> A2A;
    typedef typename Grid::NaiveStaggeredFermionD::FermionField SparseFermionField;
    
public:
    // constructor
    TStagSparseA2AVectors(const std::string name);
    // destructor
    virtual ~TStagSparseA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagSparseA2AVectors,
                    ARG(TStagSparseA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),MSolver);

/******************************************************************************
 *                       TStagSparseA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagSparseA2AVectors<FImpl, Pack>::TStagSparseA2AVectors(const std::string name)
: Module<StagSparseA2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagSparseA2AVectors<FImpl, Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    in.push_back(par().eigenPack);
    in.push_back(par().gauge);
    in.push_back(par().action);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagSparseA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() +"_v", getName() +"_w0", getName() +"_w1", getName() +"_w2"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagSparseA2AVectors<FImpl, Pack>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    
    auto &epack = envGet(Pack, par().eigenPack);
    Nl_ = epack.evec.size();
    envTmp(A2A, "a2a", 1, action);
    
    // Sparse Grid
    std::vector<int> blocksize(4);
    blocksize[0] = par().inc;
    blocksize[1] = par().inc;
    blocksize[2] = par().inc;
    blocksize[3] = par().tinc;
    envCreate(std::vector<SparseFermionField>, getName() + "_v", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w0", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w1", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w2", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagSparseA2AVectors<FImpl, Pack>::execute(void)
{
    auto    &action = envGet(FMat, par().action);
    auto    &U      = envGet(LatticeGaugeField, par().gauge);
    auto    &epack  = envGet(Pack, par().eigenPack);
    double  mass    = par().mass;
    uint64_t     nt      = env().getDim(Tp);
    uint64_t     ns      = env().getDim(Xp);
    uint64_t     glbsize = ns*ns*ns * nt;
    envGetTmp(A2A, a2a);
    
    int orthogdim=env().getNd()-1; // time dir
        
    auto &v = envGet(std::vector<SparseFermionField>, getName() + "_v");
    auto &w0 = envGet(std::vector<SparseFermionField>, getName() + "_w0");
    auto &w1 = envGet(std::vector<SparseFermionField>, getName() + "_w1");
    auto &w2 = envGet(std::vector<SparseFermionField>, getName() + "_w2");
    ScidacWriter binWriter(U.Grid()->IsBoss());
    std::string fullFilename;
    A2AVectorsIo::Record record;
    const int traj=vm().getTrajectory();
    
    LOG(Message) << "Computing all-to-all vectors using eigenpack " << par().eigenPack << " with " << 2*Nl_ << " low modes " << std::endl;

    // Staggered Phases. Do spatial gamma only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
        
    ComplexField phases(U.Grid());
    int shift;
    FermionField temp(U.Grid());
    FermionField temp2(U.Grid());
    //SparseFermionField temp3(&sparseGrid);
    
    // step size for hypercube loop
    int step=2*par().inc;
    //std::random_device rd;  // a seed source for the random number engine
    //std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    // sparsen on time slice starting anywhere between 0 and ns-1, inclusive
    // seed with rngSerial (below)
    std::uniform_int_distribution<uint32_t> uid(0, step-1);
    std::vector<uint32_t> xshift(nt);
    std::vector<uint32_t> yshift(nt);
    std::vector<uint32_t> zshift(nt);
    if(par().inc != 1){
        for(int t=0;t<nt;t++){
            xshift[t]=uid(rngSerial()._generators[0]);
            yshift[t]=uid(rngSerial()._generators[0]);
            zshift[t]=uid(rngSerial()._generators[0]);
        }
    }else{//will miss site 0000 otherwise
        for(int t=0;t<nt;t++){
            xshift[t]=0;
            yshift[t]=0;
            zshift[t]=0;
        }
    }
    CartesianCommunicator::BroadcastWorld(0,(void *)&xshift[0],sizeof(uint32_t)*xshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&yshift[0],sizeof(uint32_t)*yshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&zshift[0],sizeof(uint32_t)*zshift.size());
    
    //save for later
    std::vector<complex<double>> evalM(2*Nl_);
    
    // global to local coord shifts
    int locx=U.Grid()->_ldimensions[0];
    int locy=U.Grid()->_ldimensions[1];
    int locz=U.Grid()->_ldimensions[2];
    int loct=U.Grid()->_ldimensions[3];
    int lstartx=U.Grid()->_lstart[0];
    int lstarty=U.Grid()->_lstart[1];
    int lstartz=U.Grid()->_lstart[2];
    int lstartt=U.Grid()->_lstart[3];
    printf("my node %d local tsites= %d\n", U.Grid()->_processor_coor[3],loct);
    printf("my node %d global t start= %d\n", U.Grid()->_processor_coor[3], lstartt);
    LOG(Message) << "xshift" << xshift << std::endl;
    LOG(Message) << "yshift" << yshift << std::endl;
    LOG(Message) << "zshift" << zshift << std::endl;

    for (unsigned int il = 0; il < 2*Nl_; il++)
    {
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il/2]-mass*mass));
        
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low modes)" << std::endl;
        // don't divide by lambda. Do it in contraction since it is complex
        a2a.makeLowModeW(temp, epack.evec[il/2], eval, il%2);
        stopTimer("W low mode");
        
        il%2 ? eval=conjugate(eval) : eval ;
        evalM[il]=eval;
        
        for (int mu=0;mu<3;mu++){
            
            phases=1.0;
            if(mu==1){
                phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            } else if(mu==2){
                phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            }
            LatticeColourMatrix Umu(U.Grid());
            Umu = PeekIndex<LorentzIndex>(U,mu);
            Umu *= phases;
        
            // v vec is shifted and * link for conserved current
            temp2 = Umu*Cshift(temp, mu, 1);
                        
            thread_for(t,loct,{
            //for(int t=0; t<loct;t++){//debug
                int tglb=t+lstartt;
                // same random shift for t, t+1 in same hypercube
                if(t%2 == 1) continue;
                
                Coordinate site(Nd);
                Coordinate sparseSite(Nd);
                ColourVector vec;
                
                
                // loop over hypercubes on time slice
                // and sparsen, keep all sites in hypercube
                for(int z=0;z<ns;z+=step){
                    int zg=(zshift[tglb]+z)%ns;
                    for(int zl=0;zl<locz;zl++){
                        int zgp=zl+lstartz;
                        if(zgp==zg || zgp==(zg+1)%ns){
                            site[2]=zl;
                            if(par().inc==1){
                                sparseSite[2]=site[2];
                            }else if(zshift[tglb]!=step-1){
                                sparseSite[2]=2*int(site[2]/step) + (site[2]+zshift[tglb])%2;
                            }else{
                                sparseSite[2]=2*int(site[2]/step) + (site[2])%2;
                            }
                            for(int y=0;y<ns;y+=step){
                                int yg=(yshift[tglb]+y)%ns;
                                for(int yl=0;yl<locy;yl++){
                                    int ygp=yl+lstarty;
                                    if(ygp==yg || ygp==(yg+1)%ns){
                                        site[1]=yl;
                                        if(par().inc==1){
                                            sparseSite[1]=site[1];
                                        }else if(yshift[tglb]!=step-1){
                                            sparseSite[1]=2*int(site[1]/step) + (site[1]+yshift[tglb])%2;
                                        }else{
                                            sparseSite[1]=2*int(site[1]/step) + (site[1])%2;
                                        }
                                        for(int x=0;x<ns;x+=step){
                                            int xg=(xshift[tglb]+x)%ns;
                                            for(int xl=0;xl<locx;xl++){
                                                int xgp=xl+lstartx;
                                                if(xgp==xg || xgp==(xg+1)%ns){
                                                    site[0]=xl;
                                                    if(par().inc==1){
                                                        sparseSite[0]=site[0];
                                                    }else if(xshift[tglb]!=step-1){
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0]+xshift[tglb])%2;
                                                    }else{
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0])%2;
                                                    }
                                                    for(int that=0;that<2;that++){
                                                        
                                                        site[3]=t+that;
                                                        sparseSite[3]=site[3];
                                                        
                                                        if(mu==0){// do v once
                                                            peekLocalSite(vec,temp,site);
                                                            pokeLocalSite(vec,v[il],sparseSite);
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w0[il],sparseSite);
                                                        }else if(mu==1){
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w1[il],sparseSite);
                                                        }else if(mu==2){
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w2[il],sparseSite);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }// end mu
    }// end evecs
    
    std::string dir = dirname(par().output);
    if (!par().output.empty())
    {
        int status = mkdir(dir);
        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' ( " + std::strerror(errno) + ")");
        }
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w0", w0, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w1", w1, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w2", w2, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
    
    if ( env().getGrid()->IsBoss() ) {
        
        std::string eval_filename;
        if (!par().output.empty()){
            eval_filename=A2AVectorsIo::evalFilename(par().output,vm().getTrajectory());
        }else{
            eval_filename=A2AVectorsIo::evalFilename("evals",vm().getTrajectory());
        }
        A2AVectorsIo::initEvalFile(eval_filename,
                                   evalM.size());// total size
        A2AVectorsIo::saveEvalBlock(eval_filename,
                                    evalM.data(),
                                    0,// start of chunk
                                    2*Nl_);// size of chunk saved
    }
}

template <typename FImpl, typename Pack>
class TStagA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchurStaggered<FImpl> A2A;
public:
    // constructor
    TStagA2AVectors(const std::string name);
    // destructor
    virtual ~TStagA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagA2AVectors,
                    ARG(TStagA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),
                    MSolver);

/******************************************************************************
 *                       TStagA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagA2AVectors<FImpl, Pack>::TStagA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;
    
    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver);
    in.push_back(par().noise);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(ColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    int         Ls          = env().getObjectLs(par().action);
    
    
    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1,
              2*Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1,
              2*Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagA2AVectors<FImpl, Pack>::execute(void)
{
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver);
    auto        &noise     = envGet(ColorDiagonalNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);
    double      mass       = par().mass;
    
    envGetTmp(A2A, a2a);
    
    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
        << " using eigenpack '" << par().eigenPack << "' ("
        << 2*Nl_ << " low modes) and noise '"
        << par().noise << "' (" << noise.fermSize()
        << " noise vectors)" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
        << " using noise '" << par().noise << "' (" << noise.size()
        << " noise vectors)" << std::endl;
    }
    // Low modes
    auto &epack  = envGet(Pack, par().eigenPack);
    for (unsigned int il = 0; il < Nl_; il++)
    {
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il]-mass*mass));
        
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[2*il], epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeV(v[2*il+1], epack.evec[il], eval, 1);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[2*il], f5, epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeV5D(v[2*il+1], f5, epack.evec[il], eval, 1);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[2*il], epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeW(w[2*il+1], epack.evec[il], eval, 1);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[2*il], f5, epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeW5D(w[2*il+1], f5, epack.evec[il], eval, 1);
        }
        stopTimer("W low mode");
    }
    
    // High modes
#if 1
    
    FermionField sub(env().getGrid());
    FermionField noise_ih(env().getGrid());
    
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << 2*Nl_ + ih
        << " (" << ((Nl_ > 0) ? "high " : "")
        << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[2*Nl_ + ih], noise.getFerm(ih));
            // subtract the low modes
            sub = Zero();
            noise_ih = noise.getFerm(ih);
            for (int i=0;i<2*Nl_;i++) {
              const FermionField& tmp = v[i];
              // eval of unpreconditioned Dirac op
              std::complex<double> eval(mass,sqrt(epack.eval[i/2]-mass*mass));
              eval = i%2 ? eval : conjugate(eval);
              // need to subtract |l><l|/lambda_l |src>, so
              // mult, do not / by eval since v already has 1/eval
              axpy(sub,TensorRemove(innerProduct(tmp,noise_ih)) * eval,tmp,sub);
            }
            v[2*Nl_ + ih] -= sub;
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[2*Nl_ + ih], f5, noise.getFerm(ih));
        }
        std::cout << "v high " << ih << std::endl;
        //std::cout << v[2*Nl_+ih] << std::endl;
        
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << 2*Nl_ + ih
        << " (" << ((Nl_ > 0) ? "high " : "")
        << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[2*Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[2*Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }
#endif
    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
    //printMem("End StagA2AVectors execute() ", env().getGrid()->ThisRank());
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
