/*
 * LoadA2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Vaishakhi Moningi <vaishu.moningi@gmail.com>
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
#ifndef Hadrons_MIO_LoadA2AVectors_hpp_
#define Hadrons_MIO_LoadA2AVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadA2AVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadA2AVectorsPar,
                                    std::string,              filestem,
                                    std::vector<std::string>, filestems,
                                    bool,                     multiFile,
                                    unsigned int,             size,
                                    int,                      inc,
                                    int,                      tinc);
    // inc/tinc: sparse blocking factors; 0 or 1 means use the full (fine)
    // lattice grid.  Set inc>1 to load vectors living on a coarse grid
    // (written by StagSparseA2AVectorsGridIo with the same inc/tinc).
    // filestems: if non-empty, load each stem as a chunk of 'size' vectors
    // and concatenate them in order. 'filestem' is ignored in that case.
};

template <typename FImpl>
class TLoadA2AVectors: public Module<LoadA2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadA2AVectors(const std::string name);
    // destructor
    virtual ~TLoadA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadA2AVectors,     TLoadA2AVectors<FIMPL>,    MIO);
MODULE_REGISTER_TMP(StagLoadA2AVectors, TLoadA2AVectors<STAGIMPL>, MIO);


/******************************************************************************
 *                      TLoadA2AVectors implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadA2AVectors<FImpl>::TLoadA2AVectors(const std::string name)
: Module<LoadA2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadA2AVectors<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadA2AVectors<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadA2AVectors<FImpl>::setup(void)
{
    int inc  = (par().inc  <= 0) ? 1 : par().inc;
    int tinc = (par().tinc <= 0) ? 1 : par().tinc;
    int nchunks = par().filestems.empty() ? 1 : (int)par().filestems.size();

    GridBase *grid;
    if (inc > 1 || tinc > 1)
    {
        std::vector<int> bs = {inc, inc, inc, tinc};
        grid = envGetCoarseGrid(FermionField, bs);
    }
    else
        grid = envGetGrid(FermionField);

    envCreate(std::vector<FermionField>, getName(), 1,
              nchunks * par().size, grid);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadA2AVectors<FImpl>::execute(void)
{
    auto      &vec  = envGet(std::vector<FermionField>, getName());
    const int  traj = vm().getTrajectory();

    if (par().filestems.empty())
    {
        // single file — original behaviour
        A2AVectorsIo::read(vec, par().filestem, par().multiFile, traj);
    }
    else
    {
        // multi-chunk: load into a temporary vector sized per chunk so that
        // the record.index check in A2AVectorsIo::read (which expects indices
        // 0..size-1) stays satisfied, then copy into the correct slice of vec.
        int        chunkSz = (int)par().size;
        GridBase  *grid    = vec[0].Grid();

        for (int c = 0; c < (int)par().filestems.size(); c++)
        {
            LOG(Message) << "Loading chunk " << c
                         << " from " << par().filestems[c] << std::endl;
            std::vector<FermionField> tmp(chunkSz, grid);
            A2AVectorsIo::read(tmp, par().filestems[c], par().multiFile, traj);
            for (int i = 0; i < chunkSz; i++)
                vec[c * chunkSz + i] = tmp[i];
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadA2AVectors_hpp_
