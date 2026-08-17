/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/StagA2AMesonFieldOpt.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_StagA2AMesonFieldOpt_hpp_
#define Hadrons_MContraction_StagA2AMesonFieldOpt_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/algorithms/blas/A2ASpatialSum.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  All-to-all staggered meson field creation, A2ASpatialSum-based.
 *
 *  Same physics/output as StagSparseA2AMesonField (single mu, zero momentum,
 *  one output file), but replaces StagMesonField's scalar thread_for
 *  reduction + single monolithic GlobalSumVector (A2Autils.h, StagMesonField
 *  -- commented there as taking up to 50% of time at 16 nodes) with a
 *  batched GEMM (A2ASpatialSum::PackLeftConj/PackRight + SumCacheBlocked),
 *  which decouples the GlobalSumVector message size from the block size via
 *  cacheBlock instead of doing one big reduction per kernel call.
 *
 *  No GammaRight step and no momentum phase loop: staggered mu is a choice
 *  of pre-built A2A vector set (selected upstream of this module), not a
 *  Dirac gamma multiplication applied here, and this module -- like
 *  StagSparseA2AMesonField -- only ever computes the zero-momentum field.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class StagA2AMesonFieldOptPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagA2AMesonFieldOptPar,
                                    int, cacheBlock,
                                    int, block,
                                    int, mu,
                                    std::string, left,
                                    std::string, right,
                                    std::string, output);
};

class StagA2AMesonFieldOptMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagA2AMesonFieldOptMetadata,
                                    std::string,momstr,
                                    std::string,gamstr);
};

template <typename FImpl>
class TStagA2AMesonFieldOpt : public Module<StagA2AMesonFieldOptPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TStagA2AMesonFieldOpt(const std::string name);
    // destructor
    virtual ~TStagA2AMesonFieldOpt(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    A2ASpatialSum<typename FImpl::SiteSpinor> spatial_sum_;
};

MODULE_REGISTER(StagA2AMesonFieldOpt, ARG(TStagA2AMesonFieldOpt<STAGIMPL>), MContraction);

/******************************************************************************
*                  TStagA2AMesonFieldOpt implementation                      *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagA2AMesonFieldOpt<FImpl>::TStagA2AMesonFieldOpt(const std::string name)
: Module<StagA2AMesonFieldOptPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagA2AMesonFieldOpt<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    return in;
}

template <typename FImpl>
std::vector<std::string> TStagA2AMesonFieldOpt<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagA2AMesonFieldOpt<FImpl>::setup(void)
{
    // A2ASpatialSum owns its own device-side buffers directly (no Hadrons
    // env allocation needed, unlike A2AMatrixBlockComputation's mCache_/mBuf_).
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagA2AMesonFieldOpt<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = envGetGrid(FermionField);

    int nt    = env().getDim().back();
    int N_i   = left.size();
    int N_j   = right.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    LOG(Message) << "Computing all-to-all Sparse meson fields (Opt)" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE))
                 << "/momentum/bilinear)" << std::endl;

    auto ionameFn = [this](const unsigned int m, const unsigned int g)
    {
        std::stringstream ss;

        ss << "mf000_mu"+std::to_string(par().mu);

        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const unsigned int m, const unsigned int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory())
               + "/" + ionameFn(1,1) + ".h5";
    };

    auto metadataFn = [this](const unsigned int m, const unsigned int g)
    {
        StagA2AMesonFieldOptMetadata md;

        md.momstr= "0 0 0";
        md.gamstr = "mu" + std::to_string(par().mu);

        return md;
    };

    // Every rank creates the output directory itself -- see the identical
    // comment in A2AMesonFieldOpt.hpp for why (parallel filesystem /
    // node-local storage safety).
    std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());
    Hadrons::mkdir(dirBase);
    grid->Barrier();

    unsigned int myRank = grid->ThisRank();

    // Single output file (one mu, zero momentum): only rank 0 creates it
    // and writes to it, unlike A2AMesonFieldOpt where next_*nstr_ > 1 lets
    // writes spread across ranks. GlobalSumVector (inside SumCacheBlocked)
    // is still collective -- every rank must call it -- only the IO is
    // rank-0-only.
    if (myRank == 0)
    {
        A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filenameFn(0, 0), ionameFn(0, 0), nt, N_i, N_j);
        io.initFile(metadataFn(0, 0), block);
    }
    grid->Barrier();

    // Output buffer: one (nt, Nii, Njj) block at a time.
    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

    // Pre-allocated result buffer, reused across all blocks -- see the
    // identical comment in A2AMesonFieldOpt.hpp on why RowMajor here.
    Eigen::Tensor<ComplexD, 3, Eigen::RowMajor> result(nt, block, block);

    startTimer("Allocate");
    spatial_sum_.AllocateRight(block, grid);
    spatial_sum_.AllocateLeft(block);
    stopTimer("Allocate");

    double                fillTime   = 0.;
    double                writeTime  = 0.;
    std::array<double, 7> ioTimings  = {};
    std::array<double, 5> sumTimings = {};
    std::array<double, 5> sumBytes   = {};

    for (int jb = 0; jb < N_j; jb += block)
    {
        int Njj = std::min(N_j - jb, block);

        startTimer("Allocate");
        spatial_sum_.AllocateRight(Njj, grid);
        stopTimer("Allocate");

        // No GammaRight step: staggered has no Dirac gamma structure here,
        // pack the right vectors directly (unlike A2AMesonFieldOpt).
        startTimer("Pack vectors");
        spatial_sum_.PackRight(right, jb, Njj);
        stopTimer("Pack vectors");

        for (int ib = 0; ib < N_i; ib += block)
        {
            int Nii = std::min(N_i - ib, block);

            startTimer("Allocate");
            spatial_sum_.AllocateLeft(Nii);
            stopTimer("Allocate");

            startTimer("Pack vectors");
            spatial_sum_.PackLeftConj(left, ib, Nii);
            stopTimer("Pack vectors");

            // No momentum phase: this module only ever computes p=0, so
            // there's nothing analogous to A2AMesonFieldOpt's ApplyPhaseRight.
            startTimer("Sum");
            spatial_sum_.SumCacheBlocked(result, cacheBlock, &sumTimings, &sumBytes);
            stopTimer("Sum");

            startTimer("IO");
            double dt = -usecond();
            A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);
            thread_for_collapse(3, t, nt, {
                for (int ii = 0; ii < Nii; ii++)
                for (int jj = 0; jj < Njj; jj++)
                    mf(0, 0, (int)t, ii, jj) = result((int)t, ii, jj);
            });
            dt += usecond();
            fillTime += dt;

            if (myRank == 0)
            {
                double wt = -usecond();
                A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filenameFn(0, 0), ionameFn(0, 0), nt, N_i, N_j);
                io.saveBlock(mf, 0, 0, ib, jb, &ioTimings);
                wt += usecond();
                writeTime += wt;
            }
            stopTimer("IO");
        } // ib
    } // jb

    // Throughput of the post-GEMM Sum() stages -- see the identical comment
    // in A2AMesonFieldOpt.hpp. GlobalSumVector here is the direct
    // replacement for StagMesonField's t_gsum, which StagSparseA2AMesonField
    // computes but never surfaces to the log.
    auto gbps = [](double bytes, double us)
    {
        return (us > 0.) ? bytes / us * 1.e6 / 1024. / 1024. / 1024. : 0.;
    };
    LOG(Message) << "Sum detail (us), rank " << myRank << ":" << std::endl;
    LOG(Message) << "  GEMM            = " << sumTimings[0] << std::endl;
    LOG(Message) << "  device->host    = " << sumTimings[1]
                 << " (" << gbps(sumBytes[1], sumTimings[1]) << " GB/s)" << std::endl;
    LOG(Message) << "  transpose-1     = " << sumTimings[2]
                 << " (" << gbps(sumBytes[2], sumTimings[2]) << " GB/s)" << std::endl;
    LOG(Message) << "  GlobalSumVector = " << sumTimings[3]
                 << " (" << gbps(sumBytes[3], sumTimings[3]) << " GB/s)" << std::endl;
    LOG(Message) << "  transpose-2     = " << sumTimings[4]
                 << " (" << gbps(sumBytes[4], sumTimings[4]) << " GB/s)" << std::endl;
    LOG(Message) << "IO detail (us), rank " << myRank << ":" << std::endl;
    LOG(Message) << "  fill            = " << fillTime  << std::endl;
    LOG(Message) << "  write (rank 0)  = " << writeTime << std::endl;
    if (myRank == 0)
    {
        LOG(Message) << "  open            = " << ioTimings[0]  << std::endl;
        LOG(Message) << "  push/group      = " << ioTimings[1]  << std::endl;
        LOG(Message) << "  openDataSet     = " << ioTimings[2]  << std::endl;
        LOG(Message) << "  getSpace        = " << ioTimings[3]  << std::endl;
        LOG(Message) << "  selectHyperslab = " << ioTimings[4]  << std::endl;
        LOG(Message) << "  write           = " << ioTimings[5]  << std::endl;
        LOG(Message) << "  close(fsync)    = " << ioTimings[6]  << std::endl;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StagA2AMesonFieldOpt_hpp_
