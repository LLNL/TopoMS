/*
Copyright (c) 2018, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Attila G Gyulassy
(jediati@sci.utah.edu).
LLNL-CODE-745278. All rights reserved.

This file is part of TopoMS, Version 1.0. For details, see
https://github.com/LLNL/TopoMS. Please also read this link – Additional BSD
Notice.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

• Redistributions of source code must retain the above copyright notice, this
list of conditions and the disclaimer below.
• Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the disclaimer (as noted below) in the
documentation and/or other materials provided with the distribution.
• Neither the name of the LLNS/LLNL nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE
U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
Department of Energy (DOE). This work was produced at Lawrence Livermore
National Laboratory under Contract No.  DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
Security, LLC nor any of their employees, makes any warranty, express or
implied, or assumes any liability or responsibility for the accuracy,
completeness, or usefulness of any information, apparatus, product, or process
disclosed, or represents that its use would not infringe privately-owned
rights.

3. Also, reference herein to any specific commercial products, process, or
services by trade name, trademark, manufacturer or otherwise does not
necessarily constitute or imply its endorsement, recommendation, or favoring by
the United States Government or Lawrence Livermore National Security, LLC. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or Lawrence Livermore National
Security, LLC, and shall not be used for advertising or product endorsement
purposes.
*/

/**
 *  @file    TopoMS.h
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *  @version 1.0
 *
 *  @brief This file provides the core functionality for TopoMS
 *
 *  @section DESCRIPTION
 *
 *  This file provides the core functionality for TopoMS
 *
 */

#ifndef _TOPOMS_H_
#define _TOPOMS_H_

#include <string>
#include <vector>
#include <unordered_set>

#include "ConfigParser.h"

#include "Vec3.h"
#include "MolecularSystem.h"
#include "Material.h"
#include "InputFormats.h"

/// -----------------------------------------------------------------------
/// Forward declaration of the required components of MSC library
/// -----------------------------------------------------------------------
#include "basic_types.h"
#include "vectors.h"

namespace MSC {

    class kdtree;
    class RegularGrid;
    class RegularGridTrilinearFunction;
    class TopologicalRegularGridRestricted;
    class TopologicalRegularMaskedRestrictedGrid;
    class DiscreteGradientLabeling;
    class TopologicalExplicitDenseMeshFunction;
    class TopologicalRegularGrid;
    class IndexCompareLessThan;
    class TerminateNearExtrema;

    template<typename LABEL_TYPE>
    class DenseLabeling;

    template<typename INPUT_LABEL_TYPE>
    class VertexLabelingToBoundaryLabeling;

    template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
    class MorseSmaleComplexBasic;

    template< class Advector, class Comparer>
    class NumericIntegratorExpandingRegionStopFiltered2;

    template< int StepMultiplier>
    class AdaptiveEulerAdvector;

    template< class Advector, class Comparer>
    class NumericIntegratorRegionStop;

    template< class Advector, class Comparer>
    class NumericIntegratorExpandingRegionStop;

    template< class Advector>
    class NumericStreamlineIntegrator;

    template< class Comparer>
    class IsolatedRegionRemover;

    template< class Comparer>
    class IsolatedRegionRemoverMasked;

    template<class MSCType>
    class MSCSelectorLivingNodes;

    template<class MSCType>
    class MSCSelectorNodeIndex;

    template<class MSCType>
    class MSCSelectorRepresentative1Saddle;

    template <class Mesh, class MeshFunction>
    class RobinsLabelingAlgorithm;

    template <class GridType, class MeshFunction, class LabelingType>
    class TopologicalGradientUsingAlgorithms;

    typedef IndexCompareLessThan Comparer;
    typedef TopologicalRegularMaskedRestrictedGrid GridType;

    typedef MorseSmaleComplexBasic<FLOATTYPE, GridType, TopologicalExplicitDenseMeshFunction, DiscreteGradientLabeling> MSCType;
    typedef NumericIntegratorExpandingRegionStopFiltered2<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType;
    typedef NumericIntegratorRegionStop<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType2;
    typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
    typedef IsolatedRegionRemoverMasked<Comparer> RegionRemoverType;
}

/// -----------------------------------------------------------------------
/// Forward declaration done!
/// -----------------------------------------------------------------------

/**
  *  @brief This class provides the core functionality for TopoMS
  */
class TopoMS {

private:

    // ----------------------------------------------------------------------
    // input related
    std::string m_datadir;
    FLOATTYPE *m_func;

    // bader-related variables
    std::vector<FLOATTYPE> chg_atoms, vol_atoms;
    std::vector<FLOATTYPE> chg_extrema, vol_extrema;
    std::vector<int> labels_atoms, labels_extrema;

    std::map<size_t, size_t> extrema2atoms;

    // ----------------------------------------------------------------------
    // msc-related variables
    MSC::RegularGrid* m_grid;
    MSC::GridType *m_tgrid;

    MSC::RegularGridTrilinearFunction* m_gridfunc;
    MSC::TopologicalExplicitDenseMeshFunction* m_topofunc;

    MSC::IntegratorType* m_integrator;
    MSC::RegionRemoverType* m_integrator2;
    MSC::MSCType* m_msc;
    MSC::kdtree* m_kdtree_minima;
    MSC::kdtree* m_kdtree_atoms;

    float persistence_val;
    float filter_val;

public:

    // TODO: make this private
    enum INPUT_TYPE {IT_UNKNOWN = 0, IT_VASP = 1, IT_CUBE = 2};
    INPUT_TYPE m_inputtype;

    Config *m_config;
    MS::SystemInfo m_metadata;
    std::vector<Material> m_atoms;

    std::vector<std::vector<MSC::Vec3d> > m_paths;
    bool m_negated;
    std::vector<INT_TYPE> m_nodes;

    // ----------------------------------------------------------------------
    // basic functions
    TopoMS() : m_config(0), m_inputtype(IT_UNKNOWN), m_negated(false),
                   m_grid(0), m_tgrid(0), m_func(0), m_gridfunc(0), m_topofunc(0),
                   m_integrator(0), m_integrator2(0), m_msc(0),
                   m_kdtree_minima(0), m_kdtree_atoms(0),
                    persistence_val(1.0), filter_val(1.0) {
        std::cout << "\n TopoMS v1.0\n\n";
    }

    bool load(const std::string &configfilename);
    bool init();

    bool kdtree_add(MSC::kdtree* kt, double x, double y, double z, int data) const;
    int kdtree_query(MSC::kdtree* kt, const double pos[3]) const;

    // main analysis
    bool bader();
    bool msc();
    void simplify_msc(FLOATTYPE pvalue, FLOATTYPE fvalue);

    // basic queries
    FLOATTYPE* get_func() {                     return m_func;              }
    const FLOATTYPE* get_func() const {         return m_func;              }

    const size_t* get_gridDims() const {    return m_config->grid_dims; }
    size_t get_gridSize() const {           return m_config->grid_dims[0]*m_config->grid_dims[1]*m_config->grid_dims[2];    }

    std::pair<FLOATTYPE, FLOATTYPE> get_frange() const {

        FLOATTYPE minval = *std::min_element(m_func, m_func+get_gridSize());
        FLOATTYPE maxval = *std::max_element(m_func, m_func+get_gridSize());

        return (this->m_negated) ? std::make_pair(-1*maxval, -1*minval) :
                                   std::make_pair(minval, maxval);

    }
    FLOATTYPE get_fsum() const {
        FLOATTYPE s = std::accumulate(m_func, m_func+get_gridSize(), 0.0);
        return (this->m_negated) ? -1*s : s;
    }
    FLOATTYPE get_fint() const {
        return get_fsum() / FLOATTYPE(get_gridSize());
    }

    float get_persistence() const {         return persistence_val;    }
    float set_persistence(float _) {        persistence_val = _;   }

    float get_filterval() const {           return filter_val;          }
    float set_filterval(float _)  {         filter_val = _;          }

    float pers_val2perc(float val) const {
        static std::pair<FLOATTYPE, FLOATTYPE> vr = this->get_frange();
        return (100.0*val / (vr.second - vr.first));
    }
    float pers_perc2val(float val) {
        static std::pair<FLOATTYPE, FLOATTYPE> vr = this->get_frange();
        return (val/100.0 * (vr.second - vr.first) + vr.first);
    }

    // query msc
    const std::vector<long long> get_extrema() const;
    void cellid2Coords(const INDEX_TYPE &cellIdx, MSC::Vec3l &ncoords) const;

    bool msc_available() const;
    size_t get_msc_nnodes() const;
    bool get_msc_isNodeAlive(size_t idx) const;
    bool get_msc_isNodeFiltered(size_t idx) const;
    bool get_msc_node(size_t idx, int &dim, INDEX_TYPE &cellIdx, MSC::Vec3l &ncoords) const;

    void print_node(size_t idx) const;

    // ----------------------------------------------------------------------
    // I/O functions

    void write_bader() const {
        //write_bader_max2vol(m_datadir+"Maxima2Vol.vti");
        //write_bader_max2chgvol(m_datadir+"Maxima2ChgVol.dat");

        write_bader_atoms2vol(m_datadir+"Atoms2Vol.vti");
        write_bader_atoms2chgvol(m_datadir+"Atoms2ChgVol.dat");

        return;

        if (this->m_inputtype == IT_VASP) {
            MS::VASP::write_CHGCAR<int>(m_datadir+"Atoms2Vol", this->m_metadata, this->m_atoms, labels_atoms.data());
        }
    }

    void write_msc() const {

        char fpers[10];
        sprintf(fpers, "%.4f", get_persistence());
        write_mgraph(m_datadir+"MolecularGraph_"+std::string(fpers)+".vtp");
    }

private:

    void write_bader_max2vol(const std::string &filename) const;
    void write_bader_atoms2vol(const std::string &filename) const;
    void write_bader_atoms2chgvol(const std::string &filename) const;
    void write_bader_max2chgvol(const std::string &filename, bool filter_small = true) const;

    void write_mgraph(const std::string &filename) const;

    /// ----------------------------------------------------------------------
};

#endif
