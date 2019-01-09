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

#include "version.h"
#include "InputFormats.h"
#include "ConfigParser.h"

#include "Vec3.h"
#include "Material.h"
#include "MSCBond.h"
#include "MolecularSystem.h"

#include "basic_types.h"
#include "vectors.h"

/// -----------------------------------------------------------------------
/// Forward declaration of the required components of MSC library
/// -----------------------------------------------------------------------

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

    template<int StepMultiplier>
    class AdaptiveEulerAdvector;

    template<class Advector, class Comparer>
    class NumericIntegratorRegionStop;

    template<class Advector, class Comparer>
    class NumericIntegratorExpandingRegionStop;

    template<class Advector>
    class NumericStreamlineIntegrator;

    template<class Comparer>
    class IsolatedRegionRemover;

    template<class Comparer>
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

#ifdef USE_VTK
    class vtkImageData;
    class vtkVolumeSlicer;
#endif

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

    // the function used for msc analysis
    FLOATTYPE *m_mscfunc;

    // the function used for bader anlaysis
        // if ref charge is provided
            // m_mscfunc = ref charge
            // m_baderfunc = main charge
    FLOATTYPE *m_baderfunc;

#ifdef USE_VTK
    vtkImageData *m_vtkFunction;
    vtkImageData *m_vtkVolLabeling;
#endif

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
    MSC::DenseLabeling<int> *m_volumelabeling;
    MSC::MSCType* m_msc;
    MSC::kdtree* m_kdtree_minima;
    MSC::kdtree* m_kdtree_atoms;

    FLOATTYPE persistence_val;
    FLOATTYPE filter_val;

    FLOATTYPE vacthreshold_in_fileUnits;

    // ----------------------------------------------------------------------
public:

    // TODO: make this private
    enum INPUT_TYPE {IT_UNKNOWN = 0, IT_VASP = 1, IT_CUBE = 2};
    enum FIELD_TYPE {FT_UNKNOWN = 0, FT_CHG = 1, FT_POT = 2};

    INPUT_TYPE m_inputtype;
    FIELD_TYPE m_fieldtype;

    Config *m_config;
    MS::SystemInfo m_metadata;
    std::vector<Material> m_atoms;

    bool m_negated;
    bool slice_labels;

#ifdef USE_VTK
    vtkVolumeSlicer *m_slicer_function;
    vtkVolumeSlicer *m_slicer_label;

    vtkVolumeSlicer* slicer() { return (slice_labels) ? m_slicer_label : m_slicer_function; }
    bool slice_color_by_label() const { return slice_labels;                                }
#endif


    std::map<INT_TYPE, MSCBond> m_mscbonds;

    // ----------------------------------------------------------------------
    // initialization functions

    TopoMS() : m_config(0), m_inputtype(IT_UNKNOWN), m_fieldtype(FT_UNKNOWN), m_negated(false),
                   m_grid(0), m_tgrid(0), m_mscfunc(0), m_gridfunc(0), m_topofunc(0),
                   m_integrator(0), m_integrator2(0), m_msc(0),
                   m_kdtree_minima(0), m_kdtree_atoms(0),
                   persistence_val(1.0), filter_val(1.0)
#ifdef USE_VTK
                   , m_slicer_function(0), m_slicer_label(0), slice_labels(0)
#endif
    {
        std::cout << "\n TopoMS v" << TOPOMS_VERSION_STR << " [released on " << TOPOMS_VERSION_RELEASE_DATE << "]\n\n";
    }

    bool load(const std::string &configfilename, std::string datafile = "");
    bool init();

    // ----------------------------------------------------------------------
    // basic interface

    FLOATTYPE* get_func() {                     return m_mscfunc;              }
    const FLOATTYPE* get_func() const {         return m_mscfunc;              }

    size_t get_gridSize() const {               return m_metadata.grid_sz();    }

    FLOATTYPE get_persistence() const {         return persistence_val;    }
    FLOATTYPE get_filterval() const {           return filter_val;         }

    FLOATTYPE set_persistence(float _) {        persistence_val = _;       }
    FLOATTYPE set_filterval(float _)  {         filter_val = _;            }

    FLOATTYPE pers_val2perc(FLOATTYPE val) const {
        static std::pair<FLOATTYPE, FLOATTYPE> vr = this->get_frange();
        return (100.0*val / (vr.second - vr.first));
    }
    FLOATTYPE pers_perc2val(float val) {
        static std::pair<FLOATTYPE, FLOATTYPE> vr = this->get_frange();
        return (val/100.0 * (vr.second - vr.first) + vr.first);
    }

    float get_periodic_cutoff(bool world=false) const {

        float cut[3] = {float(m_metadata.m_grid_dims[0]), float(m_metadata.m_grid_dims[1]), float(m_metadata.m_grid_dims[2])};
        if (world) {
            float gcut[3] = {cut[0],cut[1],cut[2]};
            m_metadata.grid_to_world(gcut, cut);
        }
        return 0.3*std::min(cut[0], std::min(cut[1], cut[2]));
    }

    // ----------------------------------------------------------------------
    // show_raw:    ignore negated!
    std::pair<FLOATTYPE, FLOATTYPE> get_frange(const FLOATTYPE *func, bool show_raw) const {

        const bool negate = !show_raw && this->m_negated;
        const size_t sz = get_gridSize();

        FLOATTYPE minval = *std::min_element(func, func+sz);
        FLOATTYPE maxval = *std::max_element(func, func+sz);
        return (negate) ? std::make_pair(-1*maxval, -1*minval) : std::make_pair(minval, maxval);
    }

    std::pair<FLOATTYPE, FLOATTYPE> get_frange() const {
        return get_frange(m_mscfunc, false);
    }

    // ----------------------------------------------------------------------
    // analysis interface
    bool bader();
    bool msc();

    // ----------------------------------------------------------------------
    // bader and msc analysis

    void extract_mgraph(FLOATTYPE pvalue, FLOATTYPE fvalue);
    void extract_lpot_nbrhood(FLOATTYPE pvalue, FLOATTYPE fvalue, float gpos[], unsigned cp_idx);
    void extract_lpot_nbrhood_li(FLOATTYPE pvalue, FLOATTYPE fvalue);

    void analyze_bonds();

    // ----------------------------------------------------------------------
    // interface with bader analysis
    // ----------------------------------------------------------------------
    const std::vector<INDEX_TYPE>& bader_get_extrema() const;

    int bader_get_atomLabeling(INDEX_TYPE vIdx) const;
    int bader_get_atomLabeling(size_t x, size_t y, size_t) const;
    int bader_get_atomLabeling(double pos[3]) const;

    // ----------------------------------------------------------------------
    // interface with the msc library
    // ----------------------------------------------------------------------
    bool msc_is_available() const;
    size_t msc_get_nnodes() const;
    size_t msc_get_narcs() const;

    bool msc_is_nodeAlive(size_t idx) const;
    bool msc_is_arcAlive(size_t idx) const;
    bool msc_is_nodeFiltered(size_t idx) const;

    int msc_get_nodeDim(size_t idx) const;
    int msc_get_arcDim(size_t idx) const;

    bool msc_get_node(size_t idx, int &dim, INDEX_TYPE &cellIdx, MSC::Vec3d &ncoords) const;
    void msc_get_arcGeometry(size_t idx, std::vector<INDEX_TYPE> &arc_geom) const;

    void msc_cellid_to_gcoords(const INDEX_TYPE &cellIdx, MSC::Vec3d &ncoords) const;
    void msc_cellid_to_wcoords(const INDEX_TYPE &cellIdx, float wcoords[]) const;
    void msc_gcoords_to_wcoords(const MSC::Vec3d &ncoords, float wcoords[]) const;
    void msc_print_node(size_t idx) const;

    bool refresh_orthogonalSlice(int saddleNodeId, int param=0, int minp=-100, int maxp=100);

private:
    // ----------------------------------------------------------------------
    // kd tree handling
    // ----------------------------------------------------------------------
    bool kdtree_add(MSC::kdtree* kt, double x, double y, double z, int data) const;
    int kdtree_query(MSC::kdtree* kt, const double pos[3]) const;

    // ----------------------------------------------------------------------
    // vtk functionalities
    // ----------------------------------------------------------------------

#ifdef USE_VTK
    vtkImageData *create_vtkImagedata(const std::string &fname) const;
    vtkImageData* create_vtkImagedata(const double *volume, const std::string &fname) const;

    static void compute_slice(const MSC::Vec3d &origin, const std::vector<MSC::Vec3d> &nbrs, vtkVolumeSlicer *slicer);

    bool filter_slice(vtkVolumeSlicer *slicer, const std::vector<size_t> &atomids, bool overwrite = false) const;
    std::pair<double, double> integrate_slice(const vtkVolumeSlicer *slicer) const;
#endif

    // ----------------------------------------------------------------------
    // I/O functions
    // ----------------------------------------------------------------------

    void write_labeling_vti(const std::string &filename, const std::string &fieldname, const std::vector<int> &labels) const;
    void write_mgraph_vtp(const std::string &filename) const;

    void write_bader_atoms2chgvol(const std::string &filename) const;
    void write_bader_max2chgvol(const std::string &filename, bool filter_small = true) const;
    void write_msc_bond_stats(const std::string &filename) const;

public:

    void write_bader() const {

        write_labeling_vti(m_config->infilename+"-Atoms2Vol.vti", "atom_labeling", this->labels_atoms);
        write_bader_atoms2chgvol(m_config->infilename+"-Atoms2ChgVol.dat");
        return;
        if (this->m_inputtype == IT_VASP)
            MS::VASP::write_CAR<int>(m_config->infilename+"-Atoms2Vol", this->m_metadata, this->m_atoms, labels_atoms.data());
    }

    void write_msc() const {

        char fpers[10], ffilt[10];
        sprintf(fpers, "%.4f", get_persistence());
        sprintf(ffilt, "%.4f", get_filterval());
        write_mgraph_vtp(m_config->infilename+"-MolecularGraph_"+std::string(fpers)+"_"+std::string(ffilt)+".vtp");
        write_msc_bond_stats(m_config->infilename+"-MolecularBondStats_"+std::string(fpers)+"_"+std::string(ffilt)+".txt");
    }

    // ----------------------------------------------------------------------
    // ----------------------------------------------------------------------
};

#endif
