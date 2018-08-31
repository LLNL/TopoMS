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
 *  @file    TopoMS.cpp
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

//#define OUTPUT_DEBUG_FIELDS

// -----------------------------------------------------------------------
// TopoMS headers
#include "ConfigParser.h"
#include "InputFormats.h"
#include "MolecularSystem.h"
#include "TopoMS.h"

// -----------------------------------------------------------------------
// MSC headers
#include "strictly_numeric_integrator.h"
#include "numeric_integrator_region_stop.h"
#include "numeric_integrator_expanding_region_stop_filtered2.h"
#include "timing.h"
#include "labeling_to_bounary_labeling.h"
#include "topological_explicit_mesh_function.h"
#include "topological_region_growing_simple_gradient_builder.h"
#include "topological_convergent_gradient_builder.h"
#include "robin_labeling.h"
#include "adaptive_in_quad_euler_advector.h"
#include "numeric_integrator_2d_restricted_expanding_region.h"
#include "topological_2d_restricted_expanding_regions.h"
#include "topological_gradient_using_algorithms.h"
#include "topological_regular_grid_restricted.h"
#include "topological_regular_masked_restricted_grid.h"
#include "isolated_region_remover.h"
#include "isolated_region_remover2.h"
#include "topological_utility_functions.h"
#include "morse_smale_complex_basic.h"
#include "morse_smale_complex_restricted.h"
#include "kdtree.h"
#include "msc_selectors.h"
#include "numeric_streamline_integrator.h"

#include <functional>

// -----------------------------------------------------------------------------------------
// these wrappers are needed so rest of the code remains independent of the MSC headers
// -----------------------------------------------------------------------------------------

bool TopoMS::msc_available() const {
    return this->m_msc != 0;
}

void TopoMS::cellid2Coords(const INDEX_TYPE &cellIdx, MSC::Vec3l &ncoords) const {
    m_tgrid->cellid2Coords(cellIdx, ncoords);
}

size_t TopoMS::get_msc_nnodes() const {
    return this->m_msc->numNodes();
}

bool TopoMS::get_msc_isNodeAlive(size_t idx) const {
    return this->m_msc->isNodeAlive(idx);
}

bool TopoMS::get_msc_isNodeFiltered(size_t idx) const {
    MSC::node<FLOATTYPE> & n = this->m_msc->getNode(idx);
    return ( fabs(n.value) > fabs(this->filter_val) );
}

bool TopoMS::get_msc_node(size_t idx, int &dim, INDEX_TYPE &cellIdx, MSC::Vec3l &ncoords) const {
    MSC::node<FLOATTYPE> & n = this->m_msc->getNode(idx);
    dim = n.dim;
    if(this->m_negated){
        dim = 3-dim;
    }
    cellIdx = n.cellindex;
    m_tgrid->cellid2Coords(cellIdx, ncoords);
	return true;
}

const std::vector<INDEX_TYPE> TopoMS::get_extrema() const {
    if(m_integrator == 0)
        return std::vector<INDEX_TYPE>();
    return m_integrator->GetExtrema();
}

void TopoMS::print_node(size_t idx) const {

    const MSC::node<FLOATTYPE> &n = this->m_msc->getNode(idx);
    int dim = n.dim;
    FLOATTYPE val = m_topofunc->cellValue(n.cellindex);

    if(this->m_negated){
        dim = 3-dim;
        val = -val;
    }

    std::string type = (dim == 0) ? "min" :
                       (dim == 1) ? "1-sad" :
                       (dim == 2) ? "2-sad" :
                       (dim == 3) ? "max" : "unknown";

    MSC::Vec3l ncoords;
    m_tgrid->cellid2Coords(n.cellindex, ncoords);
    printf(" %s (node %d) at (%.1f %.1f %.1f) val = %f\n", type.c_str(), idx,
                        0.5*ncoords[0], 0.5*ncoords[1], 0.5*ncoords[2], val
            );
}

// -----------------------------------------------------------------------------------------

bool TopoMS::kdtree_add(MSC::kdtree* kt, double x, double y, double z, int data) const {

    static const INDEX_TYPE X = m_config->grid_dims[0];
    static const INDEX_TYPE Y = m_config->grid_dims[1];
    static const INDEX_TYPE Z = m_config->grid_dims[2];

    static const int sx = m_config->is_periodic[0] ? -1 : 0;
    static const int sy = m_config->is_periodic[1] ? -1 : 0;
    static const int sz = m_config->is_periodic[2] ? -1 : 0;

    static const int ex = m_config->is_periodic[0] ?  2 : 1;
    static const int ey = m_config->is_periodic[1] ?  2 : 1;
    static const int ez = m_config->is_periodic[2] ?  2 : 1;

    int* kdata = new int[1];
    kdata[0] = data;

    for (int px = sx; px < ex; px++)
    for (int py = sy; py < ey; py++)
    for (int pz = sz; pz < ez; pz++) {

        double kpos[3];
        kpos[0] = x + X * px;
        kpos[1] = y + Y * py;
        kpos[2] = z + Z * pz;

        kd_insert(kt, kpos, kdata);
    }
	return true;
}

int TopoMS::kdtree_query(MSC::kdtree *kt, const double pos[3]) const {

    MSC::kdres* res = kd_nearest(kt, pos);
    if (res == NULL) {
        //printf("TopoMS::kdtree_query(). Could not find position (%f %f %f)\n", pos[0], pos[1], pos[2]);
        return -1;
    }

    int* cc = (int*)kd_res_item(res, 0);
    return cc[0];
}

// -----------------------------------------------------------------------------------------
// load
// -----------------------------------------------------------------------------------------
bool TopoMS::load(const std::string &configfilename) {

    Utils::print_separator();
    std::cout << " =====>> Loading data using " << configfilename << "...\n";

    m_config = new Config( configfilename );
    m_config->parse();

    std::string extn = Utils::get_extn( m_config->infilename );
    m_datadir = Utils::get_directory(m_config->infilename);

    // -----------------------------------------------------------------------

    if (m_config->infiletype == "VASP") {

        m_func = MS::VASP::read_CHGCAR<FLOATTYPE>(m_config->infilename, m_metadata, m_atoms);
        this->m_inputtype = IT_VASP;
    }
    else if (m_config->infiletype == "CUBE") {

        m_func = MS::Cube::read<FLOATTYPE>(m_config->infilename, m_metadata, m_atoms);
        this->m_inputtype = IT_CUBE;
    }
    for(int i = 0; i < 3; i++)
        m_config->grid_dims[i] = m_metadata.m_grid_dims[i];

    m_metadata.print ();

    // ------------------------------------
    {
        const size_t nsz = m_metadata.grid_sz();

        const FLOATTYPE minval = *std::min_element(m_func, m_func+nsz);
        const FLOATTYPE maxval = *std::max_element(m_func, m_func+nsz);
        const FLOATTYPE sumval = std::accumulate(m_func, m_func+m_metadata.grid_sz(), 0.0);

        printf("\n num atoms = %d\n", m_atoms.size());
        printf(" function range [%f %f]\n", minval, maxval);

        // integration as a charge
        const FLOATTYPE chg_factor = (this->m_inputtype == IT_CUBE) ? m_metadata.volume() / (m_metadata.m_l_unit*m_metadata.m_l_unit*m_metadata.m_l_unit) / m_metadata.m_e_unit : 1.0;
        const FLOATTYPE totchg = sumval * chg_factor / (FLOATTYPE) m_metadata.grid_sz();

        printf(" integration as charge = %f\n", totchg);
        printf(" vol = %f, lunit = %f, eunit = %f : chg_factor = %f\n", m_metadata.volume(), m_metadata.m_l_unit, m_metadata.m_e_unit, chg_factor);
    }

    // ------------------------------------
    return true;
}

// use configuration file to initialize the program
bool TopoMS::init() {

    printf("\n\n =====>> Initializing TopoMS\n");

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    m_grid = new MSC::RegularGrid(
                MSC::Vec3i(m_config->grid_dims[0], m_config->grid_dims[1], m_config->grid_dims[2]),
                MSC::Vec3b(m_config->is_periodic[0], m_config->is_periodic[1], m_config->is_periodic[2]));

    m_gridfunc = new MSC::RegularGridTrilinearFunction(m_grid, m_func);
    m_gridfunc->ComputeGradFromImage(1);

    printf(" --- now negating the value!\n");
    m_negated = true;
    m_gridfunc->Negate();


    size_t nsz = m_config->grid_dims[0]*m_config->grid_dims[1]*m_config->grid_dims[2];
    const FLOATTYPE minval = *std::min_element(m_func, m_func+nsz);
    const FLOATTYPE maxval = *std::max_element(m_func, m_func+nsz);
    printf(" function range [%f %f]\n", minval, maxval);

    persistence_val = m_config->sval_threshold;
    filter_val = m_config->fval_threshold;

    printf("\n =====>> TopoMS initialized!");
    timer.EndGlobal ();
    timer.PrintAll ();
	return true;
}

/**
  *   @brief  Perform Bader analysis
  *
  *   @return success flag
  */

bool TopoMS::bader() {

    if(!m_config->bader)
        return false;

    // ---------------------------------------------------------------------
    // constant values required
    // ---------------------------------------------------------------------

    const size_t num_atoms = m_atoms.size();
    const size_t num_gridPts = m_metadata.grid_sz();

    const FLOATTYPE vol_box = m_metadata.volume();
    const FLOATTYPE vol_voxel = vol_box / (FLOATTYPE) num_gridPts;

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    // handle different types of units
        // l-unit and e-unit are set properly by file reading module

    const FLOATTYPE vol_Bohr2Angs = 1.0 / std::pow(m_metadata.m_l_unit, 3.0);
    const FLOATTYPE chg_hartree2e = 1.0 / m_metadata.m_e_unit;

    printf(" vol(box) = %E %s^3, %E Angs^3\n vol(voxel) = %E %s^3, %E Angs^3\n",
            vol_box, m_metadata.m_coordinate_unit.c_str(), vol_box*vol_Bohr2Angs,
            vol_voxel, m_metadata.m_coordinate_unit.c_str(), vol_voxel*vol_Bohr2Angs);

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    // the function value in VASP file is the number of electrons per voxel
        // it is the charge density multiplied by the total volume
        // therefore, to get the total number of electrons, multiply by the unit voxel volume ( / number of voxels)
            // charge is in e
            // voxel volume is in Angs

    // the function value in Cube file is the charge density per voxel
            // charge is in e or hartree
            // voxel volume is in unit grid volume (since the value is divided by total volume)
        // therefore, to get the total number of electrons, multiply by the voxel volume (total volume / number of voxels)

    // to convert from the charge density to actual charge (in e)
    const FLOATTYPE chgDens_fileUnit2e = (this->m_inputtype == IT_CUBE ? vol_box * chg_hartree2e * vol_Bohr2Angs : 1.0) / (FLOATTYPE) num_gridPts;

    //printf(" factor to convert chg density [file units] to chg [electrons] = %E\n", chgDens_fileUnit2e);

    // ---------------------------------------------------------------------
    // sum of function values
#ifdef USE_KAHAN_SUM
    const FLOATTYPE sum_func = (this->m_negated ? -1 : 1) * Utils::Kahan::sum(m_func, m_metadata.grid_sz());
#else
    const FLOATTYPE sum_func = (this->m_negated ? -1 : 1) * Utils::sum(m_func, m_metadata.grid_sz());
#endif
    const FLOATTYPE total_chg = sum_func * chgDens_fileUnit2e;

    // ---------------------------------------------------------------------
    // vacuum threshold

    // for fast thresholding, convert the vacuum threshold to the same unit as the file (so i can do a point-wise comparison later)
    const FLOATTYPE vacthreshold_in_fileUnits = m_config->vacuum_threshold *
                                                    (this->m_inputtype == IT_CUBE ? 1.0 / (chg_hartree2e * vol_Bohr2Angs) : vol_box);

    printf("\n");
    Utils::print_separator();
    printf(" =====>> Performing Bader Analysis\n\t total volume = %f Ang^3, total chg = %f e\n\t vacuum threshold = %E (in file units) = %E e per Angs^3\n\n",
                     vol_box*vol_Bohr2Angs, total_chg, vacthreshold_in_fileUnits, m_config->vacuum_threshold);

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    // ---------------------------------------------------------------------
    // volume decomposition
    // ---------------------------------------------------------------------

    const MSC::DenseLabeling<int> *volumelabeling = 0;
    {

        bool verbose = false;

        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Performing numeric integration for volume assignment");
        if (verbose) {
            printf("\n");
        }
        else {
            printf("...");
            fflush(stdout);
        }

        m_integrator = new MSC::IntegratorType(m_gridfunc, m_grid, m_config->error_threshold, m_config->grad_threshold, m_config->numiter);
        m_integrator->set_filter(vacthreshold_in_fileUnits);                         // the integrator will not seed a streamline for smaller values

        m_integrator->BeginIntegration(verbose);

#ifdef OUTPUT_DEBUG_FIELDS
        {
        std::string fname = m_datadir+"classes.raw";
        m_integrator->GetOutputLabels()->OutputToIntFile(fname.c_str());
        }
#endif

        // cleanup noise: do this only if you need msc as well
        if (this->m_config->msc ) {
            m_integrator2 = new MSC::RegionRemoverType(m_gridfunc, m_integrator->GetOutputLabels(), &(m_integrator->GetExtrema()));
            m_integrator2->ComputeOutput(verbose);
            volumelabeling = m_integrator2->GetOutputLabelsUnmasked();
        }
        else {
            volumelabeling = m_integrator->GetOutputLabels();
        }

#ifdef OUTPUT_DEBUG_FIELDS
        std::string fname = m_datadir+"classes_clean.raw";
        m_integrator2->GetOutputLabels()->OutputToIntFile(fname.c_str());
#endif

        if(verbose) {
            printf(" -- Done numeric integration of 3d volume!");
        }
        else {
            printf("Done!");
        }
        ltimer.EndGlobal ();
        ltimer.PrintAll ();
    }

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    const std::vector<INDEX_TYPE> &extrema = m_integrator->GetExtrema();

    size_t numextrema = extrema.size();
    INDEX_TYPE numlabels = volumelabeling->GetNumLabels();

    // ---------------------------------------------------------------------
    // atom--extrema mapping
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        // ---------------------------------------------------------------------
        printf(" -- Mapping %d extrema to %d atoms...", extrema.size(), num_atoms);
        fflush(stdout);

        m_kdtree_atoms = MSC::kd_create(3);

        // add all atom positions (in grid coordinates) to a kdtree
        for(int i = 0; i < num_atoms; i++) {

            double gpos[3];
            m_metadata.world_to_grid (m_atoms[i].m_pos, gpos);

            kdtree_add(m_kdtree_atoms, gpos[0], gpos[1], gpos[2], i+1);    // VASP atom ids start with 1
        }

        // now, use kdtree to find nearest atom to every extremum
        for(int extIdx = 0; extIdx < extrema.size(); extIdx++){

            INDEX_TYPE vertIdx = extrema[extIdx];

            double fcoords[3];
            m_metadata.idx_to_grid(vertIdx, fcoords);

            int nearestAtomIdx = kdtree_query(m_kdtree_atoms, fcoords);

            if (nearestAtomIdx == -1) {
                printf("TopoMS::bader(). Could not find nearest atom to (%f %f %f)\n", fcoords[0], fcoords[1], fcoords[2]);
            }
            else {
                extrema2atoms.insert(std::make_pair( extIdx, nearestAtomIdx ));
            }
        }

        printf(" Done!");
        ltimer.EndGlobal ();
        ltimer.PrintAll ();
    }

    // ---------------------------------------------------------------------
    // integrate the volume labels
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Computing Bader Volumes...");
        fflush(stdout);

        labels_atoms.resize(numlabels, -1);
        chg_atoms.resize(num_atoms+1, 0.0);       // zeroth element will be vacuum
        vol_atoms.resize(num_atoms+1, 0.0);       // zeroth element will be vacuum

        labels_extrema.resize(numlabels, -1);
        chg_extrema.resize(numextrema+1, 0.0);
        vol_extrema.resize(numextrema+1, 0.0);

#ifdef USE_KAHAN_SUM
        // temporary for kahan
        Utils::Kahan::KahanObject<FLOATTYPE> kinit {0};

        std::vector<Utils::Kahan::KahanObject<FLOATTYPE> > kchg_atoms(num_atoms+1, kinit);
        std::vector<Utils::Kahan::KahanObject<FLOATTYPE> > kchg_extrema(numextrema+1, kinit);
#endif

        for(INDEX_TYPE vIdx = 0; vIdx < numlabels; vIdx++){

            const FLOATTYPE value = (this->m_negated ? -1 : 1) * m_gridfunc->SampleImage(vIdx);

            int extIdx = (*volumelabeling)[vIdx];

            int atomIdx = ( fabs(value) <= vacthreshold_in_fileUnits || extIdx == -2 ) ? 0 : extrema2atoms[ extIdx ];
            extIdx = (extIdx == -2) ? 0 : extIdx+1;

            labels_atoms[vIdx] = atomIdx;
            vol_atoms[atomIdx] += 1.0;

            labels_extrema[vIdx] = extIdx;
            vol_extrema[extIdx] += 1.0;

#ifdef USE_KAHAN_SUM
            kchg_atoms[atomIdx] = Utils::Kahan::KahanSum( kchg_atoms[atomIdx], value );
            kchg_extrema[extIdx] = Utils::Kahan::KahanSum( kchg_extrema[extIdx], value );
#else
            chg_atoms[atomIdx] += value;
            chg_extrema[maxIdx] += value;
#endif
        }

#ifdef USE_KAHAN_SUM
        // convert from Kahan to normal
        for(int k = 0; k < chg_atoms.size(); k++){      chg_atoms[k] = kchg_atoms[k].sum;       }
        for(int k = 0; k < chg_extrema.size(); k++){    chg_extrema[k] = kchg_extrema[k].sum;   }
#endif
        // ---------------------------------------------------------------------
        // normalize

        {
#ifdef USE_KAHAN_SUM
            FLOATTYPE vsm = Utils::Kahan::sum(vol_atoms);
            FLOATTYPE sum_atoms = Utils::Kahan::sum(chg_atoms);
            FLOATTYPE sum_extrema = Utils::Kahan::sum(chg_extrema);
#else
            FLOATTYPE vsm = Utils::sum(vol_atoms);
            FLOATTYPE sum_atoms = Utils::sum(chg_atoms);
            FLOATTYPE sum_extrema = Utils::sum(chg_extrema);
#endif
            if( size_t(vsm) != numlabels) {
                printf(" \n ---> error in vol integration: %d != %f\n", numlabels, vsm);
                //exit(1);
            }

            if ( fabs(sum_atoms - sum_func) > powf(10, -6) || fabs(sum_extrema - sum_func) > powf(10, -6)) {
                printf(" \n ---> error in chg integration: %f != %f, %f (%f, %f)\n", sum_func,
                       sum_atoms, sum_extrema, (sum_atoms-sum_func), (sum_extrema-sum_func));
                //exit(1);
            }
        }

        std::transform(vol_atoms.begin(), vol_atoms.end(), vol_atoms.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), vol_voxel));
        std::transform(chg_atoms.begin(), chg_atoms.end(), chg_atoms.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), chgDens_fileUnit2e));

        std::transform(vol_extrema.begin(), vol_extrema.end(), vol_extrema.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), vol_voxel));
        std::transform(chg_extrema.begin(), chg_extrema.end(), chg_extrema.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), chgDens_fileUnit2e));

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    // ---------------------------------------------------------------------
    // print summary
    // ---------------------------------------------------------------------

#ifdef USE_KAHAN_SUM
    FLOATTYPE vsum = Utils::Kahan::sum(vol_atoms);
    FLOATTYPE csum = Utils::Kahan::sum(chg_atoms);
#else
    FLOATTYPE vsum = Utils::sum(vol_atoms);
    FLOATTYPE csum = Utils::sum(chg_atoms);
#endif
    printf("\tTotal volume = %f, vaccum volume = %f.\n", vsum, vol_atoms.front());
    printf("\tTotal charge = %f, vaccum charge = %f.\n", csum, chg_atoms.front());
    //printf(" sum of function over all grid points = %f, %f\n", get_fsum(), get_fint());

    printf("\n =====>> Bader analysis finished!");
    timer.EndGlobal ();
    timer.PrintAll ();
    Utils::print_separator();

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------
    return true;
}

/**
  *   @brief  Extract the molecular graph using MSC library
  *
  *   @return success flag
  */

bool TopoMS::msc() {

    if(!m_config->msc)
        return false;

    Utils::print_separator();
    printf(" =====>> Extracting Molecular Graph\n\n");

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    m_tgrid = new MSC::GridType(m_grid);
    //RunMeshConsistencyChecks(m_tgrid);

    // ---------------------------------------------------------------------
    // compute a dense function (defined on every cell of the topo mesh)
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Creating TopologicalExplicitDenseMeshFunction...");
        fflush(stdout);

        m_topofunc = new MSC::TopologicalExplicitDenseMeshFunction();
        m_topofunc->setMeshAndAllocate(m_tgrid);
        m_topofunc->copyVertexValuesFromGridFunction(m_gridfunc);
        m_topofunc->setCellValuesMaxOfVerts();

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    // ---------------------------------------------------------------------
    // restrict the discrete gradient based on numerical integration
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Computing volume restriction for a consistent MSC...");
        fflush(stdout);

        MSC::VertexLabelingToBoundaryLabeling<INDEX_TYPE>* edgemap = new MSC::VertexLabelingToBoundaryLabeling<INDEX_TYPE>(m_integrator2->GetOutputLabels(), m_tgrid);
        edgemap->ComputeBoundary();
#ifdef OUTPUT_DEBUG_FIELDS
        edgemap->GetOutputLabels()->OutputToFile("boundary_labels.raw");
        edgemap->OutputEdgesToFile("surfin.raw");
#endif

        m_tgrid->set_restriction (edgemap->GetOutputLabels ());

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    const FLOATTYPE vacthreshold = m_config->vacuum_threshold * m_metadata.volume();

    // ---------------------------------------------------------------------
    // restrict the discrete gradient only in nonvacuum regions
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Creating masked mesh for a sparse MSC...");
        fflush(stdout);

        MSC::DenseLabeling<char>* cellmask = new MSC::DenseLabeling<char>(m_tgrid->numCells());
        cellmask->SetAll(1);

        #pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            int thread_num = omp_get_thread_num();

            std::vector<INDEX_TYPE> partition;
            MSC::ArrayIndexPartitioner::EvenChunkSplit(m_tgrid->numCells(), num_threads, partition);

            MSC::GridType::TopologicalRegularGrid::DCellsIterator vit(m_tgrid, 0, partition[thread_num], partition[thread_num + 1]);
            INDEX_TYPE vid;
            for (vit.begin(); vit.valid(); vit.advance()) {

                vid = vit.value();
                if ( fabs( m_topofunc->cellValue(vid) ) > vacthreshold ) {
                    continue;
                }

                cellmask->SetLabel(vid, 0);

                MSC::GridType::TopologicalRegularGrid::AdjacentCellsIterator cviter(m_tgrid);
                for (cviter.begin(vid); cviter.valid(); cviter.advance()) {
                    cellmask->SetLabel(cviter.value(), 0);
                }
            }
        }

        //GInt::CellTesterDefaultTrue* cellmasktester = new GInt::CellTesterDefaultTrue();
        MSC::CellTesterLaberInput* cellmasktester = new MSC::CellTesterLaberInput();
        cellmasktester->SetLabeling(cellmask);
        m_tgrid->SetTester(cellmasktester);

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();

#ifdef OUTPUT_DEBUG_FIELDS
        cellmask->OutputToIntFile("cellmask.raw");
#endif
    }

    // ---------------------------------------------------------------------
    // based on the above, create the discrete gradient
    // ---------------------------------------------------------------------
    MSC::DiscreteGradientLabeling *labeling = new MSC::DiscreteGradientLabeling(m_tgrid);
    labeling->ClearAllGradient();

    MSC::RobinsLabelingAlgorithm<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction> *robin =
            new MSC::RobinsLabelingAlgorithm<MSC::GridType,MSC::TopologicalExplicitDenseMeshFunction>(m_topofunc, m_tgrid, labeling);
    robin->compute_output();
    robin->summarize();

    // ---------------------------------------------------------------------
    // compute full MSC
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Setting dimensionality of ascending manifolds...");
        fflush(stdout);

        MSC::TopologicalGradientUsingAlgorithms<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction, MSC::DiscreteGradientLabeling>
                * topo_algs = new MSC::TopologicalGradientUsingAlgorithms<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction, MSC::DiscreteGradientLabeling>(m_topofunc, m_tgrid, labeling);
        topo_algs->setAscendingManifoldDimensions();

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
#if 0
        ltimer.StartGlobal();

        printf(" -- Checking gradient consistency...");
        fflush(stdout);

        topo_algs->CheckGradientConsistency();

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
#endif
    }

#ifdef OUTPUT_DEBUG_FIELDS
    {
        char gradname2[1024];
        sprintf(gradname2, "%s.grad", m_config->infilename.c_str ());

        printf(" -- Writing gradient to file %s...", gradname2);
        fflush(stdout);

        labeling->outputToFile(gradname2);

        printf("Done!");
    }
#endif
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();
        printf("\n -- Creating MSC...\n");

        m_msc = new MSC::MSCType(labeling, m_tgrid, m_topofunc);
        m_msc->ComputeFromGrad();

        printf(" -- MSC creation done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    // ---------------------------------------------------------------------
    // restrict cancellation to non-atom minima only
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Ensuring consistent simplification...");
        fflush(stdout);

        // add in selector shit here
        m_kdtree_minima = MSC::kd_create(3);

        size_t cnt = 0;

        // add all my minima to the kdtree
        for (int i = 0; i < m_msc->numNodes(); i++) {

            MSC::node<FLOATTYPE>& n = m_msc->getNode(i);
            if (n.dim != 0) continue;

            if ( fabs( n.value ) < vacthreshold ) {
                continue;
            }

            MSC::Vec3l ncoords;
            m_tgrid->cellid2Coords(n.cellindex, ncoords);

            MSC::Vec3d ndcoords = ncoords;
            ndcoords *= 0.5;

            kdtree_add(m_kdtree_minima, ndcoords[0], ndcoords[1], ndcoords[2], i);
            cnt++;
        }

        // now i added all the "atoms"
        // for each node find closest "atom"
        // because i dont have atoms file, simply use a perturbed version of each node

        // i save the atom positions as well as the index so that i can render later the atom position.
        vector<pair<MSC::Vec3d, INT_TYPE> > restrictnodes;

        for(size_t i = 0; i < m_atoms.size (); i++) {

            double gpos[3];
            m_metadata.world_to_grid (m_atoms[i].m_pos, gpos);

            int nearestIdx = kdtree_query(m_kdtree_minima, gpos);

            if (nearestIdx == -1) {
                printf("TopoMS::msc(). Could not find nearest node to atom (%f %f %f)\n", gpos[0], gpos[1], gpos[2]);
            }
            else {
                restrictnodes.push_back(std::pair<MSC::Vec3d, INT_TYPE>(MSC::Vec3d(gpos[0], gpos[1], gpos[2]), nearestIdx));
            }
        }
        printf(" (by forbidding cancellation of %d nodes)...", restrictnodes.size());
        fflush(stdout);

        // now this code below actually forces the MSC to avoid cancelling the right atoms
        for (int i = 0; i < restrictnodes.size(); i++) {
            INT_TYPE nodeid = restrictnodes[i].second;
            m_msc->getNode(nodeid).boundary = 128; // make boundary type
        }

        printf("  Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    // ---------------------------------------------------------------------
    // simplify MSC now
    // ---------------------------------------------------------------------
    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf("\n -- Simplifying MSC...\n");
        //m_msc->ComputeHierarchy( this->pers_perc2val(2));
        //m_msc->ComputeHierarchy( this->pers_perc2val(10.0) );
        m_msc->ComputeHierarchy(100.0);
        printf(" -- MSC Simplification done!");
        fflush(stdout);
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    // ---------------------------------------------------------------------
    // print total time
    // ---------------------------------------------------------------------

    printf("\n =====>> Molecular Graph Extracted!");
    timer.EndGlobal ();
    timer.PrintAll ();
    Utils::print_separator();

    // ---------------------------------------------------------------------
    simplify_msc( persistence_val, filter_val );
    return true;
}

/**
  *   @brief  Simplify a precomputed molecular graph
  *
  *   @param  pvalue is the persistence value
  *   @param  fvalue is the persistence value
  *   @return void
  */
void TopoMS::simplify_msc(FLOATTYPE pvalue, FLOATTYPE fvalue) {

    persistence_val = pvalue;
    filter_val = fvalue;

    // ----------------------------------------------------------
    printf("\n -- Extracting simplified topology (persistence = %f, filtered at %f)\n", persistence_val, filter_val);
    m_msc->SetSelectPersAbs(persistence_val);

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    this->m_paths.clear();
    this->m_nodes.clear();

    // now do the selection to show 0-1 arcs
    MSC::MSCSelectorLivingNodes<MSC::MSCType> *s_nodes = new MSC::MSCSelectorLivingNodes<MSC::MSCType>(m_msc);

    MSC::MSCSelectorNodeIndex<MSC::MSCType> *s_1saddles = new MSC::MSCSelectorNodeIndex<MSC::MSCType>(m_msc, 1);
    s_1saddles->add_parent(s_nodes);

    MSC::MSCSelectorRepresentative1Saddle<MSC::MSCType> *s_r1saddles = new MSC::MSCSelectorRepresentative1Saddle<MSC::MSCType>(m_msc);
    s_r1saddles->add_parent(s_1saddles);

    // here is how to actually extract arcs
    MSC::StreamlineIntegratorType *sintegrator = new MSC::StreamlineIntegratorType(m_grid, m_gridfunc, m_config->error_threshold, m_config->grad_threshold, m_config->numiter);
    MSC::TerminateNearExtrema *extremumtermination = new MSC::TerminateNearExtrema(m_integrator->GetExtrema(), m_grid);
    sintegrator->SetAdvectionChecker(extremumtermination);

    s_r1saddles->compute_output();

    // ----------------------------------------------------------
    // compute paths
    static const float eps = 0.000001;
    for (auto it = s_r1saddles->output.begin(); it != s_r1saddles->output.end(); it++) {

        // iterate over arcs around the saddle
        MSC::MSCType::SurroundingArcsIterator sit(m_msc);
        for (sit.begin(*it); sit.valid(); sit.advance()) {

            INT_TYPE arcid = sit.value();
            MSC::arc<FLOATTYPE>& arc = m_msc->getArc(arcid);

            // only look at saddle-minimum arcs
            if (arc.upper != *it) {
                continue;
            }

            INDEX_TYPE lowercellid = m_msc->getNode(arc.lower).cellindex;       // combinatorial says this is where we should end up
            INDEX_TYPE uppercellid = m_msc->getNode(arc.upper).cellindex;       // combinatorial says this is where we should end up

            INDEX_TYPE lowervertid = m_tgrid->VertexNumberFromCellID(lowercellid);  // get vertex id in grid coordinates from combinatorial min from topological coordinates


            // filtering
            if (fabs(m_topofunc->cellValue(lowercellid)) <= filter_val || fabs(m_topofunc->cellValue(uppercellid)) <= filter_val){
                continue;
            }

            vector<INDEX_TYPE> geom_combinatorial;
            m_msc->fillArcGeometry(arcid, geom_combinatorial);              // combinatorial arc geometry

            std::vector<MSC::Vec3d> geom_numerical;                        // structure to hold numerical path

            // now try integrating down starting from cominatorial geomery
            for (size_t j = 0; j < geom_combinatorial.size(); j++) {

                MSC::Vec3l c;
                m_tgrid->cellid2Coords(geom_combinatorial[j], c);
                MSC::Vec3d cf = c;
                cf *= 0.5;


                if (j > 0) {

                    // ignore overlapping points as they cause problem while rendering
                    // tube computation requires nondegenerate sequence of points
                    if ((geom_numerical.back() - cf).Mag() < eps)
                        continue;
                }

                geom_numerical.push_back(cf); // record combinatorial step as part of geometric path

                std::vector<MSC::Vec3d> sline_num;
                MSC::ADVECTION_EVENT event = sintegrator->IntegrateStreamline(cf, sline_num); // do numerical streamline integration
                if (event == MSC::ADVECTION_EVENT::LOW_GRADIENT) {               //glColor3f(0.5, .5, .9);
                } else if (event == MSC::ADVECTION_EVENT::HIT_PREASSIGNED) {     //glColor3f(0.5, .9, .2);
                } else if (event == MSC::ADVECTION_EVENT::OVER_MAX_ITERATIONS) { //glColor3f(0.9, 0.5, 0.2);
                } else {                                                         //glColor3f(.2, .2, .2);
                }

                // test if this is located near a critical vertex that was used as a stopping criteria in region integration
                if (extremumtermination->WhichExtremum(sline_num.back()) == lowervertid) {

                    // this numerical path is good! so paste it onto rest of geometry and we're done with this arc
                    geom_numerical.insert(geom_numerical.end(), sline_num.begin(), sline_num.end());
                    break;
                }
            }

            this->m_paths.push_back(geom_numerical);
            this->m_nodes.push_back(arc.lower);
            this->m_nodes.push_back(arc.upper);
        }
    }

    // count number of critical points
    size_t counts[4] = {0,0,0,0};
    for (auto it = s_nodes->output.begin(); it != s_nodes->output.end(); it++) {

        MSC::node<FLOATTYPE>& node = m_msc->getNode(*it);
        counts[node.dim]++;
    }

    printf(" -- Simplified topology extracted! Found %d bond paths!", this->m_paths.size());
    printf("    -- # crit pts = [%d, %d, %d, %d] = %d\n", counts[0], counts[1], counts[2], counts[3],
                                                          counts[0] + counts[1] + counts[2] + counts[3] );
    timer.EndGlobal ();
    timer.PrintAll ();
    Utils::print_separator();
}


// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
/**
  *   @brief  Write bader charges per atom
  */

void TopoMS::write_bader_atoms2chgvol(const std::string &filename) const {
    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    FILE *outfile = fopen(filename.c_str(), "w");
    if(!outfile){
        printf(" could not open file\n");
        return;
    }
    fprintf(outfile, " %4s %11s %11s %11s %11s %13s\n",
                        "#", "X", "Y", "Z", "CHARGE", "VOLUME");
    fprintf(outfile, "------------------------------------------------------------------------------\n");

    for(int i = 0; i < m_atoms.size (); i++){
        fprintf(outfile, " %4d %11.7f %11.7f %11.7f %11.7f %13.7f\n",
                i+1, m_atoms[i].m_pos[0], m_atoms[i].m_pos[1], m_atoms[i].m_pos[2],
                     chg_atoms[i+1], vol_atoms[i+1]);
    }
    fprintf(outfile, " %6s %33s %11.7f %13.7f\n", "Vacuum", "", chg_atoms[0], vol_atoms[0]);
    fprintf(outfile, "------------------------------------------------------------------------------\n");
    fprintf(outfile, " %6s %33s %11.7f %13.7f\n", "Total", "",
            std::accumulate(chg_atoms.begin(), chg_atoms.end(), 0.0),
            std::accumulate(vol_atoms.begin(), vol_atoms.end(), 0.0)
           );
    fprintf(outfile, "------------------------------------------------------------------------------\n");
    fclose(outfile);
    printf(" Done!\n");
}


/**
  *   @brief  Write bader charges per maximum
  */

void TopoMS::write_bader_max2chgvol(const std::string &filename, bool filter_small) const {
#if 0
    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    FILE *outfile = fopen(filename.c_str(), "w");
    if(!outfile){
        printf(" could not open file\n");
        return;
    }

    fprintf(outfile, " %4s %11s %11s %11s %13s %8s\n",
                        "#", "X", "Y", "Z", "CHARGE", "ATOM");
    fprintf(outfile, "---------------------------------------------------------------------\n");

    int c = 0;
    for(int i = 0; i < minima.size (); i++){

        if (filter_small && chg_extrema[i] < powf(10, -6))
            continue;

        int maxCellIdx = minima.at(i);

        double xyzIdx[3];
        m_metadata.idx_to_grid(maxCellIdx, xyzIdx);

        Vec3f maxPos;
        for(int d = 0; d < 3; d++)
            maxPos[d] = m_metadata.grid_to_world ( xyzIdx[d], d );


        fprintf(outfile, " %4d %11.7f %11.7f %11.7f %13.7f %8d\n",
                ++c, maxPos[0], maxPos[1], maxPos[2],
                chg_extrema[i], extrema2atoms.at (i));
    }
    fprintf(outfile, "---------------------------------------------------------------------\n");
    fclose(outfile);
    printf(" Done!\n");
#endif
}


// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// vtk headers
#ifdef USE_VTK
    #include <vtkSmartPointer.h>
    #include <vtkIntArray.h>
    #include <vtkFloatArray.h>
    #include <vtkPoints.h>
    #include <vtkVertex.h>
    #include <vtkPolyLine.h>
    #include <vtkCellArray.h>
    #include <vtkPointData.h>
    #include <vtkCellData.h>
    #include <vtkFieldData.h>
    #include <vtkPolyData.h>
    #include <vtkImageData.h>
    #include <vtkXMLPolyDataWriter.h>
    #include <vtkXMLImageDataWriter.h>
#endif

/**
  *   @brief  Write bader maximum labels as volumes (VTK file)
  */
void TopoMS::write_bader_max2vol(const string &filename) const {

#ifndef USE_VTK
    printf("TopoMS::write_bader_max2vol. VTK not available. Writing as raw file instead.\n");

    FILE* fout = fopen(filename.c_str(), "wb");
    fwrite(labels_extrema.data(), sizeof(int), labels_extrema.size(), fout);
    fclose(fout);
#else

    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

    imageData->SetOrigin(0.0, 0.0, 0.0);
    imageData->SetSpacing( m_metadata.m_lattice_vectors[0][0] / float(m_metadata.m_grid_dims[0]-1),
                           m_metadata.m_lattice_vectors[1][1] / float(m_metadata.m_grid_dims[1]-1),
                           m_metadata.m_lattice_vectors[2][2] / float(m_metadata.m_grid_dims[2]-1)
                         );
    imageData->SetDimensions( m_metadata.m_grid_dims[0], m_metadata.m_grid_dims[1], m_metadata.m_grid_dims[2] );
    imageData->AllocateScalars(VTK_INT, 1);
    imageData->GetPointData()->GetScalars()->SetName("max_labeling");

    // Fill every entry of the image data with "2.0"
    for (size_t z = 0; z < m_metadata.m_grid_dims[2]; z++) {
    for (size_t y = 0; y < m_metadata.m_grid_dims[1]; y++) {
    for (size_t x = 0; x < m_metadata.m_grid_dims[0]; x++) {

        int* pixel = static_cast<int*>(imageData->GetScalarPointer(x,y,z));
        pixel[0] = labels_extrema[ m_metadata.grid_to_idx(x,y,z) ];
    }
    }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(imageData);
    writer->Write();
    printf(" Done!\n");
#endif
}

/**
  *   @brief  Write bader atom labels as volumes (VTK file)
  */
void TopoMS::write_bader_atoms2vol(const string &filename) const {

#ifndef USE_VTK
    printf("TopoMS::write_bader_atoms2vol. VTK not available. Writing as raw file instead.\n");

    FILE* fout = fopen(filename.c_str(), "wb");
    fwrite(labels_atoms.data(), sizeof(int), labels_atoms.size(), fout);
    fclose(fout);
#else
    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

    imageData->SetOrigin(0.0, 0.0, 0.0);
    imageData->SetSpacing( m_metadata.m_lattice_vectors[0][0] / float(m_metadata.m_grid_dims[0]-1),
                           m_metadata.m_lattice_vectors[1][1] / float(m_metadata.m_grid_dims[1]-1),
                           m_metadata.m_lattice_vectors[2][2] / float(m_metadata.m_grid_dims[2]-1)
                         );
    imageData->SetDimensions( m_metadata.m_grid_dims[0], m_metadata.m_grid_dims[1], m_metadata.m_grid_dims[2] );
    imageData->AllocateScalars(VTK_INT, 1);
    imageData->GetPointData()->GetScalars()->SetName("atom_labeling");

    // Fill every entry of the image data with "2.0"
    for (size_t z = 0; z < m_metadata.m_grid_dims[2]; z++) {
    for (size_t y = 0; y < m_metadata.m_grid_dims[1]; y++) {
    for (size_t x = 0; x < m_metadata.m_grid_dims[0]; x++) {

        int* pixel = static_cast<int*>(imageData->GetScalarPointer(x,y,z));
        pixel[0] = labels_atoms[ m_metadata.grid_to_idx(x,y,z) ];
    }
    }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(imageData);
    writer->Write();
    printf(" Done!\n");
#endif
}

/**
  *   @brief  Write molecular graph (VTK file)
  */
void TopoMS::write_mgraph(const std::string &filename) const {

#ifndef USE_VTK
    printf("TopoMS::write_mgraph. VTK not available\n");
#else

    printf(" - Writing Molecular Graph to %s...", filename.c_str());
    fflush(stdout);

    static const float periodic_cutoff = 0.3*this->get_gridDims()[0];
    static const float eps = 0.000001;

    // ------------------------------------------------------------------------------
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> critPts = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkIntArray> pDims = vtkSmartPointer<vtkIntArray>::New();
    pDims->SetName("cp_dim");
    pDims->SetNumberOfComponents(1);

    /*vtkSmartPointer<vtkIntArray> cDims = vtkSmartPointer<vtkIntArray>::New();
    cDims->SetName("cell_dims");
    cDims->SetNumberOfComponents(1);*/

    // ------------------------------------------------------------------------------
    // Create a cell array to store the critical points
    for(size_t i = 0; i < this->m_nodes.size(); i++) {

        int ndim;
        INDEX_TYPE ncidx;
        MSC::Vec3l coord;

        this->get_msc_node(m_nodes[i], ndim, ncidx, coord);

        // get_msc_node returns MSC grid coordinates (twice the actual value)
        points->InsertNextPoint( 0.5*this->m_metadata.grid_to_world(coord[0], 0),
                                 0.5*this->m_metadata.grid_to_world(coord[1], 1),
                                 0.5*this->m_metadata.grid_to_world(coord[2], 2)
                                );

        vtkSmartPointer<vtkVertex> cPnt = vtkSmartPointer<vtkVertex>::New();
        cPnt->GetPointIds()->InsertNextId(points->GetNumberOfPoints()-1);
        critPts->InsertNextCell(cPnt);
        pDims->InsertNextValue( i%2 == 0 ? 3 : 2 );
        //cDims->InsertNextValue( i%2 == 0 ? 3 : 2 );
    }

    // ------------------------------------------------------------------------------
    // Create a cell array to store the lines in and add the lines to it
    for(size_t i = 0; i < this->m_paths.size(); i++) {

        const std::vector<MSC::Vec3d> &path = this->m_paths[i];

        vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
        for(unsigned int pidx = 0; pidx < path.size(); pidx++) {

            if (pidx > 0) {

                double disp = (path[pidx] - path[pidx-1]).Mag();

                if (disp < eps){
                    continue;
                }

                if (disp > periodic_cutoff) {
                    lines->InsertNextCell(polyLine);
                    //cDims->InsertNextValue(6);
                    polyLine = vtkSmartPointer<vtkPolyLine>::New();
                }
            }

            points->InsertNextPoint( this->m_metadata.grid_to_world(path[pidx][0], 0),
                                     this->m_metadata.grid_to_world(path[pidx][1], 1),
                                     this->m_metadata.grid_to_world(path[pidx][2], 2)
                                    );
            pDims->InsertNextValue(6);
            polyLine->GetPointIds()->InsertNextId(points->GetNumberOfPoints()-1);
        }

        lines->InsertNextCell(polyLine);
        //cDims->InsertNextValue(6);
    }

    // -----------------------------------------------------------------------------
    // add persistence and filtering value to field data
    vtkSmartPointer<vtkFloatArray> pValue = vtkSmartPointer<vtkFloatArray>::New();
    pValue->SetNumberOfComponents(1);
    pValue->SetName("persistence");
    pValue->InsertNextValue( this->get_persistence() );

    vtkSmartPointer<vtkFloatArray> fValue = vtkSmartPointer<vtkFloatArray>::New();
    fValue->SetNumberOfComponents(1);
    fValue->SetName("filter");
    fValue->InsertNextValue( this->get_filterval() );

    // ------------------------------------------------------------------------------
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

    polyData->SetPoints(points);
    polyData->SetLines(lines);
    polyData->SetVerts(critPts);
    polyData->GetPointData()->AddArray(pDims);
    //polyData->GetCellData()->AddArray(cDims);
    polyData->GetFieldData()->AddArray(pValue);
    polyData->GetFieldData()->AddArray(fValue);

    // ------------------------------------------------------------------------------
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetInputData(polyData);
    writer->SetFileName(filename.c_str());
    writer->Write();

    printf(" Done!\n");
#endif
}
