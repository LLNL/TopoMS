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
 *
 *  @brief This file provides the core functionality for TopoMS
 *
 *  @section DESCRIPTION
 *
 *  This file provides the core functionality for TopoMS
 *
 */

#include <cassert>
#include <functional>

// -----------------------------------------------------------------------
// TopoMS headers
#include "ConfigParser.h"
#include "InputFormats.h"
#include "MolecularSystem.h"
#include "TopoMS.h"
#include "MSCBond.h"

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
#include "msc_iterators.h"
#include "numeric_streamline_integrator.h"

#ifdef USE_VTK
    #include "vtkVolumeSlicer.h"
#endif

// -----------------------------------------------------------------------------------------
// these wrappers are needed so rest of the code remains independent of the MSC headers
// -----------------------------------------------------------------------------------------

bool TopoMS::msc_is_available() const {              return this->m_msc != nullptr;         }
size_t TopoMS::msc_get_nnodes() const {              return this->m_msc->numNodes();        }
size_t TopoMS::msc_get_narcs() const {               return this->m_msc->numArcs();         }

bool TopoMS::msc_is_nodeAlive(size_t idx) const {   return this->m_msc->isNodeAlive(idx);   }
bool TopoMS::msc_is_arcAlive(size_t idx) const {    return this->m_msc->isArcAlive(idx);    }
bool TopoMS::msc_is_nodeFiltered(size_t idx) const {
    return (fabs(this->m_msc->getNode(idx).value) > fabs(this->filter_val));
}

int TopoMS::msc_get_nodeDim(size_t idx) const {
    const MSC::node<FLOATTYPE> &n = this->m_msc->getNode(idx);
    return (this->m_negated) ? 3-n.dim : n.dim;
}
int TopoMS::msc_get_arcDim(size_t idx) const {
    if(this->m_negated)     return this->msc_get_nodeDim(this->m_msc->getArc(idx).upper);
    else                    return this->msc_get_nodeDim(this->m_msc->getArc(idx).lower);
}

void TopoMS::msc_cellid_to_gcoords(const INDEX_TYPE &cellIdx, MSC::Vec3d &ncoords) const {

    MSC::Vec3l lcoord;
    m_tgrid->cellid2Coords(cellIdx, lcoord);
    ncoords = MSC::Vec3d(lcoord);
    ncoords *= 0.5;
}
void TopoMS::msc_gcoords_to_wcoords(const MSC::Vec3d &ncoords, float wcoords[]) const {
    float gcoords[3] = {ncoords[0], ncoords[1], ncoords[2]};
    m_metadata.grid_to_world(gcoords, wcoords);
}
void TopoMS::msc_cellid_to_wcoords(const INDEX_TYPE &cellIdx, float wcoords[]) const {

    MSC::Vec3d ncoords;
    this->msc_cellid_to_gcoords(cellIdx, ncoords);

    float gcoords[3] = {ncoords[0], ncoords[1], ncoords[2]};
    m_metadata.grid_to_world(gcoords, wcoords);
}

bool TopoMS::msc_get_node(size_t idx, int &dim, INDEX_TYPE &cellIdx, MSC::Vec3d &ncoords) const {

    const MSC::node<FLOATTYPE> &n = this->m_msc->getNode(idx);
    cellIdx = n.cellindex;

    dim = (this->m_negated) ? 3-n.dim : n.dim;
    this->msc_cellid_to_gcoords(cellIdx, ncoords);
    return true;
}

void TopoMS::msc_get_arcGeometry(size_t idx, std::vector<INDEX_TYPE> &arc_geom) const {
    this->m_msc->fillArcGeometry(idx, arc_geom);
}

void TopoMS::msc_print_node(size_t idx) const {

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

    //MSC::Vec3d ncoords;
    //this->msc_cellid_to_gcoords(n.cellindex, ncoords);

    float ncoords[3];
    this->msc_cellid_to_wcoords(n.cellindex, ncoords);

    printf(" %s (node %d, cell %d) at (%.3f %.3f %.3f) val = %f\n",
                  type.c_str(), idx, n.cellindex,
                  ncoords[0], ncoords[1], ncoords[2], val);
}

// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
const std::vector<long long> &TopoMS::bader_get_extrema() const {
    if(m_integrator == 0)
        return std::vector<INDEX_TYPE>();
    return m_integrator->GetExtrema();
}

int TopoMS::bader_get_atomLabeling(INDEX_TYPE vIdx) const {
    int extIdx = (*m_volumelabeling)[vIdx];
	int atomIdx;
	if (extIdx == -2) {
		atomIdx = 0;
	}
	else {
		atomIdx = extrema2atoms.at(extIdx);
	}
    return atomIdx;
}
int TopoMS::bader_get_atomLabeling(size_t x, size_t y, size_t z) const {
    return this->bader_get_atomLabeling(this->m_metadata.grid_to_idx(x,y,z));
}
int TopoMS::bader_get_atomLabeling(double pos[3]) const {
    return this->bader_get_atomLabeling(size_t(pos[0]), size_t(pos[1]), size_t(pos[2]));
}

// -----------------------------------------------------------------------------------------
// functions for kd-tree handling
// -----------------------------------------------------------------------------------------

bool TopoMS::kdtree_add(MSC::kdtree* kt, double x, double y, double z, int data) const {

    static const INDEX_TYPE X = m_metadata.m_grid_dims[0];
    static const INDEX_TYPE Y = m_metadata.m_grid_dims[1];
    static const INDEX_TYPE Z = m_metadata.m_grid_dims[2];

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

bool TopoMS::load(const std::string &configfilename, std::string datafile) {

    Utils::print_separator();
    std::cout << "\n =====>> Loading data using " << configfilename << "...\n";

    m_config = new Config(configfilename);
    m_config->parse();

    // overwrite the file
    if (datafile.length() > 0) {
        std::cout << " Overwriting the input file ("<<m_config->infilename<<") with argument ("<<datafile<<")\n";
        m_config->infilename = datafile;
    }

    //std::string extn = Utils::get_extn(m_config->infilename);
    m_datadir = Utils::get_directory(m_config->infilename);

    // -----------------------------------------------------------------------

    if (m_config->infiletype == "VASP") {
        m_mscfunc = MS::VASP::read_CAR<FLOATTYPE>(m_config->infilename, m_metadata, m_atoms);
        this->m_inputtype = IT_VASP;
    }
    else if (m_config->infiletype == "CUBE") {
        m_mscfunc = MS::Cube::read<FLOATTYPE>(m_config->infilename, m_metadata, m_atoms);
        this->m_inputtype = IT_CUBE;
    }

    if (m_metadata.grid_sz() == 0) {
        std::cerr << " Failed to read the input file correctly! got 0 grid values!\n";
        exit(1);
    }
    for(uint8_t i = 0; i < 3; i++){
        if (m_metadata.m_grid_dims[i] < 2) {
            std::cerr << "Incorrect grid! [" <<
                         m_metadata.m_grid_dims[0] << ", " <<
                         m_metadata.m_grid_dims[1] << ", " <<
                         m_metadata.m_grid_dims[2] << "]\n";
            exit(1);
        }
    }

    m_baderfunc = m_mscfunc;
    if (m_config->infiletype == "VASP" && m_config->reffilename.length() > 0) {

        MS::SystemInfo ref_metadata;
        std::vector<Material> ref_atoms;

        m_baderfunc = MS::VASP::read_CAR<FLOATTYPE>(m_config->reffilename, ref_metadata, ref_atoms);

        // make sure the two files give the same atoms
        if (ref_atoms.size() != m_atoms.size()) {
            std::cerr << "Mismatch in number of atoms given by "
                      << "(" << m_config->infilename << ") [" << m_atoms.size() << "] and "
                      << "(" << m_config->reffilename << ") [" << ref_atoms.size() << "] \n";
            exit(1);
        }

        for(size_t i = 0; i < m_atoms.size(); i++) {

            if (m_atoms[i] == ref_atoms[i])
                continue;

            std::cerr << "Mismatch in atom " << i << " : \n";
            m_atoms[i].print();
            ref_atoms[i].print();
            exit(1);
        }

        // make sure the two files give the same metadata
        if (!(m_metadata == ref_metadata)) {
            std::cerr << "Mismatch in metadata\n";
            m_metadata.print();
            ref_metadata.print();
            exit(1);
        }


        // if the user specifies a ref function,
        // the ref function should be used for graph computation
        // but function for bader charges

        // since the code is written completely to work with m_mscfunc,
        // we will swap the pointers
        // at the time of integrating bader charges, we will use m_baderfunc
        std::swap(m_mscfunc, m_baderfunc);
    }

    this->m_fieldtype = FT_UNKNOWN;
    if (m_config->fieldtype == "CHG") {             this->m_fieldtype = FT_CHG;         }
    else if (m_config->fieldtype == "POT") {        this->m_fieldtype = FT_POT;         }

    m_metadata.print ();

    // -----------------------------------------------------------------------
    std::cout << "    # atoms = " << m_atoms.size() << "\n";

    if (m_baderfunc == m_mscfunc) {
        std::pair<FLOATTYPE, FLOATTYPE> vr = get_frange(m_mscfunc, true);
        std::cout << "    function range = [" << vr.first << ", " << vr.second<< "]\n";
    }
    else {

        std::pair<FLOATTYPE, FLOATTYPE> vr = get_frange(m_baderfunc, true);
        std::cout << "    bader function range     = [" << vr.first << ", " << vr.second << "]\n";

        vr = get_frange(m_mscfunc, true);
        std::cout << "    msc (ref) function range = [" << vr.first << ", " << vr.second<< "]\n";
    }

    // -----------------------------------------------------------------------

#ifdef USE_VTK
    m_vtkFunction = this->create_vtkImagedata(m_baderfunc, "function");
    m_slicer_function = new vtkVolumeSlicer();
    m_slicer_function->set_volume(m_vtkFunction);
#endif

    return true;
}

// use configuration file to initialize the program
bool TopoMS::init() {

    Utils::print_separator();
    std::string ftype = "unknown";
    switch(this->m_fieldtype) {
    case FT_CHG:    ftype = "charge_density";   break;
    case FT_POT:    ftype = "local_potential";  break;
    }

    std::cout << "\n =====>> Initializing TopoMS! (field_type = " << ftype << ")\n";

    // -----------------------------------------------------------------------
    // compute bader for CHG, if msc needs to be computed
    if (m_fieldtype == FT_CHG && m_config->do_msc) {
        m_config->do_bader = true;
        std::cout << " -- Enabling Bader anlaysis since you asked for Molecular Graph!\n";
    }

    // bader must be computed only for CHG
    if (m_fieldtype != FT_CHG) {
        m_config->do_bader = false;
        std::cout << " -- Disabling Bader anlaysis since the input is not charge density!\n";
    }

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    m_grid = new MSC::RegularGrid(MSC::Vec3i(m_metadata.m_grid_dims[0], m_metadata.m_grid_dims[1], m_metadata.m_grid_dims[2]),
                                  MSC::Vec3b(m_config->is_periodic[0], m_config->is_periodic[1], m_config->is_periodic[2]));

    m_gridfunc = new MSC::RegularGridTrilinearFunction(m_grid, m_mscfunc);
    m_gridfunc->ComputeGradFromImage(1);

    if (1) {
        std::cout << " -- Negating the function.\n";
        m_negated = true;
        m_gridfunc->Negate();

        if (m_baderfunc == m_mscfunc) {

            std::pair<FLOATTYPE, FLOATTYPE> vr = get_frange(m_mscfunc, true);
            std::cout << "    -- negated function range = [" << vr.first << ", " << vr.second<< "]\n";
        }
        else {

            const size_t nsz = m_metadata.grid_sz();
#pragma omp parallel for schedule(static)
            for (INDEX_TYPE i = 0; i < nsz; i++) {
                m_baderfunc[i] *= -1;
            }

            // remember, we have swapped the two functions in the load()
            // since, this print is only for display, we use the *correct* names of the functions
            // to not confuse the user

            std::pair<FLOATTYPE, FLOATTYPE> vr = get_frange(m_baderfunc, true);
            std::cout << "    -- negated function range = [" << vr.first << ", " << vr.second<< "]\n";

            vr = get_frange(m_mscfunc, true);
            std::cout << "    -- negated ref_function range = [" << vr.first << ", " << vr.second<< "]\n";
        }
    }

    persistence_val = m_config->threshold_simp;
    filter_val = m_config->threshold_filt;

    std::cout << "\n =====>> TopoMS initialized!";
    timer.EndGlobal ();
    timer.PrintAll ();
    return true;
}


// -----------------------------------------------------------------------------------------
// bader analysis
// -----------------------------------------------------------------------------------------
/**
  *   @brief  Perform Bader analysis
  *
  *   @return success flag
  */

bool TopoMS::bader() {

    if(!m_config->do_bader)
        return false;

    // ---------------------------------------------------------------------
    // constant values required
    // ---------------------------------------------------------------------

    const size_t num_atoms = m_atoms.size();
    const size_t num_gridPts = m_metadata.grid_sz();

    const FLOATTYPE vol_box = m_metadata.volume_box();
    const FLOATTYPE vol_voxel = m_metadata.volume_voxel();
    const FLOATTYPE vol_box_angs = m_metadata.volume_box_in_Angs();
    const FLOATTYPE vol_voxel_angs = m_metadata.volume_voxel_in_Angs();

    const FLOATTYPE file_to_ae = m_metadata.volume_file2Angs() * m_metadata.charge_file2electrons();

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
    const FLOATTYPE chgDens_fileUnit2e = (this->m_inputtype == IT_CUBE ? vol_box * file_to_ae : 1.0) / (FLOATTYPE) num_gridPts;

    // ---------------------------------------------------------------------
    // for fast thresholding, convert the vacuum threshold to the same unit as the file (so i can do a point-wise comparison later)
    vacthreshold_in_fileUnits = m_config->threshold_vacuum * (this->m_inputtype == IT_CUBE ? 1.0 / file_to_ae : vol_box);

    // ---------------------------------------------------------------------
    // sum of function values
    const FLOATTYPE sum_func = (this->m_negated ? -1 : 1) * Utils::Kahan::sum(m_baderfunc, m_metadata.grid_sz());
    const FLOATTYPE total_chg = sum_func * chgDens_fileUnit2e;

    Utils::print_separator();
    printf("\n =====>> Performing Bader Analysis\n");
    if (m_metadata.is_lunit_Angstrom()) {
        printf("\t total volume = %f Ang^3\n", vol_box_angs);
        printf("\t voxel volume = %f Ang^3\n", vol_voxel_angs);
    }
    else {
        printf("\t total volume = %f Ang^3 (= %f %s^3)\n", vol_box_angs, vol_box, m_metadata.m_length_unit.c_str());
        printf("\t voxel volume = %f Ang^3 (= %f %s^3)\n", vol_voxel_angs, vol_voxel, m_metadata.m_length_unit.c_str());
    }
    printf("\t total chg = %f e\n", total_chg);
    printf("\t vacuum threshold = %E e per Ang^3 (= %E, in file units)\n\n", m_config->threshold_vacuum, vacthreshold_in_fileUnits);


    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    // ---------------------------------------------------------------------
    // volume decomposition
    // ---------------------------------------------------------------------

    m_volumelabeling = 0;
    {

        bool verbose = false;

        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Performing numeric integration for volume assignment");
        if (verbose) {  printf("\n");                   }
        else {          printf("...");  fflush(stdout); }

        // vacuum is applied on the m_gridfunc = m_mscfunc
        m_integrator = new MSC::IntegratorType(m_gridfunc, m_grid, m_config->threshold_error, m_config->threshold_grad, m_config->numiter);
        m_integrator->set_filter(vacthreshold_in_fileUnits);                         // the integrator will not seed a streamline for smaller values

        m_integrator->BeginIntegration(verbose);

#ifdef OUTPUT_DEBUG_FIELDS
        {
        std::string fname = m_datadir+"classes.raw";
        m_integrator->GetOutputLabels()->OutputToIntFile(fname.c_str());
        }
#endif

        // cleanup noise: do this only if you need msc as well
        if (this->m_config->do_msc) {
            m_integrator2 = new MSC::RegionRemoverType(m_gridfunc, m_integrator->GetOutputLabels(), &(m_integrator->GetExtrema()));
            m_integrator2->ComputeOutput(verbose);
            m_volumelabeling = m_integrator2->GetOutputLabelsUnmasked();
        }
        else {
            m_volumelabeling = m_integrator->GetOutputLabels();
        }

#ifdef OUTPUT_DEBUG_FIELDS
        std::string fname = m_datadir+"classes_clean.raw";
        m_integrator2->GetOutputLabels()->OutputToIntFile(fname.c_str());
#endif

        if(verbose) {   printf(" -- Done numeric integration of 3d volume!");   }
        else {          printf(" Done!");                                       }
        ltimer.EndGlobal ();
        ltimer.PrintAll ();
    }

    // ---------------------------------------------------------------------
    // ---------------------------------------------------------------------

    const std::vector<INDEX_TYPE> &extrema = m_integrator->GetExtrema();

    const size_t numextrema = extrema.size();
    const INDEX_TYPE numlabels = m_volumelabeling->GetNumLabels();

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
        for(size_t atomIdx = 0; atomIdx < num_atoms; atomIdx++) {

            float gpos[3];
            m_metadata.world_to_grid (m_atoms[atomIdx].m_pos, gpos);
            kdtree_add(m_kdtree_atoms, gpos[0], gpos[1], gpos[2], atomIdx+1);    // VASP atom ids start with 1
        }

        // now, use kdtree to find nearest atom to every extremum
        for(size_t extIdx = 0; extIdx < extrema.size(); extIdx++){

            INDEX_TYPE vertIdx = extrema[extIdx];

            double fcoords[3];
            m_metadata.idx_to_grid(vertIdx, fcoords);
            int nearestAtomIdx = kdtree_query(m_kdtree_atoms, fcoords);

            if (nearestAtomIdx == -1) {
                printf("TopoMS::bader(). Could not find nearest atom to (%f %f %f)\n", fcoords[0], fcoords[1], fcoords[2]);
            }
            else {
                extrema2atoms.insert(std::make_pair(extIdx, size_t(nearestAtomIdx)));
            }
        }

        printf(" Done!");
        ltimer.EndGlobal ();
        ltimer.PrintAll ();

        // check if all atoms have been accounted for!
        {
            std::set<size_t> test_atoms;
            for(auto iter = extrema2atoms.begin(); iter != extrema2atoms.end(); ++iter) {
                test_atoms.insert(iter->second);
            }
            if (test_atoms.size() != num_atoms) {
                std::cerr << " ERROR: mapping failed! " << num_atoms-test_atoms.size() << " fewer atoms mapped!\n";
                for(auto iter = extrema2atoms.begin(); iter != extrema2atoms.end(); ++iter) {
                    std::cout << "\t> extrema " << iter->first << " --> atom " << iter->second << std::endl;
                }
                exit(1);
            }
        }
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

        // temporary for kahan
        Utils::Kahan::KahanObject<FLOATTYPE> kinit {0};
        std::vector<Utils::Kahan::KahanObject<FLOATTYPE> > kchg_atoms(num_atoms+1, kinit);
        std::vector<Utils::Kahan::KahanObject<FLOATTYPE> > kchg_extrema(numextrema+1, kinit);

        for(INDEX_TYPE vIdx = 0; vIdx < numlabels; vIdx++){

            // 2018/10/10:  use "bader_func" which could be m_func or m_reffunc
            //const FLOATTYPE value = (this->m_negated ? -1 : 1) * m_gridfunc->SampleImage(vIdx);
            const FLOATTYPE value = (this->m_negated ? -1 : 1) * m_baderfunc[vIdx];

            int extIdx = (*m_volumelabeling)[vIdx];
            int atomIdx = (fabs(value) <= vacthreshold_in_fileUnits || extIdx == -2 ) ? 0 : extrema2atoms[extIdx];
            extIdx = (extIdx == -2) ? 0 : extIdx+1;

            labels_atoms[vIdx] = atomIdx;
            vol_atoms[atomIdx] += 1.0;

            labels_extrema[vIdx] = extIdx;
            vol_extrema[extIdx] += 1.0;

            kchg_atoms[atomIdx] = Utils::Kahan::KahanSum(kchg_atoms[atomIdx], value);
            kchg_extrema[extIdx] = Utils::Kahan::KahanSum(kchg_extrema[extIdx], value);
        }

        // convert from Kahan to normal
        for(size_t k = 0; k < chg_atoms.size(); k++){      chg_atoms[k] = kchg_atoms[k].sum;       }
        for(size_t k = 0; k < chg_extrema.size(); k++){    chg_extrema[k] = kchg_extrema[k].sum;   }

        // ---------------------------------------------------------------------
        // normalize
        {
            const FLOATTYPE sum_avol = Utils::Kahan::sum(vol_atoms);
            const FLOATTYPE sum_achg = Utils::Kahan::sum(chg_atoms);
            const FLOATTYPE sum_echg = Utils::Kahan::sum(chg_extrema);

            if(size_t(sum_avol) != numlabels) {
                std::cerr << " \n ---> error in vol integration: " << numlabels << " != " << sum_avol << "\n";
                exit(1);
            }

            if (fabs(sum_achg - sum_func) > 0.000001 || fabs(sum_echg - sum_func) > 0.000001) {
                std::cerr << " \n ---> error in chg integration: " << sum_func << " != " << sum_achg << ", " <<
                                sum_echg << " (" << (sum_achg-sum_func) << ", " << (sum_echg-sum_func) << std::endl;
                exit(1);
            }
        }

        std::transform(vol_atoms.begin(), vol_atoms.end(), vol_atoms.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), vol_voxel));
        std::transform(vol_extrema.begin(), vol_extrema.end(), vol_extrema.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), vol_voxel));

        std::transform(chg_atoms.begin(), chg_atoms.end(), chg_atoms.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), chgDens_fileUnit2e));
        std::transform(chg_extrema.begin(), chg_extrema.end(), chg_extrema.begin(), std::bind1st(std::multiplies<FLOATTYPE>(), chgDens_fileUnit2e));

        printf("Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    printf("\n =====>> Bader analysis finished!");
    timer.EndGlobal ();
    timer.PrintAll ();

    // ---------------------------------------------------------------------
    // print summary
    // ---------------------------------------------------------------------

    size_t zeroAtoms = 0;
    FLOATTYPE vsum = Utils::Kahan::sum(vol_atoms);
    FLOATTYPE csum = Utils::Kahan::sum(chg_atoms);

    std::cout << "\t total : chg = " << csum << ", vol = " << vsum << std::endl;
    std::cout << "\t vacuum: chg = " << chg_atoms[0] <<", vol = " << vol_atoms[0] << "\n";
    for(size_t i = 1; i < vol_atoms.size(); i++) {
        std::cout << "\t atom " << i << ": chg = " << chg_atoms[i] <<", vol = " << vol_atoms[i] << "\n";

        if (fabs(vol_atoms[i]) < 0.000001 || fabs(chg_atoms[i]) < 0.000001){
            std::cout << "\t --> atom " << i << " has zero chg (" << chg_atoms[i] <<") or vol (" << vol_atoms[i] << ")\n";
            zeroAtoms ++;
        }
    }
    if (zeroAtoms > 0) {
        std::cout << " Bader analysis failed for " << zeroAtoms << " atoms!\n";
        exit(1);
    }

#ifdef USE_VTK
    m_vtkVolLabeling = this->create_vtkImagedata("volLabeling");
    m_slicer_label = new vtkVolumeSlicer();
    m_slicer_label->set_volume(m_vtkVolLabeling);
#endif
    return true;
}

/**
  *   @brief  Extract the complete topologuical graph using MSC library
  *
  *   @return success flag
  */

bool TopoMS::msc() {

    if(!m_config->do_msc)
        return false;

    Utils::print_separator();
    printf("\n =====>> Extracting Topological Graph\n\n");

    const FLOATTYPE vacthreshold = m_config->threshold_vacuum * m_metadata.volume_box();

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Creating TopologicalRegularGrid...");
        fflush(stdout);

        m_tgrid = new MSC::GridType(m_grid);
#ifdef DEBUG
        //RunMeshConsistencyChecks(m_tgrid);
#endif

        printf(" Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

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

        printf(" Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }

    // ---------------------------------------------------------------------
    // restrict the discrete gradient based on numerical integration
    // ---------------------------------------------------------------------
    if (this->m_fieldtype == FT_CHG) {

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

        printf(" Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();
    }
    // else, no restriction

    // ---------------------------------------------------------------------
    // restrict the discrete gradient only in nonvacuum regions
    // ---------------------------------------------------------------------
    if (this->m_fieldtype == FT_CHG) {

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
                if ( fabs(m_topofunc->cellValue(vid) ) > vacthreshold ) {
                    continue;
                }

                cellmask->SetLabel(vid, 0);

                MSC::GridType::TopologicalRegularGrid::AdjacentCellsIterator cviter(m_tgrid);
                for (cviter.begin(vid); cviter.valid(); cviter.advance()) {
                    cellmask->SetLabel(cviter.value(), 0);
                }
            }
        }

        MSC::CellTesterLaberInput* cellmasktester = new MSC::CellTesterLaberInput();
        cellmasktester->SetLabeling(cellmask);
        m_tgrid->SetTester(cellmasktester);

        printf(" Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();

#ifdef OUTPUT_DEBUG_FIELDS
        cellmask->OutputToIntFile("cellmask.raw");
#endif
    }
    else {
        // no restriction here
        MSC::CellTesterDefaultTrue* cellmasktester = new MSC::CellTesterDefaultTrue();
        m_tgrid->SetTester(cellmasktester);
    }

    // ---------------------------------------------------------------------
    // based on the above, create the discrete gradient
    // ---------------------------------------------------------------------
    MSC::DiscreteGradientLabeling *labeling = new MSC::DiscreteGradientLabeling(m_tgrid);
    labeling->ClearAllGradient();

    MSC::RobinsLabelingAlgorithm<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction> *robin =
            new MSC::RobinsLabelingAlgorithm<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction>(m_topofunc, m_tgrid, labeling);

    robin->compute_output();
    robin->summarize();

    // ---------------------------------------------------------------------
    // compute full MSC
    // ---------------------------------------------------------------------


    // Oct 22, 2018. Attila suggested setting the restriction to null
    // since discrete gradient has already been computed, we can remove the restriction
    // to handle some cases of mismatch between numerical and discrete
    // otherwise, we see different boundary values for some critical point pairs
    // that persistence simplification cannot remove
    m_tgrid->set_restriction (nullptr);

    {
        MSC::ThreadedTimer ltimer(1);
        ltimer.StartGlobal();

        printf(" -- Setting dimensionality of ascending manifolds...");
        fflush(stdout);

        MSC::TopologicalGradientUsingAlgorithms<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction, MSC::DiscreteGradientLabeling>
                * topo_algs = new MSC::TopologicalGradientUsingAlgorithms<MSC::GridType, MSC::TopologicalExplicitDenseMeshFunction, MSC::DiscreteGradientLabeling>(m_topofunc, m_tgrid, labeling);
        topo_algs->setAscendingManifoldDimensions();

        printf(" Done!");
        ltimer.EndGlobal();
        ltimer.PrintAll();

#ifdef DEBUG
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
    if (this->m_fieldtype == FT_CHG) {

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

            if ( fabs(n.value) < vacthreshold ) {
                continue;
            }

            MSC::Vec3d ndcoords;
            this->msc_cellid_to_gcoords(n.cellindex, ndcoords);
            kdtree_add(m_kdtree_minima, ndcoords[0], ndcoords[1], ndcoords[2], i);
            cnt++;
        }

        // now i added all the "atoms"
        // for each node find closest "atom"
        // because i dont have atoms file, simply use a perturbed version of each node

        // i save the atom positions as well as the index so that i can render later the atom position.
        vector<pair<MSC::Vec3d, INT_TYPE> > restrictnodes;

        for(size_t i = 0; i < m_atoms.size (); i++) {

            float gpos[3];
            m_metadata.world_to_grid (m_atoms[i].m_pos, gpos);

            double gdpos[3] = {gpos[0], gpos[1], gpos[2]};
            int nearestIdx = kdtree_query(m_kdtree_minima, gdpos);

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

    printf("\n =====>> Topological Graph Extracted!");
    timer.EndGlobal ();
    timer.PrintAll ();

    // ---------------------------------------------------------------------
    switch(this->m_fieldtype) {
    case FT_CHG:    extract_mgraph( persistence_val, filter_val );          break;
    case FT_POT:    extract_lpot_nbrhood_li( persistence_val, filter_val);  break;
    }
    return true;
}

/**
  *   @brief  Simplify a precomputed molecular graph
  *
  *   @param  pvalue is the persistence value
  *   @param  fvalue is the persistence value
  *   @return void
  */
void TopoMS::extract_mgraph(FLOATTYPE pvalue, FLOATTYPE fvalue) {

    assert(this->m_fieldtype == FT_CHG);

    persistence_val = pvalue;
    filter_val = fvalue;

    this->m_mscbonds.clear();

    // ----------------------------------------------------------
    Utils::print_separator();
    printf("\n -- Extracting simplified molecular graph (persistence = %f, filtered at %f)\n", persistence_val, filter_val);
    m_msc->SetSelectPersAbs(persistence_val);

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    // now do the selection to show 0-1 arcs
    MSC::MSCSelectorLivingNodes<MSC::MSCType> *s_nodes = new MSC::MSCSelectorLivingNodes<MSC::MSCType>(m_msc);

    MSC::MSCSelectorNodeIndex<MSC::MSCType> *s_1saddles = new MSC::MSCSelectorNodeIndex<MSC::MSCType>(m_msc, 1);
    s_1saddles->add_parent(s_nodes);

    // Oct 22, 2018. Harsh commented this filter
    // since it seems to obstract extracting bonds when two atoms are bonded in two directions
    #if 0
      MSC::MSCSelectorRepresentative1Saddle<MSC::MSCType> *s_r1saddles = new MSC::MSCSelectorRepresentative1Saddle<MSC::MSCType>(m_msc);
      s_r1saddles->add_parent(s_1saddles);
    #else
      auto s_r1saddles = s_1saddles;
    #endif
    s_r1saddles->compute_output();

    // here is how to actually extract arcs
    MSC::StreamlineIntegratorType *sintegrator = new MSC::StreamlineIntegratorType(m_grid, m_gridfunc, m_config->threshold_error, m_config->threshold_grad, m_config->numiter);

    MSC::TerminateNearExtrema *extremumtermination = new MSC::TerminateNearExtrema(m_integrator->GetExtrema(), m_grid);
    sintegrator->SetAdvectionChecker(extremumtermination);

    // ----------------------------------------------------------
    // compute paths
    size_t pfixes = 0;
    static const float eps = 0.000001;
    for (auto it = s_r1saddles->output.begin(); it != s_r1saddles->output.end(); it++) {

        bool debug = false;

        MSCBond bond (*it);

        int sdim;
        INDEX_TYPE scidx;
        this->msc_get_node(bond.saddle, sdim, scidx, bond.scoords);

        if (debug)
            printf("\n\t > looking at saddle %d\n", *it);

        // iterate over arcs around the saddle
        MSC::MSCIteratorSurroundingArcs<MSC::MSCType> sit(m_msc);
        for (sit.begin(*it); sit.valid(); sit.advance()) {

            INT_TYPE arcid = sit.value();
            MSC::arc<FLOATTYPE>& arc = m_msc->getArc(arcid);

            // only look at saddle-minimum arcs
            if (arc.upper != *it) {
                continue;
            }

            INDEX_TYPE lowercellid = m_msc->getNode(arc.lower).cellindex;       // combinatorial says this is where we should end up
            INDEX_TYPE uppercellid = m_msc->getNode(arc.upper).cellindex;       // combinatorial says this is where we should end up

            if (debug){

                size_t i1, i2;
                int d1, d2;
                INDEX_TYPE c1, c2;
                MSC::Vec3d n1, n2;

                msc_get_node(arc.lower, d1, c1, n1);
                msc_get_node(arc.upper, d2, c2, n2);

                //bool TopoMS::msc_get_node(size_t idx, int &dim, INDEX_TYPE &cellIdx, MSC::Vec3d &ncoords) const {
                printf(" - arc: %d ([%f %f %f] %d %f) -> %d ([%f %f %f] %d %f)\n",
                       arc.lower, n1[0],n1[1],n1[2], lowercellid, m_topofunc->cellValue(lowercellid),
                       arc.upper, n2[0],n2[1],n2[2], uppercellid, m_topofunc->cellValue(uppercellid));
            }
            // filtering
            if (fabs(m_topofunc->cellValue(lowercellid)) <= filter_val ||
                fabs(m_topofunc->cellValue(uppercellid)) <= filter_val){
                continue;
            }


#if 1
            // 05/13/2018
            // Harsh put in a failsafe here
            // because, somehow, persistence simplification is not working properly

            FLOATTYPE pval = fabs(m_topofunc->cellValue(uppercellid)-m_topofunc->cellValue(lowercellid));
            if (pval < pvalue) {
                pfixes ++;
                //continue;
                printf(" -- WARNING: Explicitely filtered a small persistence value (%f) for nodes (%d: %f) and (%d: %f)\n",
                          pval,
                          arc.lower, m_topofunc->cellValue(lowercellid),
                          arc.upper, m_topofunc->cellValue(uppercellid));
                continue;
            }
#endif


            INDEX_TYPE lowervertid = m_tgrid->VertexNumberFromCellID(lowercellid);  // get vertex id in grid coordinates from combinatorial min from topological coordinates

            vector<INDEX_TYPE> geom_combinatorial;
            m_msc->fillArcGeometry(arcid, geom_combinatorial);              // combinatorial arc geometry

            std::vector<MSC::Vec3d> geom_numerical;                        // structure to hold numerical path

            // now try integrating down starting from combinatorial geomery
            for (size_t j = 0; j < geom_combinatorial.size(); j++) {

                MSC::Vec3d cf;
                this->msc_cellid_to_gcoords(geom_combinatorial[j], cf);

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


            // for this saddle, find the atom it bonds to!
            // use the kdtree to find nearest atom to saddle's minima
            {
                int ndim;
                INDEX_TYPE ncidx;
                MSC::Vec3d coord;
                this->msc_get_node(arc.lower, ndim, ncidx, coord);

                double mpos[3] = {coord[0], coord[1], coord[2]};
                int nearestAtomIdx = kdtree_query(m_kdtree_atoms, mpos);

                bond.paths.push_back(geom_numerical);
                bond.extrema.push_back(arc.lower);
                bond.ecoords.push_back(coord);
                bond.atomIds.push_back(nearestAtomIdx);

                if (debug)
                    printf(" - adding extrema %d and atom %d, path = %d\n", arc.lower, nearestAtomIdx, geom_numerical.size());
            }
        }

        // Harsh: used this explicit condition for proper definition of molecular graph
        if (bond.extrema.size() == 2){
            this->m_mscbonds[*it] = bond;
        }
        else {
          //printf(" ignoring a bond with %d extrema!\n", bond.extrema.size());
          //exit(1);
        }
    }

    // count number of critical points
    size_t counts[4] = {0,0,0,0};
    for (auto it = s_nodes->output.begin(); it != s_nodes->output.end(); it++) {

        MSC::node<FLOATTYPE>& node = m_msc->getNode(*it);
        counts[node.dim]++;
    }

    printf(" -- Simplified topology extracted! Found %d bonds!", this->m_mscbonds.size());

    timer.EndGlobal ();
    timer.PrintAll ();
    if (pfixes > 0)
        printf("   -- WARNING: explicitely filtered %d connections for persistence %f\n", pfixes, pvalue);

    printf("   -- # crit pts = [%d, %d, %d, %d] = %d\n", counts[0], counts[1], counts[2], counts[3],
                                                         counts[0] + counts[1] + counts[2] + counts[3] );

    Utils::print_separator();
    printf("\n -- Parameterizing MSC bonds...");
    fflush(stdout);
    for(auto iter = this->m_mscbonds.begin(); iter != this->m_mscbonds.end(); iter++)
        iter->second.parameterize(m_metadata);
    printf(" Done!\n");

    // use vtk to compute cross sectional areas for Bader regions
    this->analyze_bonds();
    Utils::print_separator();
}

void TopoMS::extract_lpot_nbrhood(FLOATTYPE pvalue, FLOATTYPE fvalue, float gpos[], unsigned cp_idx) {

    assert(this->m_fieldtype == FT_POT);
    assert(cp_idx == 0 || cp_idx == 3);

    if (cp_idx != 0 && cp_idx != 3) {
        printf(" extract_lpot_nbrhood works only for cp_idx = 0 or cp_idx = 3\n");
        return;
    }

    bool debug = false;

    persistence_val = pvalue;
    filter_val = fvalue;

    // --------------------------------------------------------------------------------
    // if we want to look at the maximum, we are interested in 3--2 arcs
    // otherwise, for minimum, we find 0-1 arcs!
    // the arcs we are interested in, are therefore, 3--2 arcs (max--2saddle)
    unsigned arc_idx = (cp_idx == 3) ? 2 : 0;

    // if we negated the function, we need to negate the indices as well
    if (this->m_negated) {
        cp_idx = 3-cp_idx;
        arc_idx = 2-arc_idx;
    }

    // --------------------------------------------------------------------------------
    printf("\n -- Extracting simplified potential neighborhood of a cp (%d) with outgoing arcs (%d--%d) lithium (persistence = %f, filtered at %f)\n",
           cp_idx, arc_idx, arc_idx+1, persistence_val, filter_val);
    m_msc->SetSelectPersAbs(persistence_val);

    MSC::ThreadedTimer timer(1);
    timer.StartGlobal();

    this->m_mscbonds.clear();

    // select all nodes living at current level of simplification
    MSC::MSCSelectorLivingNodes<MSC::MSCType> *s_nodes = new MSC::MSCSelectorLivingNodes<MSC::MSCType>(m_msc);
    s_nodes->compute_output();
    if(debug)
        printf(" have %d living nodes!\n", s_nodes->output.size());

    // --------------------------------------------------------------------------------

    // select all critical points of required type
    MSC::MSCSelectorNodeIndex<MSC::MSCType> *s_cpIdx = new MSC::MSCSelectorNodeIndex<MSC::MSCType>(m_msc, cp_idx);
    s_cpIdx->add_parent(s_nodes);
    s_cpIdx->compute_output();
    if(debug)
        printf(" have %d nodes of indx %d!\n", s_cpIdx->output.size(), cp_idx);

    // select the nearest critical point to lithium
    MSC::MSCSelectorNearestNode<MSC::MSCType> *s_cpAtom = new MSC::MSCSelectorNearestNode<MSC::MSCType>(m_msc, MSC::Vec3d(gpos[0], gpos[1], gpos[2]));
    s_cpAtom->add_parent(s_cpIdx);
    s_cpAtom->compute_output();
    if(debug)
        printf(" have %d nearest nodes!\n", s_cpAtom->output.size());

    // should return a single node
    assert(s_cpAtom->output.size() == 1);

    INDEX_TYPE nearest_node = *s_cpAtom->output.begin();
    int sdim;
    INDEX_TYPE scidx;
    MSC::Vec3d ncoords;
    this->msc_get_node(nearest_node, sdim, scidx, ncoords);

    if (this->m_negated)
        sdim = 3-sdim;

    if(debug)
        printf(" nearest node is %d, dim = %d, cidx = %d, pos = (%f %f %f)\n", nearest_node, sdim, scidx, ncoords[0], ncoords[1], ncoords[2]);
    assert(sdim == cp_idx);

    // collect all arcs and nodes within 2 nbrhood
    // (atom -- sad) and (sad -- XX)
    std::set<INDEX_TYPE> arcs;

    // collect Li-saddle arcs
    MSC::MSCSelectorIncidentArcs2<MSC::MSCType> *s_arcs01 = new MSC::MSCSelectorIncidentArcs2<MSC::MSCType>(m_msc, arc_idx);
    s_arcs01->add_parent(s_cpAtom);
    s_arcs01->compute_output();
    if(debug)
        printf(" have %d incident arcs of type (%d -- %d)!\n", s_arcs01->output.size(), arc_idx, 1+arc_idx);

#if 0
    // collect higher critical point of all these arcs (i.e., all 1-saddles connected to the initial minimum)
    MSC::MSCSelectorIncidentNodes<MSC::MSCType> *s_nodes01 = new MSC::MSCSelectorIncidentNodes<MSC::MSCType>(m_msc, 1);
    s_nodes01->add_parent(s_arcs01);
    s_nodes01->compute_output();
    printf(" have %d incident nodes!\n", s_nodes01->output.size());

    // collect 1-0 arcs
    MSC::MSCSelectorIncidentArcs<MSC::MSCType> *s_arcs12 = new MSC::MSCSelectorIncidentArcs<MSC::MSCType>(m_msc, 0);
    s_arcs12->add_parent(s_nodes01);
    s_arcs12->compute_output();
    printf(" have %d incident arcs!\n", s_arcs12->output.size());

    //MSC::MSCSelectorIncidentNodes<MSC::MSCType> *s_nodes12 = new MSC::MSCSelectorIncidentNodes<MSC::MSCType>(m_msc);
    //s_nodes12->add_parent(s_arcs12);
    //s_nodes12->compute_output();
    //printf(" have %d incident nodes!\n", s_nodes12->output.size());

#endif
    // collect all of these arcs in a single set
    for(auto it = s_arcs01->output.begin(); it != s_arcs01->output.end(); it++) {     arcs.insert(*it);     }
    //for(auto it = s_arcs12->output.begin(); it != s_arcs12->output.end(); it++) {     arcs.insert(*it);     }
    //for(auto iter = s_nodes01->output.begin(); iter != s_nodes01->output.end(); iter++) {   nodes.insert(*iter);    }
    //for(auto iter = s_nodes12->output.begin(); iter != s_nodes12->output.end(); iter++) {   nodes.insert(*iter);    }

    // now, for create bond paths
    for(auto it = arcs.begin(); it != arcs.end(); it++) {

        INT_TYPE arcid = *it;
        MSC::arc<FLOATTYPE>& arc = m_msc->getArc(arcid);
        //printf(" found an arc: (%d %d) %d %f\n", arc.lower, arc.upper, arc.dim, arc.persistence);

        INDEX_TYPE lowercellid = m_msc->getNode(arc.lower).cellindex;       // combinatorial says this is where we should end up
        INDEX_TYPE uppercellid = m_msc->getNode(arc.upper).cellindex;       // combinatorial says this is where we should end up
        if (fabs(m_topofunc->cellValue(lowercellid)) <= filter_val || fabs(m_topofunc->cellValue(uppercellid)) <= filter_val){
            continue;
        }

        if (this->m_fieldtype == FT_POT && arc.lower != nearest_node) {
            printf(" incorrect arc: (%d %d) %d %f\n", arc.lower, arc.upper, arc.dim, arc.persistence);
        }

        vector<INDEX_TYPE> geom_combinatorial;
        m_msc->fillArcGeometry(arcid, geom_combinatorial);              // combinatorial arc geometry

        std::vector<MSC::Vec3d> points;
        points.resize(geom_combinatorial.size());
        for(int k = 0; k < geom_combinatorial.size(); k++) {
            this->msc_cellid_to_gcoords(geom_combinatorial[k], points[k]);
        }

        // i am using this datastructure temporarily
        // these do not represent bonds, but for now, i am using them to avoid writing new ones
        // create temp bonds
        MSCBond bond (arc.lower);   // check!
        bond.scoords = ncoords;
        bond.extrema.push_back(arc.lower);
        bond.extrema.push_back(arc.upper);
        bond.paths.push_back(points);
        this->m_mscbonds[this->m_mscbonds.size()] = bond;
    }

    timer.EndGlobal ();
    timer.PrintAll ();
    Utils::print_separator();
}

void TopoMS::extract_lpot_nbrhood_li(FLOATTYPE pvalue, FLOATTYPE fvalue) {

/*
    {   //OLDER tool
        /* --------------------------------------------------------------------------------
        // Potential field
        // H and P are minima, all others are maxima (in original, non-negated field)
        * /

        MSC::MSCSelectorNodeIndex<MSC::MSCType> *s_cp0 = new MSC::MSCSelectorNodeIndex<MSC::MSCType>(m_msc, 0);
        MSC::MSCSelectorNodeIndex<MSC::MSCType> *s_cp3 = new MSC::MSCSelectorNodeIndex<MSC::MSCType>(m_msc, 3);
        s_cp0->add_parent(s_nodes);
        s_cp3->add_parent(s_nodes);

        double tpos[3];
        for(int i = 0; i < this->m_atoms.size(); i++) {

            m_atoms[i].print();
            m_metadata.world_to_grid (m_atoms[i].m_pos, tpos);

            MSC::MSCSelectorNearestNode<MSC::MSCType> *s_cpAtom = new MSC::MSCSelectorNearestNode<MSC::MSCType>(m_msc, MSC::Vec3d(tpos[0], tpos[1], tpos[2]));
            s_cpAtom->add_parent(s_cp0);
            s_cpAtom->add_parent(s_cp3);
            s_cpAtom->compute_output();

            assert(s_cpAtom->output.size() == 1);

            INDEX_TYPE nearest_node = *s_cpAtom->output.begin();
            int sdim;
            INDEX_TYPE scidx;
            MSC::Vec3d ncoords;
            this->get_msc_node(nearest_node, sdim, scidx, ncoords);

            if (this->m_negated)
                sdim = 3-sdim;

            printf(" nearest node to atom %d %s is %d, dim = %d\n", i+1, m_atoms[i].m_symbol.c_str(), nearest_node, sdim);
        }
    }
*/
    // --------------------------------------------------------------------------------
    // code to analyze lithium's neighborhood

    // Li should show up as a maximum (idx 3) in a potential field
    // the arcs we are interested in, are therefore, 3--2 arcs (max--2saddle)
    unsigned int cp_idx = 3;

    float gpos[3];
    // POT

    if (this->m_fieldtype == FT_POT) {
        unsigned atomidx = m_atoms.size()-1;    // last atom!
        m_metadata.world_to_grid (m_atoms[atomidx].m_pos, gpos);
        std::cout << " WARNING: picking up the element " << int(atomidx) << "'s position ("<<gpos[0]<<","<<gpos[1]<<","<<gpos[2]<<")!\n";
        this->m_atoms[atomidx].print();
    }
    // sans-L
    else {
        gpos[0] = 82.7324;
        gpos[1] = 188.419;
        gpos[2] = 26.0935;
        std::cout << " WARNING: picking up hardcoded position ("<<gpos[0]<<","<<gpos[1]<<","<<gpos[2]<<")!\n";
    }

    extract_lpot_nbrhood( pvalue, fvalue, gpos, cp_idx);
}

// -----------------------------------------------------------------------------------------
// I/O functions
// -----------------------------------------------------------------------------------------
/**
  *   @brief  Write value along arcs
  */
void TopoMS::write_msc_bond_stats(const std::string &filename) const {

    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    FILE *outfile = fopen(filename.c_str(), "w");
    if(!outfile){
        printf(" could not open file\n");
        return;
    }

    size_t scnt = 0;
    for(auto iter = this->m_mscbonds.begin(); iter != this->m_mscbonds.end(); iter++) {

        const MSCBond &bond = iter->second;

        if(!bond.check2())
            continue;

        float wpos[3];
        float gpos[3] = {bond.scoords[0],bond.scoords[1],bond.scoords[2]};
        m_metadata.grid_to_world(gpos, wpos);

        fprintf(outfile, "\nbond_critical_point [id = %d] at (%f %f %f)\n", ++scnt, wpos[0], wpos[1], wpos[2]);

        for(unsigned k = 0; k < bond.atomIds.size(); k++) {
            const size_t &a = bond.atomIds[k];
            fprintf(outfile, "   atom [id = %d] %s (%f %f %f)\n", a, m_atoms[a-1].m_symbol.c_str(),
                        m_atoms[a-1].m_pos[0], m_atoms[a-1].m_pos[1], m_atoms[a-1].m_pos[2]);
        }

        fprintf(outfile, "   iarea %f, ichg %f\n", bond.iarea, bond.ichg);

        const size_t nparams = bond.parameterization.size();
        for(size_t i = 0; i < nparams; i++){

            const double &x = bond.parameterization[i].first;
            const double &y = bond.bchg[i];
            fprintf(outfile, "% 08.5f, %.6E\n", x,y);
            /*const double ic = (bond.bichg.empty())  ? 0 : bond.bichg[i];
            const double ia = (bond.biarea.empty()) ? 0 : bond.biarea[i];
            fprintf(outfile, "% 08.5f, %.6E, %.6E, %.6E\n", x,y,ic,ia);*/
        }
    }

    fclose(outfile);
    printf(" Done!\n");
}

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
    fprintf(outfile, " %4s %11s %11s %11s %11s %13s\n", "#", "X", "Y", "Z", "CHARGE", "VOLUME");
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
    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    FILE *outfile = fopen(filename.c_str(), "w");
    if(!outfile){
        printf(" could not open file\n");
        return;
    }

    fprintf(outfile, " %4s %11s %11s %11s %13s %8s\n", "#", "X", "Y", "Z", "CHARGE", "ATOM");
    fprintf(outfile, "---------------------------------------------------------------------\n");

    const std::vector<INDEX_TYPE> &extrema = this->bader_get_extrema();

    int c = 0;
    for(size_t i = 0; i < extrema.size (); i++){

        if (filter_small && chg_extrema[i] < 0.000001)
            continue;

        INDEX_TYPE maxCellIdx = extrema.at(i);

        float worldPos[3], gridPos[3];

        m_metadata.idx_to_grid(maxCellIdx, gridPos);
        m_metadata.grid_to_world(gridPos, worldPos);

        fprintf(outfile, " %4d %11.7f %11.7f %11.7f %13.7f %8d\n",
                ++c, worldPos[0], worldPos[1], worldPos[2],
                chg_extrema[i], extrema2atoms.at(i));
    }
    fprintf(outfile, "---------------------------------------------------------------------\n");
    fclose(outfile);
    printf(" Done!\n");
}
