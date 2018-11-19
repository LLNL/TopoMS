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
 *  @file    TopoMS-vtk.cpp
 *  @author  Harsh Bhatia (hbhatia@llnl.gov)
 *  @date    10/01/2017
 *
 *  @brief This file provides the vtk functionality for TopoMS
 *
 *  @section DESCRIPTION
 *
 *  This file provides the vtk functionality for TopoMS
 *
 */

/// --------------------------------------------------------------------------------------
#include <set>
#include <string>

#include "MultilinearInterpolator.h"
#include "MSCBond.h"
#include "TopoMS.h"

// vtk headers
#ifdef USE_VTK
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkMatrix4x4.h>
#include <vtkPoints.h>
#include <vtkVertex.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkImageReslice.h>
#include <vtkMarchingCubes.h>
#include <vtkColorTransferFunction.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>

#include "vtkVolumeSlicer.h"
#endif

/// --------------------------------------------------------------------------------------
/**
  *   @brief  Write pixel labeling as volumes (VTK file)
  */
void TopoMS::write_labeling_vti(const std::string &filename, const std::string &fieldname, const std::vector<int> &labels) const {

#ifndef USE_VTK
    printf("TopoMS::write_labeling_vti. VTK not available. Writing as raw file instead.\n");

    FILE* fout = fopen(filename.c_str(), "wb");
    fwrite(labels.data(), sizeof(int), labels.size(), fout);
    fclose(fout);
#else

    printf(" - Writing %s...", filename.c_str());
    fflush(stdout);

    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

    float spacings[3];
    for(uint8_t i = 0; i < 3; i++) {
        spacings[i] = m_metadata.m_lattice.v[i][i] / float(m_config->is_periodic[i] ? m_metadata.m_grid_dims[i]-1 : m_metadata.m_grid_dims[i]);
    }

    imageData->SetDimensions(m_metadata.m_grid_dims[0], m_metadata.m_grid_dims[1], m_metadata.m_grid_dims[2]);
    imageData->SetOrigin(m_metadata.m_lattice_origin[0], m_metadata.m_lattice_origin[1], m_metadata.m_lattice_origin[2]);
    imageData->SetSpacing(spacings[0], spacings[1], spacings[2]);
    imageData->AllocateScalars(VTK_INT, 1);
    imageData->GetPointData()->GetScalars()->SetName(fieldname.c_str());

    for (size_t z = 0; z < m_metadata.m_grid_dims[2]; z++) {
    for (size_t y = 0; y < m_metadata.m_grid_dims[1]; y++) {
    for (size_t x = 0; x < m_metadata.m_grid_dims[0]; x++) {
        int* pixel = static_cast<int*>(imageData->GetScalarPointer(x,y,z));
        pixel[0] = labels[m_metadata.grid_to_idx(x,y,z)];
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
void TopoMS::write_mgraph_vtp(const std::string &filename) const {

#ifndef USE_VTK
    printf("TopoMS::write_mgraph. VTK not available\n");
#else

    printf(" - Writing Molecular Graph to %s...", filename.c_str());
    fflush(stdout);


    static const float periodic_cutoff = get_periodic_cutoff();
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
    std::set<INT_TYPE> nodes;
    for(auto iter = m_mscbonds.begin(); iter != m_mscbonds.end(); iter++) {

        const MSCBond &bond = iter->second;
        nodes.insert(bond.saddle);
        for(auto e = bond.extrema.begin(); e != bond.extrema.end(); e++)
            nodes.insert(*e);
    }

    for(auto n = nodes.begin(); n != nodes.end(); n++) {

        int ndim;
        INDEX_TYPE ncidx;
        MSC::Vec3d coord;

         // could also pick directly from the mscbond struct
        this->msc_get_node(*n, ndim, ncidx, coord);

        float gcoord[3] = {coord[0], coord[1], coord[2]};
        float wcoord[3] = {0,0,0};
        this->m_metadata.grid_to_world(gcoord, wcoord);
        points->InsertNextPoint(wcoord[0], wcoord[1], wcoord[2]);

        vtkSmartPointer<vtkVertex> cPnt = vtkSmartPointer<vtkVertex>::New();
        cPnt->GetPointIds()->InsertNextId(points->GetNumberOfPoints()-1);
        critPts->InsertNextCell(cPnt);
        pDims->InsertNextValue( (this->m_negated) ? 3-ndim : ndim );
        //pDims->InsertNextValue( i%2 == 0 ? 3 : 2 );
        //cDims->InsertNextValue( i%2 == 0 ? 3 : 2 );
    }

    // ------------------------------------------------------------------------------
    // Create a cell array to store the lines in and add the lines to it
    for(auto iter = m_mscbonds.begin(); iter != m_mscbonds.end(); iter++) {
    for(size_t i = 0; i < iter->second.paths.size(); i++) {

        const std::vector<MSC::Vec3d> &path = iter->second.paths[i];

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

            float gcoord[3] = {path[pidx][0], path[pidx][1], path[pidx][2]};
            float wcoord[3];
            this->m_metadata.grid_to_world(gcoord, wcoord);
            points->InsertNextPoint(wcoord[0], wcoord[1], wcoord[2]);

            pDims->InsertNextValue(6);
            polyLine->GetPointIds()->InsertNextId(points->GetNumberOfPoints()-1);
        }

        lines->InsertNextCell(polyLine);
        //cDims->InsertNextValue(6);
    }
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


/// --------------------------------------------------------------------------------------
/**
  *   @brief  Create a vtkImageData
  */

#ifdef USE_VTK
vtkImageData *TopoMS::create_vtkImagedata(const double *volume, const std::string &fname) const {

    vtkImageData *vdata = vtkImageData::New();
    vdata->SetOrigin(0.0, 0.0, 0.0);
    vdata->SetSpacing(1.0,1.0,1.0);
    vdata->SetDimensions(m_metadata.m_grid_dims[0], m_metadata.m_grid_dims[1], m_metadata.m_grid_dims[2]);
    vdata->AllocateScalars(VTK_DOUBLE, 1);
    vdata->GetPointData()->GetScalars()->SetName(fname.c_str());

    for (size_t z = 0; z < m_metadata.m_grid_dims[2]; z++) {
    for (size_t y = 0; y < m_metadata.m_grid_dims[1]; y++) {
    for (size_t x = 0; x < m_metadata.m_grid_dims[0]; x++) {
        double* pixel = static_cast<double*>(vdata->GetScalarPointer(x,y,z));
        const size_t idx = m_metadata.grid_to_idx(x,y,z);
        pixel[0] = m_negated ? -1 * volume[idx] : volume[idx];
    }}}

    return vdata;
}

vtkImageData* TopoMS::create_vtkImagedata(const std::string &fname) const {

    vtkImageData *vdata = vtkImageData::New();
    vdata->SetOrigin(0.0, 0.0, 0.0);
    vdata->SetSpacing(1.0,1.0,1.0);
    vdata->SetDimensions( m_metadata.m_grid_dims[0], m_metadata.m_grid_dims[1], m_metadata.m_grid_dims[2] );

    vdata->AllocateScalars(VTK_INT, 1);
    vdata->GetPointData()->GetScalars()->SetName(fname.c_str());

    for (size_t z = 0; z < m_metadata.m_grid_dims[2]; z++) {
    for (size_t y = 0; y < m_metadata.m_grid_dims[1]; y++) {
    for (size_t x = 0; x < m_metadata.m_grid_dims[0]; x++) {
        int* pixel = static_cast<int*>(vdata->GetScalarPointer(x,y,z));
        pixel[0] = bader_get_atomLabeling(x,y,z);
    }}}
    return vdata;
}

/**
  *   @brief  Compute a 2D slice through a given point and direction
  */
void TopoMS::compute_slice(const MSC::Vec3d &origin, const std::vector<MSC::Vec3d> &nbrs, vtkVolumeSlicer *slicer) {

    if (slicer == nullptr) {
        printf(" TopoMS::compute_slice: null slicer!\n");
        exit(1);
    }

    // compute the normal vector (along the line l-o-r)
    MSC::Vec3d mz = nbrs[1]-nbrs[0];
    mz *= 1.0/mz.Mag();

    // find the smallest of the three components
    uint8_t m = fabs(mz[0]) < fabs(mz[1]) ? 0 : 1;
            m = fabs(mz[m]) < fabs(mz[2]) ? m : 2;

    // compute tangent and cotangent
    MSC::Vec3d w (0,0,0);
    w[m] = 1;

    MSC::Vec3d mx = w.Cross(mz);
    mx *= 1.0/mx.Mag();

    MSC::Vec3d my = mz.Cross(mx);
    my *= 1.0/my.Mag();

    double nx[3] = {mx[0], mx[1], mx[2]};
    double ny[3] = {my[0], my[1], my[2]};
    double nz[3] = {mz[0], mz[1], mz[2]};
    double oo[3] = {origin[0], origin[1], origin[2]};

    slicer->set_orientation(oo, nx, ny, nz);
    slicer->compute();
}

bool TopoMS::filter_slice(vtkVolumeSlicer *slicer, const std::vector<size_t> &atomids, bool overwrite /*= false*/) const {

    // set nan for any pixel that does not have atomids as its label
    size_t n = slicer->slice()->GetNumberOfPoints();
    const int *vdims = slicer->volume()->GetDimensions();
    const size_t gdims[3] = {vdims[0], vdims[1], vdims[2]};

    double pnt2d[3], pnt3d[3];

    for(size_t i = 0; i < n; i++) {

        // this value is coming from bader function
        float value = slicer->slice()->GetPointData()->GetScalars()->GetComponent(i, 0);
        if(std::isnan(value)) {
            continue;
        }

        // transform the slice point to 3D point
        slicer->slice()->GetPoint(i, pnt2d);
        slicer->transform(pnt2d, pnt3d);

        // wrap periodic points and set the correct value
        for(uint8_t d = 0; d < 3; d++) {
            if (pnt3d[d] < 0)               pnt3d[d] += gdims[d];
            else if (pnt3d[d] >= gdims[d])  pnt3d[d] -= gdims[d];
        }

        if (pnt3d[0] < 0 || pnt3d[0] >= gdims[0] ||
            pnt3d[1] < 0 || pnt3d[1] >= gdims[1] ||
            pnt3d[2] < 0 || pnt3d[2] >= gdims[2]) {
            //slicer->slice()->GetPointData()->GetScalars()->SetComponent(i, 0, std::nan(""));
            //printf(" skipping a point~ (%f %f %f) :: (%f %f %f)\n", pnt2[0], pnt2[1], pnt2[2], pnt[0], pnt[1], pnt[2]);
            continue;
            exit(1);
        }

        // use periodic value of mscfunction for vacum filtering
        FLOATTYPE mscval = MultilinearInterpolator::trilinear_interpolation(pnt3d, m_mscfunc, gdims);

        // filter out vaccum
        int atomIdx = (fabs(mscval) <= vacthreshold_in_fileUnits) ? 0 : this->bader_get_atomLabeling(pnt3d);

        // filter out the pixels outside these atomic regions
        if (std::find(atomids.begin(), atomids.end(), atomIdx) == atomids.end()) {
            if (this->slice_labels) {
                slicer->slice()->GetPointData()->GetScalars()->SetComponent(i, 0, -1);
            }
            else {
                slicer->slice()->GetPointData()->GetScalars()->SetComponent(i, 0, std::nan(""));
            }
            continue;
        }


        // but, to write the value, use bader function
        value = MultilinearInterpolator::trilinear_interpolation(pnt3d, m_baderfunc, gdims);

        // finally, if this pixel persists, overwrite the correct value
        if (overwrite) {
            if (this->slice_labels) {
                slicer->slice()->GetPointData()->GetScalars()->SetComponent(i, 0, atomIdx);
            }
            else {
                slicer->slice()->GetPointData()->GetScalars()->SetComponent(i, 0, value);
            }
        }
    }
}

std::pair<double, double> TopoMS::integrate_slice(const vtkVolumeSlicer *slicer) const {

    // integrate the area and function on a slice
    double ifunc = 0;
    double iarea = 0;
    Utils::Kahan::KahanObject<FLOATTYPE> kchg_slice = {0};

    // now, just go over all points and integrate
    size_t n = slicer->slice()->GetNumberOfPoints();
    const int *vdims = slicer->volume()->GetDimensions();
    const size_t gdims[3] = {vdims[0], vdims[1], vdims[2]};

    double pnt2d[3], pnt3d[3];

    for(size_t i = 0; i < n; i++) {

        float value = slicer->slice()->GetPointData()->GetScalars()->GetComponent(i, 0);
        if(std::isnan(value)) {
            continue;
        }

        slicer->slice()->GetPoint(i, pnt2d);
        slicer->transform(pnt2d, pnt3d);

        // wrap periodic points and set the correct value
        for(uint8_t d = 0; d < 3; d++) {

            if (pnt3d[d] < 0)                pnt3d[d] += gdims[d];
            else if (pnt3d[d] >= gdims[d])   pnt3d[d] -= gdims[d];
        }

        if (pnt3d[0] < 0 || pnt3d[0] >= gdims[0] ||
            pnt3d[1] < 0 || pnt3d[1] >= gdims[1] ||
            pnt3d[2] < 0 || pnt3d[2] >= gdims[2]) {
            //slicer->slice()->GetPointData()->GetScalars()->SetComponent(i, 0, std::nan(""));
            //printf(" skipping a point~ (%f %f %f) :: (%f %f %f)\n", pnt2[0], pnt2[1], pnt2[2], pnt[0], pnt[1], pnt[2]);
            continue;
            exit(1);
        }

        size_t idx = this->m_metadata.grid_to_idx(int(pnt3d[0]), int(pnt3d[1]), int(pnt3d[2]));

        iarea += 1;
        kchg_slice = Utils::Kahan::KahanSum(kchg_slice, this->m_baderfunc[idx]);
    }

    ifunc = kchg_slice.sum;

    // fix the function
    if (this->m_negated)    ifunc *= -1;
    return std::pair<double, double>(iarea, ifunc);
}
#endif

/// --------------------------------------------------------------------------------------
/**
  *   @brief  Compute a 2D slice through a given saddle and parameterized bond path
  */
bool TopoMS::refresh_orthogonalSlice(int saddleNodeId, int param/*=0*/, int minp/*=-100*/, int maxp/*=100*/) {

#ifndef USE_VTK
    printf("TopoMS::refresh_orthogonalSlice. VTK not available\n");
    return false;
#else

    // identify the saddle
    int ndim;
    INDEX_TYPE ncidx;
    MSC::Vec3d ncoord;

    this->msc_get_node(saddleNodeId, ndim, ncidx, ncoord);

    // only look at 1-saddles and 2-saddles
    if(1 != ndim && 2 != ndim){
        return false;
    }

    const MSCBond &bond = this->m_mscbonds[saddleNodeId];
    if (!bond.check2()){
        return false;
    }

    MSC::Vec3d origin;
    std::vector<MSC::Vec3d> nbrs;
    nbrs.reserve(2);

    const size_t gdims[3] = {gdims[0],gdims[1],gdims[2]};

    if (param == 0) {
        bond.get_points_idx(origin, nbrs, gdims, int(-1));
    }
    else {
        float p = 0.0;
        if (param < 0)          p = -1.0 * float(param) / float(minp);
        else if (param > 0)     p = float(param) / float(maxp);

        bond.get_points(origin, nbrs, gdims, p);
    }

    vtkVolumeSlicer *cslicer = slicer();

    this->compute_slice(origin, nbrs, cslicer);
    this->filter_slice(cslicer, bond.atomIds, true);
    return true;
#endif
}

/// --------------------------------------------------------------------------------------
/**
  *   @brief  Analyze molecular bonds
  */
void TopoMS::analyze_bonds() {

    bool debug = false;

    // ---------------------------------------------------------------------
    // upon computation, the integrated areas and function values are simply "sums"
    // now, i need to convert them into physical units
    // i am using the same conversion as I do in the bader analysis
    // that is, treating each point on the slice as a full 3D voxel
    // this is only an apporximation, since depending upon the orientation of the slice,
    // this can be more or less accurate!
    // ---------------------------------------------------------------------

    const size_t gdims[3] = {this->m_metadata.m_grid_dims[0],this->m_metadata.m_grid_dims[1],this->m_metadata.m_grid_dims[2]};
    const size_t num_gridPts = m_metadata.grid_sz();

    const FLOATTYPE vol_box = m_metadata.volume_box();
    const FLOATTYPE vol_voxel = m_metadata.volume_voxel();

    const FLOATTYPE file_to_ae = m_metadata.volume_file2Angs() * m_metadata.charge_file2electrons();
    const FLOATTYPE chgDens_fileUnit2e = (this->m_inputtype == IT_CUBE ? vol_box * file_to_ae : 1.0) / (FLOATTYPE) num_gridPts;

    printf(" -- Analyzing properties of %d bonds...", this->m_mscbonds.size());
    fflush(stdout);

    if(debug)
        printf("\n");

#ifdef USE_VTK
    // we create a single slicer that we move along the bond
    vtkVolumeSlicer *slicer = new vtkVolumeSlicer();
    slicer->set_volume(this->m_vtkFunction);

    // return value of integrate function
    std::pair<double,double> res;
#endif

    // let's go over all the bonds
    for(auto iter = this->m_mscbonds.begin(); iter != this->m_mscbonds.end(); iter++) {

        MSCBond &bond = iter->second;
        const size_t nparams = bond.parameterization.size();

        if(debug)
            printf("\n saddle %d\n", bond.saddle);

#ifdef USE_VTK
        // first, analyze the mid point (the saddle)
        MSC::Vec3d origin;
        std::vector<MSC::Vec3d> nbrs;
        nbrs.reserve(2);

        bond.get_points_idx(origin, nbrs, gdims, -1);

        // compute, filter, and integrate the slice
        TopoMS::compute_slice(origin, nbrs, slicer);
        this->filter_slice(slicer, bond.atomIds);
        res = this->integrate_slice(slicer);

        bond.iarea = vol_voxel*res.first;
        bond.ichg  = chgDens_fileUnit2e*res.second;

        if(debug)
            printf(" %d : %f %f\n", bond.saddle, bond.iarea, bond.ichg);

        bond.bichg.resize(nparams);
        bond.biarea.resize(nparams);
#endif

        bond.bchg.resize(nparams);

        for(size_t i = 0; i < nparams; i++) {

            // the bond parameterized value and point
            const float x = bond.parameterization[i].first;
            const MSC::Vec3d &p = bond.parameterization[i].second;

            /*
#ifdef USE_VTK
            // get the point and neighbors for this point on the bond
            bond.get_points_idx(origin, nbrs, gdims, int(i));
                // origin will always be the same as p

            // compute, filter, and integrate the slice
            TopoMS::compute_slice(origin, nbrs, slicer);
            this->filter_slice(slicer, bond.atomIds);
            res = this->integrate_slice(slicer);

            bond.biarea[i] = vol_voxel*res.first;
            bond.bichg[i]  = chgDens_fileUnit2e*res.second;
#endif
            */

            // ok, now, also note the actual value at this point
            double pos[3] = {p[0], p[1], p[2]};
            bond.bchg[i] = MultilinearInterpolator::trilinear_interpolation(pos, m_baderfunc, gdims);
            bond.bchg[i] *= chgDens_fileUnit2e*(this->m_negated ? -1.0 : 1.0);
        }
    }
    printf(" Done!\n");
#ifdef USE_VTK
    delete slicer;
#endif
}

/// --------------------------------------------------------------------------------------
/// --------------------------------------------------------------------------------------


#if 0
void TopoMS::refreshIsosurface() {
#ifndef USE_VTK
    print Visit not available!
#else
    /*

    //printf(" LithiumApp::refreshIsosurface(%d)\n", curr_tIndex);
    float val = ui.dsb_isosurface->value();

    // compute 3D surface
    if(ui.cb_isosurface_vol->isChecked()) {
        //vtkIsoSurfaces[curr_tIndex].insert( std::pair<float, vtkPolyData*> ( val, compute_isosurface(vtkImage.at(curr_tIndex), val) ) );
        vtkIsoSurfaces[curr_tIndex] = compute_isosurface(vtkImage.at(curr_tIndex), val);
    }
/*
    // compute 2D surface
    if(ui.cb_isosurface_slice->isChecked() && vtkSlice.at(curr_tIndex) != 0)
        vtkIsoSurfaces[curr_tIndex].insert( std::pair<float, vtkPolyData*> ( val, compute_isosurface(vtkSlice.at(curr_tIndex), val) ) );
*/

    //viewer->updateGL();
#endif
}

// -----------------------------------------------------------------------------------------
// static functionality
// -----------------------------------------------------------------------------------------
#ifdef USE_VTK
vtkPolyData* TopoMS::compute_isosurface(vtkImageData *vtkdata, float value) {

    vtkSmartPointer<vtkMarchingCubes> vtkSlicer = vtkSmartPointer<vtkMarchingCubes>::New();
    vtkSlicer->SetValue(0, value);
    vtkSlicer->ComputeNormalsOn();

    vtkSlicer->SetInputData( vtkdata );
    vtkSlicer->Update();

    vtkPolyData *slice = vtkPolyData::New();
    slice->DeepCopy( vtkSlicer->GetOutput() );
    return slice;
}
#endif
#endif
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
