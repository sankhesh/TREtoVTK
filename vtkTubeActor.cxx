/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit @ UNC
  Module:    $RCSfile: vtkTubeUtils.cxx
  Language:  C++:
  Date:      $Date:
  Version:   $Revision: 1.0

=========================================================================*/

#include "vtkTubeFilterColor.h"
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkScalarsToColors.h>
#include <vtkSphereSource.h>
//#include <vtkVectors.h>
#include "../itkTube.h"
#include "../itkTubePoint.h"
#include "vtkTubeFilter.h"
#include "vtkTubeUtils.h"

/** Constructor */
vtkTubeUtils::vtkTubeUtils() { m_ViewType = 0; }

bool vtkTubeUtils::vesselActor(vtkActor *actor, vtkActor *sphere1,
                               vtkActor *sphere2, itk::Tube *tube, int nSides) {
  int i;

  if (nSides < 3)
    nSides = 3;

  if (tube->GetPoints()->size() < 2) {
    std::cout << "Not enough points" << std::endl;
    return false;
  }

  // Step 1: copy skeleton points from a vessel into vtkPoints
  // vtkpoints assumes a triplet is coming so use pointer arithmetic
  // to jump to the next spot in a multidimensional array
  std::list<itk::TubePoint *>::iterator pnt = tube->GetPoints()->begin();
  itk::Vector<float> lT(3);
  itk::Vector<float> lX(3);

  float lR;
  lX[0] = (*(*pnt)->GetX())[0];
  lX[1] = (*(*pnt)->GetX())[1];
  lX[2] = (*(*pnt)->GetX())[2];
  lR = (*pnt)->GetR();
  pnt++;
  lT[0] = lX[0] - (*(*pnt)->GetX())[0];
  lT[1] = lX[1] - (*(*pnt)->GetX())[1];
  lT[2] = lX[2] - (*(*pnt)->GetX())[2];
  // vNorm(&lT);
  lT.Get_vnl_vector().normalize();

  pnt++;
  itk::Vector<float> cT(3);
  float tf;
  for (i = 1; pnt != tube->GetPoints()->end(); pnt++) {
    cT[0] = lX[0] - (*(*pnt)->GetX())[0];
    cT[1] = lX[1] - (*(*pnt)->GetX())[1];
    cT[2] = lX[2] - (*(*pnt)->GetX())[2];
    // vNorm(&cT);
    cT.Get_vnl_vector().normalize();
    tf = dot_product(cT.Get_vnl_vector(), lT.Get_vnl_vector());
    if (fabs(tf) < 0.975 || fabs((*pnt)->GetR() - lR) > 0.025 * lR) {
      i++;
      lX[0] = (*(*pnt)->GetX())[0];
      lX[1] = (*(*pnt)->GetX())[1];
      lX[2] = (*(*pnt)->GetX())[2];
      pnt++;
      if (pnt != tube->GetPoints()->end()) {
        lT[0] = lX[0] - (*(*pnt)->GetX())[0];
        lT[1] = lX[1] - (*(*pnt)->GetX())[1];
        lT[2] = lX[2] - (*(*pnt)->GetX())[2];
        // vNorm(&lT);
        lT.Get_vnl_vector().normalize();
        lR = (*pnt)->GetR();
      } else
        pnt--;
    }
  }
  int nPoints = i + 1;
  vtkPoints *vPoints = vtkPoints::New();
  vPoints->SetNumberOfPoints(nPoints);

  vtkScalars *vScalars = vtkScalars::New();
  vScalars->SetNumberOfScalars(nPoints);

  pnt = tube->GetPoints()->begin();
  lX[0] = (*(*pnt)->GetX())[0];
  lX[1] = (*(*pnt)->GetX())[1];
  lX[2] = (*(*pnt)->GetX())[2];

  lR = (*pnt)->GetR();
  vPoints->SetPoint(0, (float)((*(*pnt)->GetX())[0]),
                    (float)((*(*pnt)->GetX())[1]),
                    (float)((*(*pnt)->GetX())[2]));
  vScalars->SetScalar(0, (*pnt)->GetR());

  //
  vtkSphereSource *sphereSource1 = vtkSphereSource::New();
  sphereSource1->SetCenter((float)((*(*pnt)->GetX())[0]),
                           (float)((*(*pnt)->GetX())[1]),
                           (float)((*(*pnt)->GetX())[2]));
  sphereSource1->SetRadius((*pnt)->GetR() * 0.95);

  vtkPolyDataMapper *sphereMapper1 = vtkPolyDataMapper::New();
  sphereMapper1->SetInput(sphereSource1->GetOutput());

  sphere1->SetMapper(sphereMapper1);

  sphereMapper1->Delete();
  sphereSource1->Delete();
  //

  pnt++;

  lT[0] = lX[0] - (*(*pnt)->GetX())[0];
  lT[1] = lX[1] - (*(*pnt)->GetX())[1];
  lT[2] = lX[2] - (*(*pnt)->GetX())[2];
  // vNorm(&lT);
  lT.Get_vnl_vector().normalize();
  pnt++;
  for (i = 1; pnt != tube->GetPoints()->end(); pnt++) {
    cT[0] = lX[0] - (*(*pnt)->GetX())[0];
    cT[1] = lX[1] - (*(*pnt)->GetX())[1];
    cT[2] = lX[2] - (*(*pnt)->GetX())[2];

    // vNorm(&cT);
    cT.Get_vnl_vector().normalize();

    // tf = vDotProd(&cT, &lT);
    tf = dot_product(cT.Get_vnl_vector(), lT.Get_vnl_vector());
    if (fabs(tf) < 0.975 || fabs((*pnt)->GetR() - lR) > 0.025 * lR) {
      vPoints->SetPoint(i, (float)((*(*pnt)->GetX())[0]),
                        (float)((*(*pnt)->GetX())[1]),
                        (float)((*(*pnt)->GetX())[2]));
      vScalars->SetScalar(i, (*pnt)->GetR());

      i++;
      lX[0] = (*(*pnt)->GetX())[0];
      lX[1] = (*(*pnt)->GetX())[1];
      lX[2] = (*(*pnt)->GetX())[2];
      pnt++;
      if (pnt != tube->GetPoints()->end()) {
        lT[0] = lX[0] - (*(*pnt)->GetX())[0];
        lT[1] = lX[1] - (*(*pnt)->GetX())[1];
        lT[2] = lX[2] - (*(*pnt)->GetX())[2];
        // vNorm(&lT);
        lT.Get_vnl_vector().normalize();
        lR = (*pnt)->GetR();
      } else
        pnt--;
    }
  }
  pnt--;
  vPoints->SetPoint(i, (float)((*(*pnt)->GetX())[0]),
                    (float)((*(*pnt)->GetX())[1]),
                    (float)((*(*pnt)->GetX())[2]));
  vScalars->SetScalar(i, (*pnt)->GetR());

  //
  vtkSphereSource *sphereSource2 = vtkSphereSource::New();
  sphereSource2->SetCenter((float)((*(*pnt)->GetX())[0]),
                           (float)((*(*pnt)->GetX())[1]),
                           (float)((*(*pnt)->GetX())[2]));
  sphereSource2->SetRadius((*pnt)->GetR() * 0.95);

  vtkPolyDataMapper *sphereMapper2 = vtkPolyDataMapper::New();
  sphereMapper2->SetInput(sphereSource2->GetOutput());

  sphere2->SetMapper(sphereMapper2);

  sphereMapper2->Delete();
  sphereSource2->Delete();
  //

  // Step 2: create a point id list (for a polyline this is just linear)
  int *pntIds = new int[nPoints];
  for (i = 0; i < nPoints; i++)
    pntIds[i] = i;

  // Step3: create a polyline from the points and pt id list
  vtkPolyLine *vPLine = vtkPolyLine::New();
  vPLine->Initialize(nPoints, pntIds, vPoints);

  // Step 4: convert this to a cellarray (required for input to polydata)
  vtkCellArray *vCA = vtkCellArray::New();
  vCA->InsertNextCell(vPLine->MakeObject());

  // Step 5: create a scalar array that indicates the radius at each
  // skeleton point. Vtk's way of setting radius is screwy: it fails if every
  // point has the same radius. It also uses a minimum radius  (called radius)
  // and a max radius (defined by scale factor). In order to get this to work,
  // you need to find the min and max of your vessel radii--if the same, later
  // set a constant radius in the tube filter. If not the same, you need to
  // define the min radius and the ratio max/min. If you send these params,
  // the tube will end up with proper radius settings. Weird.

  // Step 6: convert to polydata (required for input to anything else)
  vtkPolyData *vPData = vtkPolyData::New();
  vPData->SetLines(vCA);
  vPData->SetPoints(vPoints);
  float range[2];
  bool use_scalars = false;
  float min_scalar, max_scalar;
  vScalars->GetRange(range);
  min_scalar = range[0];
  max_scalar = range[1];
  // std::cout << "min_scalar = " << min_scalar << " max_scalar = " <<
  // max_scalar << std::endl;
  if (min_scalar <= 0.0001)
    min_scalar = 0.0001;
  if (max_scalar < min_scalar)
    max_scalar = min_scalar;
  if (min_scalar != max_scalar) {
    use_scalars = true;
    vPData->GetPointData()->SetScalars(vScalars);
  }

  // Step 7: remove any duplicate points from polydata. The tube filter
  // fails if any duplicates are present
  vtkCleanPolyData *vClean = vtkCleanPolyData::New();
  vClean->SetInput(vPData);

  // Step 8: make tubes. The number of sides per tube is set by nsides.
  // Even an nsides of 3 looks surprisingly good.
  vtkTubeFilter *vTFilter = vtkTubeFilter::New();
  vTFilter->SetNumberOfSides(nSides);
  vTFilter->SetInput(vClean->GetOutput());
  vTFilter->CappingOff();

  if (use_scalars) {
    vTFilter->SetVaryRadiusToVaryRadiusByScalar();
    vTFilter->SetRadius(min_scalar); // this call sets min rad. Weird.
    vTFilter->SetRadiusFactor(max_scalar / min_scalar); // sets max rad. Weird
  } else {
    vTFilter->SetRadius(min_scalar);
    vTFilter->SetVaryRadiusToVaryRadiusOff();
  }

  // Step 9: create a mapper of the tube
  vtkPolyDataMapper *vMapper = vtkPolyDataMapper::New();
  vtkPolyData *vPDOutput = vtkPolyData::New();

  // vPDOutput=vTFilter->GetOutput();

  // float ColorRange[2];

  // vtkLookupTable* lut = vtkLookupTable::New();
  /*lut->SetNumberOfColors(3);
  lut->Build();
  lut
  */
  /*
  lut->SetHueRange(0.6667,0.0);
  lut->SetSaturationRange(1.0,1.0);
  lut->SetValueRange(0, 256);
  lut->SetAlphaRange(1.0,1.0);
  lut->SetNumberOfColors(256);
  lut->Build();

  vMapper->SetLookupTable(lut);
*/
  vMapper->SetInput(vTFilter->GetOutput());

  vMapper->ScalarVisibilityOff(); // do not interpret scalars as color command
  // vMapper->ScalarVisibilityOn();    //interpret scalars as color command

  // vtkObject* vScalar2Color = vtkObject::New();
  // vScalar2Color->MapScalarsThroughTable (vtkScalars *scalars)
  // vtkScalarsToColors* vScalar2Color = new vtkScalarsToColors;

  // Step 10: Add the mapper to the actor. You can now delete everything.
  // A matrix for the actor, colors, opacities, etc can be set by
  // the caller before or after this function is called.
  actor->SetMapper(vMapper);

  vPoints->Delete();
  delete[] pntIds;
  vScalars->Delete();
  vPLine->Delete();
  vClean->Delete();
  vCA->Delete();
  vPData->Delete();
  vTFilter->Delete();
  vMapper->Delete();

  return true;
}
