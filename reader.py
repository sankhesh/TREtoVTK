#!/usr/bin/env python

from __future__ import print_function
from io import StringIO

import itk
import vtk

def DowncastToVesselTubeSOPoint(soPoint):
    '''Hacky way to downcast SpatialObjectPoint.'''
    buf = StringIO()
    print(soPoint, file=buf)
    buf.seek(0)
    props = buf.read().split("\n")

    dim = len(soPoint.GetPosition())
    vesselTubePoint = itk.VesselTubeSpatialObjectPoint[dim]()

    vesselTubePoint.SetID(soPoint.GetID())
    vesselTubePoint.SetPosition(*soPoint.GetPosition())
    vesselTubePoint.SetBlue(soPoint.GetBlue())
    vesselTubePoint.SetGreen(soPoint.GetGreen())
    vesselTubePoint.SetRed(soPoint.GetRed())
    vesselTubePoint.SetAlpha(soPoint.GetAlpha())

    radius = float(props[3].strip()[len("R: "):])
    vesselTubePoint.SetRadius(radius)

    tangent = list(map(float, props[5].strip()[len("T: ["):-1].split(",")))
    vesselTubePoint.SetTangent(*tangent)

    normal1 = list(map(float, props[6].strip()[len("Normal1: ["):-1].split(",")))
    normal2 = list(map(float, props[7].strip()[len("Normal2: ["):-1].split(",")))
    vesselTubePoint.SetNormal1(*normal1)
    vesselTubePoint.SetNormal2(*normal2)

    medialness = float(props[8].strip()[len("Medialness: "):])
    vesselTubePoint.SetMedialness(medialness)

    ridgeness = float(props[9].strip()[len("Ridgeness: "):])
    vesselTubePoint.SetRidgeness(ridgeness)

    alpha1 = float(props[10].strip()[len("Alpha1: "):])
    alpha2 = float(props[11].strip()[len("Alpha2: "):])
    alpha3 = float(props[12].strip()[len("Alpha3: "):])
    vesselTubePoint.SetAlpha1(alpha1)
    vesselTubePoint.SetAlpha2(alpha2)
    vesselTubePoint.SetAlpha3(alpha3)

    mark = float(props[13].strip()[len("Mark: "):])
    vesselTubePoint.SetMark(bool(mark))

    return vesselTubePoint

def GetTubePoints(tube):
    '''Gets the points and radii associated with the tube.'''
    points = list()
    for j in range(tube.GetNumberOfPoints()):
        point = tube.GetPoint(j)
        point = DowncastToVesselTubeSOPoint(point)

        radius = point.GetRadius()
        pos = point.GetPosition()

        # I think I need to extract the values otherwise corruption occurs
        # on the itkPointD3 objects.
        points.append(((pos[0], pos[1], pos[2]), radius))
    return points

def TubeIterator(tubeGroup):
    '''Iterates over all tubes in a tube group.'''
    obj = itk.down_cast(tubeGroup)
    if isinstance(obj, itk.VesselTubeSpatialObject[3]):
        yield obj

    # otherwise, `obj` is a GroupSpatialObject
    children = obj.GetChildren()
    for i in range(obj.GetNumberOfChildren()):
        for tube in TubeIterator(children[i]):
            yield tube

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print('Usage: %s TRE' % sys.argv[0])
        sys.exit(1)
    tre_file = sys.argv[1]

    reader = itk.SpatialObjectReader[3].New()
    reader.SetFileName(tre_file)
    reader.Update()

    polyData = vtk.vtkPolyData()
    radiusArray = vtk.vtkFloatArray()
    radiusArray.SetNumberOfComponents(1)
    radiusArray.SetName("Radius")
    lines = vtk.vtkCellArray()
    points = vtk.vtkPoints()

    for index, tube in enumerate(TubeIterator(reader.GetGroup())):
        tube.ComputeObjectToWorldTransform()
        transform = tube.GetIndexToWorldTransform()
        scaling = [transform.GetMatrix()(i,i) for i in range(3)]
        scale = sum(scaling) / len(scaling)

        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(tube.GetNumberOfPoints())

        linePtId = 0
        for pt, radius in GetTubePoints(tube):
            pt = transform.TransformPoint(pt)
            # columns: TUBE_ID POINT_X POINT_Y POINT_Z RADIUS
            print('%d %f %f %f %f' % (tube.GetId(), pt[0], pt[1], pt[2], radius*scale))
            ptId = points.InsertNextPoint(pt)
            line.GetPointIds().SetId(linePtId, ptId)
            radiusArray.InsertNextTuple1(radius * scale)
            linePtId += 1
        lines.InsertNextCell(line)
        break

    polyData.SetPoints(points)
    polyData.SetLines(lines)
    polyData.GetPointData().SetScalars(radiusArray)

    xmlWriter = vtk.vtkXMLPolyDataWriter()
    xmlWriter.SetFileName(tre_file + ".vtp")
    xmlWriter.SetInputData(polyData)
    #  xmlWriter.SetInputConnection(clean.GetOutputPort())
    xmlWriter.Write()
