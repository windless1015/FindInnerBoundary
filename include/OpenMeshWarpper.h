#pragma once
#ifndef OPENMESH_WRAPPER_H
#define OPENMESH_WRAPPER_H


#ifndef M_PI_2
#define M_PI_2	1.57079632679489661923
#endif 
#pragma warning(push)
#pragma warning(disable:4244)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>


#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>

using TriMesh = OpenMesh::TriMesh_ArrayKernelT<>;
using Triangle_mesh = OpenMesh::TriMesh_ArrayKernelT<>;
using Triangle_meshF = OpenMesh::PolyMesh_ArrayKernelT<>;
using Polygon_mesh = OpenMesh::PolyMesh_ArrayKernelT<>;
using Edge = Triangle_mesh::Edge;
using Halfedge = Triangle_mesh::Halfedge;
using FaceHandle = Triangle_mesh::FaceHandle;
using VertexHandle = Triangle_mesh::VertexHandle;
using EdgeHandle = Triangle_mesh::EdgeHandle;
using HalfedgeHandle = Triangle_mesh::HalfedgeHandle;
#pragma warning(pop)

#endif