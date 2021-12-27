#pragma once
#include <map>
#include "OpenMeshWarpper.h"
/*
本算法的输入: 1. 原始牙颌数据  2. 原始牙颌数据的面片标记数据
本算法的输出: 对于每一个牙齿的牙龈线的点
*/

class FindInnerBoundary
{
public:
	struct BoundaryInfo
	{
		BoundaryInfo(const VertexHandle& givenVH, 
			const OpenMesh::Vec3f& givenPt,
			const OpenMesh::FaceHandle& givenFH)
		{
			boundaryVH = givenVH;
			boundaryPt = givenPt;
			boundaryFH = givenFH;
		};
		VertexHandle boundaryVH;
		OpenMesh::Vec3f boundaryPt;
		OpenMesh::FaceHandle boundaryFH;
	};
public:
	FindInnerBoundary(Triangle_mesh& inputMesh):originalMesh(inputMesh){};
	~FindInnerBoundary() {};
public:
	void SetFaceLabelMap(const std::string& fileName);
	void TraverseFaces();
	bool FindArbitraryBoundaryFace(int toothId);
	void TraverseFacesByEndVertex(const OpenMesh::VertexHandle&, const int&, const int&);
	bool FindBoundaryEdgeFromAFace(const OpenMesh::FaceHandle&, const OpenMesh::VertexHandle&, const int&);
	bool CounterClockTraverseEdges(const OpenMesh::FaceHandle&, OpenMesh::VertexHandle&);
	void PrintBoundaryInfo(int i);
	void DebugPts(int count);
private:
	std::map<int, int> faceLabelHashmap;//面片号 -> 面片所属的标记
	Triangle_mesh& originalMesh;
	std::vector<BoundaryInfo> innerBoundary; //边界信息数组
};


