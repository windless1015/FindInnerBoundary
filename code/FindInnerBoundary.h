#pragma once
#include <map>
#include "OpenMeshWarpper.h"
/*
���㷨������: 1. ԭʼ�������  2. ԭʼ������ݵ���Ƭ�������
���㷨�����: ����ÿһ�����ݵ������ߵĵ�
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
	std::map<int, int> faceLabelHashmap;//��Ƭ�� -> ��Ƭ�����ı��
	Triangle_mesh& originalMesh;
	std::vector<BoundaryInfo> innerBoundary; //�߽���Ϣ����
};


