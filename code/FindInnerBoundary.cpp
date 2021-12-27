#include <iostream>
#include "FindInnerBoundary.h"
#include "MeshDoctor.h"

void FindInnerBoundary::SetFaceLabelMap(const std::string& fileName)
{
	//��ȡmesh����Ƭ����
	//int faceCount = originalMesh.faces.
	//��ȡ�������label����
	std::ifstream ifile(fileName);
	//���ļ����뵽ostringstream����buf��
	//ostringstream buf;
	char ch[10];
	int rowCount = 0;
	while (!ifile.eof())
	{
		ifile.getline(ch, 10);
		faceLabelHashmap[rowCount++] = atoi(ch);;
	}
	ifile.close();
}


void FindInnerBoundary::TraverseFaces()
{

	//MeshDoctor doctor(originalMesh);
	//doctor.debugFaceColor("D:/2879.ply", faceLabelHashmap);//��ʾһ���ض�����Ƭ����ɫ



	//auto f_end = originalMesh.faces_end();
	//for (auto f_it = originalMesh.faces_begin(); f_it != f_end; ++f_it)
	//{
	//	//fhΪ��handle
	//	FaceHandle fh = *f_it;
	//	int i = fh.idx();
	//	//����fh��ʱ�������������
	//	VertexHandle vhs[3];
	//	int index = 0;
	//	auto fv_ccw_it = originalMesh.fv_ccwiter(fh);
	//	for (; fv_ccw_it.is_valid(); ++fv_ccw_it)
	//	{
	//		//��ʱ�����������ÿ������
	//		vhs[index] = *fv_ccw_it;
	//		index++;
	//	}
	//	vhs[0] = *fv_ccw_it;
	//	vhs[1] = *(++fv_ccw_it);
	//	vhs[2] = *(++fv_ccw_it);
	//	//Ҳ���Ը�������ʸ���ıߺͰ��
	//	auto fe_it = originalMesh.fe_ccwiter(fh);
	//	auto fh_it = originalMesh.fh_ccwiter(fh);
	//	for (; fh_it.is_valid(); ++fh_it)
	//	{
	//		HalfedgeHandle hh = *fh_it;
	//	}
	//}

	FindArbitraryBoundaryFace(15);

}

/*�������ݵ�ID, Ѱ������һ���߽���Ƭ,���ر߽���Ƭ������*/
/*�㷨˼·: �Ӹ�����label��ǩ��������,�ҵ� ��һ������toothid����Ƭ,
���������Ƭ�������ڽ���*/
bool FindInnerBoundary::FindArbitraryBoundaryFace(int toothId)
{

	//1. Ѱ��map�е�һ������toothId����Ƭ
	int firstToothFaceIdx = -1;
	std::map<int, int>::const_iterator it = faceLabelHashmap.begin();
	for (; it != faceLabelHashmap.end(); it++)
	{
		if (toothId == it->second)
		{
			firstToothFaceIdx = it->first;
			break;
		}
	}
	if (firstToothFaceIdx == -1)
		return false;
	FaceHandle firstBoundaryFace(firstToothFaceIdx);// ��һ���������ݵ���Ƭ
	OpenMesh::VertexHandle startPtHandle;
	bool isFindStartPt = CounterClockTraverseEdges(firstBoundaryFace, startPtHandle);
	if (!isFindStartPt)
		return false;

	bool isFindBoundaryEdge = FindBoundaryEdgeFromAFace(firstBoundaryFace, startPtHandle, toothId);
	if (!isFindBoundaryEdge) //��һ��һ�����ڱ߽�,������Ǿͳ�����
	{
		return false;
	}
	//��¼��һ��ѹ��vector�ĵ�,�Ա������һ�����ж�
	BoundaryInfo firstBoundaryInfo = innerBoundary.at(0);
	PrintBoundaryInfo(0);
	TraverseFacesByEndVertex(firstBoundaryInfo.boundaryVH, firstBoundaryFace.idx(), toothId);
	PrintBoundaryInfo(1);
	BoundaryInfo& lastBoundaryInfo = innerBoundary.back();
	int kkk = 2;
	while (firstBoundaryInfo.boundaryVH.idx() != lastBoundaryInfo.boundaryVH.idx())
	{

		TraverseFacesByEndVertex(lastBoundaryInfo.boundaryVH, lastBoundaryInfo.boundaryFH.idx(), toothId);
		PrintBoundaryInfo(kkk++);
		lastBoundaryInfo = innerBoundary.back();
		//if (kkk == 30)
		//	break;
	}
	DebugPts(kkk);
	return true;
}

// param1 �ǵ�ǰ�߽��һ���˵�,param2 ������˵����ڵ���
void FindInnerBoundary::TraverseFacesByEndVertex(
	const OpenMesh::VertexHandle& boundaryEndVH, const int& curBoundaryFaceIdx, const int&toothId)
{
	int count = 0;
	int fddd = 0;
	for (auto vf = originalMesh.vf_begin(boundaryEndVH); vf != originalMesh.vf_end(boundaryEndVH); ++vf)
	{
		count = count + 1;
		fddd = (*vf).idx();
	}

	//���������������������,�ҵ���һ���߽���������Ƭ,Ȼ��ݹ���ú���
	for (auto vf = originalMesh.vf_begin(boundaryEndVH); vf != originalMesh.vf_end(boundaryEndVH); ++vf)
	{
		FaceHandle tranverseFace = *vf; //���������������һ����Ƭ
		int tranverseFaceIdx = tranverseFace.idx();//��ȡ�����Ƭ��mesh�����е���Ƭ����id
		////������������Ƭ boundaryFH, ��Ϊ�����Ƭ���Ѿ�֪���Ǳ߽���Ƭ��
		//if (tranverseFaceIdx == curBoundaryFaceIdx) //�����Ӧ��ȥ��,��Ϊ�п���ͬһ�����������������Ǳ߽�
		//	continue;
		//������������Ƭlabel�͵�ǰlabel��һ������Ƭ, ��Ҫ�ҵ��Ǳ߽���Ƭ,�߽���Ƭ������id��������֪����Ƭ������һֱ��
		if (faceLabelHashmap[curBoundaryFaceIdx] != faceLabelHashmap[tranverseFaceIdx])
			continue;
		//��ʣ�µ���������Ƭface��Ѱ��,�ĸ��Ǳ߽�������
		//�������ν�����ѡ,�������������������,��һ���ߺ���Ա߷�����������ǩlabel�����ݼ��Ǳ߽���Ƭ
		bool isFindBoundaryEdge = FindBoundaryEdgeFromAFace(tranverseFace, boundaryEndVH,toothId);//��ʱ��tranverseFace�ǿ�����Ŀ��߽���Ƭ��,�ٴν���ϸ��ɸѡ
		if (isFindBoundaryEdge)
			break;
	}
}

bool FindInnerBoundary::FindBoundaryEdgeFromAFace(const OpenMesh::FaceHandle&givenFH, 
	const OpenMesh::VertexHandle& givenBoundaryEndVH, const int& toothId)
{
	//OpenMesh::ArrayItems::Face needTestedFace = originalMesh.face(givenFH);
	//int curFaceLabel = faceLabelHashmap[givenFH.idx()]; //��map���ҵ������Ƭ��Ӧ��label
	bool isFindInnerBoundary = false;//�Ƿ��ҵ��߽�
	//������ǰ��Ƭ��������,���Ƿ����ҵ�����"�ڱ߽�"���Ǹ���
	auto faceHalfEdge_it = originalMesh.fh_ccwiter(givenFH);//��ʱ�������Ƭ, faceHalfEdge_it ��Ƭ�İ�ߵ�����
	for (; faceHalfEdge_it.is_valid(); ++faceHalfEdge_it)
	{
		//�жϵ�ǰ��ߵ�����ǲ��Ǹ����ĵ�, �����������
		if (originalMesh.from_vertex_handle(*faceHalfEdge_it).idx() != givenBoundaryEndVH.idx())
			continue;


		//��õ�ǰ�����ߵĶԱ�(�Ա�����������Ƭ)
		HalfedgeHandle opHalfEdgeHandle = originalMesh.opposite_halfedge_handle(*faceHalfEdge_it);
		//��ȡ������handle
		FaceHandle opFaceHandle = originalMesh.face_handle(opHalfEdgeHandle);
		//�ų��Ա��Ǳ߽�����
		if (!opFaceHandle.is_valid()) //�����ж��Ǳ߽�
			continue;
		//�������Ա�������Ƭ�ı�ǩ�͸���label��ǩһ��,����,Ҫ�ұ߽�. �߽��������Ƭlabel��Ϣ�ǲ���ͬ��
		if (faceLabelHashmap[opFaceHandle.idx()] == toothId) //��ǰ��Ƭ��������Ƭͬ��һ����ǩ,��������
			continue;
		isFindInnerBoundary = true; //�����е�����,˵����ʱ�ı߾����ڱ߽�ı���
		//�����ߵ��յ�
		VertexHandle endVH = originalMesh.to_vertex_handle(*faceHalfEdge_it);
		OpenMesh::Vec3f point = originalMesh.point(endVH);
		BoundaryInfo boundaryPt(endVH, point, givenFH);
		innerBoundary.push_back(boundaryPt); //�洢�յ�

	}
	return isFindInnerBoundary;
}

//����һ����Ƭ,��ʱ����������Ƭ,�ҵ�������Ϊ�߽����ʼ��, д��retVH
bool FindInnerBoundary::CounterClockTraverseEdges(const OpenMesh::FaceHandle&givenFH, OpenMesh::VertexHandle&retVH)
{
	//��ʱ����������Ƭ,�ҵ�һ�����,��������������ߺ���Ա�����ͬһ������,����
	auto fh_it = originalMesh.fh_ccwiter(givenFH);//��ʱ�������Ƭ
	for (; fh_it.is_valid(); ++fh_it)
	{
		HalfedgeHandle opHalfEdgeHandle = originalMesh.opposite_halfedge_handle(*fh_it);//��ǰ��ߵĶԱ�(�Ա�����������Ƭ)
		FaceHandle opFaceHandle = originalMesh.face_handle(opHalfEdgeHandle);//��ȡ������handle
		//�ų��Ա��Ǳ߽�����
		if (!opFaceHandle.is_valid()) //�����ж��Ǳ߽�
			continue;
		//�������Ա�������Ƭ�ı�ǩ�͸���label��ǩһ��,����,Ҫ�ұ߽�. �߽��������Ƭlabel��Ϣ�ǲ���ͬ��
		if (faceLabelHashmap[opFaceHandle.idx()] == faceLabelHashmap[givenFH.idx()]) //��ǰ��Ƭ��������Ƭͬ��һ����ǩ,��������
			continue;
		retVH = originalMesh.from_vertex_handle(*fh_it); //�ҵ����
		return true;
	}
	return false;
}



void FindInnerBoundary::PrintBoundaryInfo(int i)
{
	if (i < 0 || i >= innerBoundary.size())
		return;
	BoundaryInfo bi = innerBoundary.at(i);
	//std::cout << i << " "<< bi.boundaryPt[0] << " " << bi.boundaryPt[1] << " " << bi.boundaryPt[2] << std::endl;
	std::cout << bi.boundaryPt[0] << " " << bi.boundaryPt[1] << " " << bi.boundaryPt[2] << std::endl;
}

void FindInnerBoundary::DebugPts(int count)
{
	std::ofstream f("D:/boundaryPt.asc");
	if (!f.is_open())
		return;
	for(int i=0;i< count;i++)
	{
		BoundaryInfo bi = innerBoundary.at(i);
		f << bi.boundaryPt[0] << " " << bi.boundaryPt[1] << " " << bi.boundaryPt[2] << std::endl;
	}
	f.close();
}
