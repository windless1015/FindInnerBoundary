#include <iostream>
#include "FindInnerBoundary.h"
#include "MeshDoctor.h"

void FindInnerBoundary::SetFaceLabelMap(const std::string& fileName)
{
	//获取mesh的面片总数
	//int faceCount = originalMesh.faces.
	//读取调整后的label数据
	std::ifstream ifile(fileName);
	//将文件读入到ostringstream对象buf中
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
	//doctor.debugFaceColor("D:/2879.ply", faceLabelHashmap);//显示一个特定的面片的颜色



	//auto f_end = originalMesh.faces_end();
	//for (auto f_it = originalMesh.faces_begin(); f_it != f_end; ++f_it)
	//{
	//	//fh为面handle
	//	FaceHandle fh = *f_it;
	//	int i = fh.idx();
	//	//根据fh逆时针访问三个顶点
	//	VertexHandle vhs[3];
	//	int index = 0;
	//	auto fv_ccw_it = originalMesh.fv_ccwiter(fh);
	//	for (; fv_ccw_it.is_valid(); ++fv_ccw_it)
	//	{
	//		//逆时针获得三角面的每个顶点
	//		vhs[index] = *fv_ccw_it;
	//		index++;
	//	}
	//	vhs[0] = *fv_ccw_it;
	//	vhs[1] = *(++fv_ccw_it);
	//	vhs[2] = *(++fv_ccw_it);
	//	//也可以根据面访问该面的边和半边
	//	auto fe_it = originalMesh.fe_ccwiter(fh);
	//	auto fh_it = originalMesh.fh_ccwiter(fh);
	//	for (; fh_it.is_valid(); ++fh_it)
	//	{
	//		HalfedgeHandle hh = *fh_it;
	//	}
	//}

	FindArbitraryBoundaryFace(15);

}

/*根据牙齿的ID, 寻找任意一个边界面片,返回边界面片的索引*/
/*算法思路: 从给定的label标签数据流中,找到 第一个属于toothid的面片,
遍历这个面片的三个邻接面*/
bool FindInnerBoundary::FindArbitraryBoundaryFace(int toothId)
{

	//1. 寻找map中第一个属于toothId的面片
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
	FaceHandle firstBoundaryFace(firstToothFaceIdx);// 第一个给定牙齿的面片
	OpenMesh::VertexHandle startPtHandle;
	bool isFindStartPt = CounterClockTraverseEdges(firstBoundaryFace, startPtHandle);
	if (!isFindStartPt)
		return false;

	bool isFindBoundaryEdge = FindBoundaryEdgeFromAFace(firstBoundaryFace, startPtHandle, toothId);
	if (!isFindBoundaryEdge) //第一个一定是内边界,如果不是就出错了
	{
		return false;
	}
	//记录第一次压入vector的点,以便与最后一次做判断
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

// param1 是当前边界的一个端点,param2 是这个端点所在的面
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

	//遍历这个顶点所有三角形,找到另一个边界三角形面片,然后递归调用函数
	for (auto vf = originalMesh.vf_begin(boundaryEndVH); vf != originalMesh.vf_end(boundaryEndVH); ++vf)
	{
		FaceHandle tranverseFace = *vf; //遍历给定点的其中一个面片
		int tranverseFaceIdx = tranverseFace.idx();//获取这个面片在mesh数据中的面片索引id
		////跳过给定的面片 boundaryFH, 因为这个面片我已经知道是边界面片了
		//if (tranverseFaceIdx == curBoundaryFaceIdx) //这个不应该去掉,因为有可能同一个三角形有两个边是边界
		//	continue;
		//跳过遍历的面片label和当前label不一样的面片, 需要找的是边界面片,边界面片的索引id和我们已知的面片索引是一直的
		if (faceLabelHashmap[curBoundaryFaceIdx] != faceLabelHashmap[tranverseFaceIdx])
			continue;
		//在剩下的三角形面片face中寻找,哪个是边界三角形
		//对三角形进行甄选,遍历这个三角形三个边,有一个边和其对边分属于两个标签label的牙齿即是边界面片
		bool isFindBoundaryEdge = FindBoundaryEdgeFromAFace(tranverseFace, boundaryEndVH,toothId);//此时的tranverseFace是可能是目标边界面片的,再次进行细致筛选
		if (isFindBoundaryEdge)
			break;
	}
}

bool FindInnerBoundary::FindBoundaryEdgeFromAFace(const OpenMesh::FaceHandle&givenFH, 
	const OpenMesh::VertexHandle& givenBoundaryEndVH, const int& toothId)
{
	//OpenMesh::ArrayItems::Face needTestedFace = originalMesh.face(givenFH);
	//int curFaceLabel = faceLabelHashmap[givenFH.idx()]; //在map中找到这个面片对应的label
	bool isFindInnerBoundary = false;//是否找到边界
	//遍历当前面片的三个边,看是否能找到属于"内边界"的那个边
	auto faceHalfEdge_it = originalMesh.fh_ccwiter(givenFH);//逆时针遍历面片, faceHalfEdge_it 面片的半边迭代器
	for (; faceHalfEdge_it.is_valid(); ++faceHalfEdge_it)
	{
		//判断当前半边的起点是不是给定的点, 如果不是跳过
		if (originalMesh.from_vertex_handle(*faceHalfEdge_it).idx() != givenBoundaryEndVH.idx())
			continue;


		//获得当前遍历边的对边(对边属于相邻面片)
		HalfedgeHandle opHalfEdgeHandle = originalMesh.opposite_halfedge_handle(*faceHalfEdge_it);
		//获取相邻面handle
		FaceHandle opFaceHandle = originalMesh.face_handle(opHalfEdgeHandle);
		//排除对边是边界的情况
		if (!opFaceHandle.is_valid()) //这种判断是边界
			continue;
		//如果这个对边所在面片的标签和给定label标签一样,跳过,要找边界. 边界的两个面片label信息是不相同的
		if (faceLabelHashmap[opFaceHandle.idx()] == toothId) //当前面片和相邻面片同属一个标签,跳出继续
			continue;
		isFindInnerBoundary = true; //能运行到这里,说明此时的边就是内边界的边了
		//存入半边的终点
		VertexHandle endVH = originalMesh.to_vertex_handle(*faceHalfEdge_it);
		OpenMesh::Vec3f point = originalMesh.point(endVH);
		BoundaryInfo boundaryPt(endVH, point, givenFH);
		innerBoundary.push_back(boundaryPt); //存储终点

	}
	return isFindInnerBoundary;
}

//给定一个面片,逆时针遍历这个面片,找到可以作为边界的起始点, 写入retVH
bool FindInnerBoundary::CounterClockTraverseEdges(const OpenMesh::FaceHandle&givenFH, OpenMesh::VertexHandle&retVH)
{
	//逆时针遍历这个面片,找到一个起点,如果遍历的这个半边和其对边属于同一个牙齿,跳过
	auto fh_it = originalMesh.fh_ccwiter(givenFH);//逆时针遍历面片
	for (; fh_it.is_valid(); ++fh_it)
	{
		HalfedgeHandle opHalfEdgeHandle = originalMesh.opposite_halfedge_handle(*fh_it);//当前半边的对边(对边属于相邻面片)
		FaceHandle opFaceHandle = originalMesh.face_handle(opHalfEdgeHandle);//获取相邻面handle
		//排除对边是边界的情况
		if (!opFaceHandle.is_valid()) //这种判断是边界
			continue;
		//如果这个对边所在面片的标签和给定label标签一样,跳过,要找边界. 边界的两个面片label信息是不相同的
		if (faceLabelHashmap[opFaceHandle.idx()] == faceLabelHashmap[givenFH.idx()]) //当前面片和相邻面片同属一个标签,跳出继续
			continue;
		retVH = originalMesh.from_vertex_handle(*fh_it); //找到起点
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
