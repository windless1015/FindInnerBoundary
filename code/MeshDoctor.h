#pragma once

#include <OpenMeshWarpper.h>

enum  QualityType
{
	idMinAngle= 0x0001, //最小角
	idMaxAngle= 0x0002, //最大角
	idAspectRatio= 0x0004, //最短边与最长边比
	idSmoothVertex= 0x0008, //顶点光滑度
	idSkewnessEquiangle= 0x0020, //偏斜等角
	idAreaRatio= 0x0040,//三角形面积比
	idSkewnessEquiarea= 0x0800, //等面积比
};
enum  CurvatureType {
	idMean, //平均曲率
	idGauss,//高斯曲率
	idMax,  //最大曲率流
	idMin,  //最小曲率流
};

enum Color {
	idRed, //需要处理的单元
	idBlue //good
};

class MeshDoctor 
{
public:
	using TriMesh = Triangle_mesh;
	using Point3f = TriMesh::Point;
	using ftype = float;
public:
	MeshDoctor(TriMesh& msh);
	~MeshDoctor();
public://质量检查工具
	/*
	 * @brief 三角形最大角与最小角计算
	 * @detail 一般认为小于10度或者大于145度的三角形为不可接受三角形
	 */
	void updateAngle();
	void calcAngle(const FaceHandle& fh,ftype& minAgle,ftype& maxAngle);
	/*
	 * @brief 纵横比计算，
	 * @datail一般接近1的三角形为质量较好的，超过3的三角形为不好的三角形
	 */
	void updateAspectRatio();
	ftype calcAspectRatio(const FaceHandle& fh);

	 /*
	  * @brief 偏斜等角
	  * @dateil 一般建议给小于0.7
	  */
	void updateSkewnessEquiangle();
	ftype calcSkewnessEquiarea(const FaceHandle& fh,int cn=3);
	/*
	 * @brief 面积比
	 * @dateil 一般适用于均匀网格情况
	 */
	void updateSkewnessEquiarea();
	ftype calcSkewnessEquiangle(const FaceHandle& fh);
	/*
	 * @brief 相邻面积比
	 * @detail
	 */
	void updateAreaRatio();
	ftype calcAreaRatio(const EdgeHandle& eh);
	/*
	 * @brief 张亮投票法判别网格特征
	 */
	void updateSmoothVertex(int n=1);
	Point3f calcSmoothVertex(const VertexHandle& vh,int n = 1);
public:
	/*
	 * @brief 亏格检测
	 * @detail 检查之前必须调用deleteIrrVerts(),网格必须有边界
	 */
	void checkGenius();
public:
	//删除孤立点、非流行边、小组件
	void deleteIrrVerts(ftype esmallcomp);
	//删除孤立点、非流行边、保留最大组件
	void deleteIrrVerts();
	//删除顶点的边界光速法
	bool smoothBoundaryDelVertex();
	//网格去重
	void repeatMesh(ftype _eps = FLT_EPSILON);
	//保存纹理坐标网格去重
	void repeatMeshTextCoord(ftype _eps = FLT_EPSILON);
public:
	//去除局部的纵横比
	void nLocalRemoveAspectRatio(ftype fc = 3);


public://曲率计算工具
	static ftype calcCotanWeight(TriMesh& mesh, const EdgeHandle& eh);
	static ftype calcVoronoiArea(TriMesh& mesh, const VertexHandle& eh);
	static ftype calcVoronoiAreaBarycentric(TriMesh& mesh, const VertexHandle& eh);
	static Point3f calcLaplace(TriMesh& mesh, const VertexHandle& vh);
	static ftype calcAngleSum(TriMesh& mesh, const VertexHandle& v);
	static ftype calcVertexCurvature(TriMesh& mesh, const VertexHandle& v,
		CurvatureType ct= CurvatureType::idGauss);
	static void calcVertexCurvature(TriMesh& mesh, const VertexHandle& v,
		ftype& minvc, ftype& maxvc, ftype& gaussvc, ftype& meanvc);
public: //统一到处ply格式文件
	void debugAngle(const std::string& filename);
	void debugAspectRatio(const std::string& filename);
	void debugSkewnessEquiangle(const std::string& filename);
	void debugGenus(const std::string& filename);
	void debugFaceColor(const std::string& filename, int faceIdx);
	void debugFaceColor(const std::string& filename, std::map<int, int>& map);
protected:
	static void nRemeshKeepBoundry(Triangle_mesh& mesh, ftype avedg);
	static void extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& caids);	
protected:
	ftype calcFaceArea(const FaceHandle& fh);
protected:
	void selectFace(const FaceHandle& fh, int n, std::vector<FaceHandle>& fhs,bool issave);
	void selectVFace(const VertexHandle& fh, int n, std::vector<FaceHandle>& fhs, bool issave);
	void selectFVace(const std::vector<FaceHandle>& fh, int n, std::vector<VertexHandle>& fhs, bool issave);
	void selectVerts(const VertexHandle& fh, int n, std::vector<VertexHandle>& fhs,bool issave);
	void selectFace(const std::vector<FaceHandle>& fh, int n, std::vector<FaceHandle>& fhs, bool issave);
	void selectVerts(Triangle_mesh& mesh,
		const std::vector<int>& that,int n,
		std::vector<int>&spts, bool issave);
protected:
	TriMesh& mesh;
public:
	static OpenMesh::FPropHandleT<ftype> MinAngle; //最小角
	static OpenMesh::FPropHandleT<ftype> MaxAngle; //最大角
	static OpenMesh::FPropHandleT<ftype> AspectRatio; //纵横比
	static OpenMesh::FPropHandleT<ftype> SkewnessEquiangle; //偏斜等角
	static OpenMesh::VPropHandleT<Point3f> SmoothVertex;//特征度
	static OpenMesh::FPropHandleT<ftype> SkewnessEquiarea; //偏斜比
	static OpenMesh::EPropHandleT<ftype> AreaRatio; //偏斜比
	static OpenMesh::VPropHandleT<bool> Genius; //亏格
};