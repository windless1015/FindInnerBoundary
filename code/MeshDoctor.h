#pragma once

#include <OpenMeshWarpper.h>

enum  QualityType
{
	idMinAngle= 0x0001, //��С��
	idMaxAngle= 0x0002, //����
	idAspectRatio= 0x0004, //��̱�����߱�
	idSmoothVertex= 0x0008, //����⻬��
	idSkewnessEquiangle= 0x0020, //ƫб�Ƚ�
	idAreaRatio= 0x0040,//�����������
	idSkewnessEquiarea= 0x0800, //�������
};
enum  CurvatureType {
	idMean, //ƽ������
	idGauss,//��˹����
	idMax,  //���������
	idMin,  //��С������
};

enum Color {
	idRed, //��Ҫ����ĵ�Ԫ
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
public://������鹤��
	/*
	 * @brief ��������������С�Ǽ���
	 * @detail һ����ΪС��10�Ȼ��ߴ���145�ȵ�������Ϊ���ɽ���������
	 */
	void updateAngle();
	void calcAngle(const FaceHandle& fh,ftype& minAgle,ftype& maxAngle);
	/*
	 * @brief �ݺ�ȼ��㣬
	 * @datailһ��ӽ�1��������Ϊ�����Ϻõģ�����3��������Ϊ���õ�������
	 */
	void updateAspectRatio();
	ftype calcAspectRatio(const FaceHandle& fh);

	 /*
	  * @brief ƫб�Ƚ�
	  * @dateil һ�㽨���С��0.7
	  */
	void updateSkewnessEquiangle();
	ftype calcSkewnessEquiarea(const FaceHandle& fh,int cn=3);
	/*
	 * @brief �����
	 * @dateil һ�������ھ����������
	 */
	void updateSkewnessEquiarea();
	ftype calcSkewnessEquiangle(const FaceHandle& fh);
	/*
	 * @brief ���������
	 * @detail
	 */
	void updateAreaRatio();
	ftype calcAreaRatio(const EdgeHandle& eh);
	/*
	 * @brief ����ͶƱ���б���������
	 */
	void updateSmoothVertex(int n=1);
	Point3f calcSmoothVertex(const VertexHandle& vh,int n = 1);
public:
	/*
	 * @brief ������
	 * @detail ���֮ǰ�������deleteIrrVerts(),��������б߽�
	 */
	void checkGenius();
public:
	//ɾ�������㡢�����бߡ�С���
	void deleteIrrVerts(ftype esmallcomp);
	//ɾ�������㡢�����бߡ�����������
	void deleteIrrVerts();
	//ɾ������ı߽���ٷ�
	bool smoothBoundaryDelVertex();
	//����ȥ��
	void repeatMesh(ftype _eps = FLT_EPSILON);
	//����������������ȥ��
	void repeatMeshTextCoord(ftype _eps = FLT_EPSILON);
public:
	//ȥ���ֲ����ݺ��
	void nLocalRemoveAspectRatio(ftype fc = 3);


public://���ʼ��㹤��
	static ftype calcCotanWeight(TriMesh& mesh, const EdgeHandle& eh);
	static ftype calcVoronoiArea(TriMesh& mesh, const VertexHandle& eh);
	static ftype calcVoronoiAreaBarycentric(TriMesh& mesh, const VertexHandle& eh);
	static Point3f calcLaplace(TriMesh& mesh, const VertexHandle& vh);
	static ftype calcAngleSum(TriMesh& mesh, const VertexHandle& v);
	static ftype calcVertexCurvature(TriMesh& mesh, const VertexHandle& v,
		CurvatureType ct= CurvatureType::idGauss);
	static void calcVertexCurvature(TriMesh& mesh, const VertexHandle& v,
		ftype& minvc, ftype& maxvc, ftype& gaussvc, ftype& meanvc);
public: //ͳһ����ply��ʽ�ļ�
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
	static OpenMesh::FPropHandleT<ftype> MinAngle; //��С��
	static OpenMesh::FPropHandleT<ftype> MaxAngle; //����
	static OpenMesh::FPropHandleT<ftype> AspectRatio; //�ݺ��
	static OpenMesh::FPropHandleT<ftype> SkewnessEquiangle; //ƫб�Ƚ�
	static OpenMesh::VPropHandleT<Point3f> SmoothVertex;//������
	static OpenMesh::FPropHandleT<ftype> SkewnessEquiarea; //ƫб��
	static OpenMesh::EPropHandleT<ftype> AreaRatio; //ƫб��
	static OpenMesh::VPropHandleT<bool> Genius; //����
};