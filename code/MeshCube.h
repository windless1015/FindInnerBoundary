#pragma once

#ifndef MESH_CUBE_H
#define MESH_CUBE_H


#include <vector>
#include <map>
#include <iostream>
#include <OpenMeshWarpper.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <vector>
#include <iostream>

//����ƽ���ı��������ռ�����������
//���ó��޲�ֵ�㷨
struct  MeshCube
{
	using Vertex3D = OpenMesh::Vec3f;

	static OpenMesh::VPropHandleT<Vertex3D> u_param;
	static OpenMesh::VPropHandleT<Vertex3D> v_param;
	
	void operator()(
		/*
		*plane:����Fsce[0~umax,(vmax~1)*umax-vmax*umax,(0~vmax)*umax����0~vmax��*umax+umax-1]
		*block: ����js��je��ks��ke��ls��le��ֵע�������˳���ǣ�jkl��
		*/
		double* ptsXX, //x-jkl�������ɢ����֮��(uv�������ɢ����֮��)
		double* ptsYY, //y-jkl�������ɢ����֮��(uv�������ɢ����֮��)
		double* ptsZZ, //z-jkl�������ɢ����֮��(uv�������ɢ����֮��) 
		int* start,//[1,1,1]��https://people.sc.fsu.edu/~jburkardt/cpp_src/tiler_3d/tiler_3d.html��
		int* iend,//face:[umax,vmax,1],Block:[jmax,kmax,lmax]
		unsigned dim1, //face:umax  block:jmax
		unsigned dim2, //face:vmax  block:kmax
		unsigned dim3, //face:1     block:lmax
		unsigned dimax);//face:max(umax,vmax),block:max(jmax,kmax,lmax)


    //! ����Coons���棨�������ڲ�������
	/*!
	  \param  udim u����ά��
	  \param  vdim v����ά��
	  \param  us,ue u����߽����
	  \param  vs,ve v����߽����
	*/
	void MeshCoonsSurface(
		int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve,
		Vertex3D * dompts);


	//! �����ı�������
	/*!
	  \param  udim u����ά��
	  \param  vdim v����ά��
	  \param  us,ue u����߽����
	  \param  vs,ve v����߽����
	  \return
	*/
	void MeshPolyMesh(Polygon_mesh& msh,
		int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve, bool isappran = true);

	//! ��������������
	/*!
	  \param  udim u����ά��
	  \param  vdim v����ά��
	  \param  us,ue u����߽����
	  \param  vs,ve v����߽����
	  \return
	*/
	void MeshTriangMesh(Triangle_mesh& msh,
		int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve,bool isappran=true,bool dir=true);

	void MeshTriangleParam(
		Triangle_mesh& mesh, int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve,
		bool isappran = true, bool dir = true,
		float umin=0,float umax=1,
		float vmin=0,float vmax=1,
		int umaxp = 11, int vmaxp = 11);


	void MeshTriangleMeshTeeth(Triangle_mesh& msh,
		int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve, bool dir = true,
		float tus=0.0f,float tue=1.0f,float tvs=0.0f,float tve=1.0f);


	void GetUVValueCir(int uvdim, double r, double z,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve);
};



#endif