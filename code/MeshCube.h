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

//生成平面四边形网格或空间六面体网格
//采用超限插值算法
struct  MeshCube
{
	using Vertex3D = OpenMesh::Vec3f;

	static OpenMesh::VPropHandleT<Vertex3D> u_param;
	static OpenMesh::VPropHandleT<Vertex3D> v_param;
	
	void operator()(
		/*
		*plane:给出Fsce[0~umax,(vmax~1)*umax-vmax*umax,(0~vmax)*umax，（0~vmax）*umax+umax-1]
		*block: 给出js，je，ks，ke，ls，le的值注意遍历的顺序是（jkl）
		*/
		double* ptsXX, //x-jkl方向的离散点数之积(uv方向的离散点数之积)
		double* ptsYY, //y-jkl方向的离散点数之积(uv方向的离散点数之积)
		double* ptsZZ, //z-jkl方向的离散点数之积(uv方向的离散点数之积) 
		int* start,//[1,1,1]（https://people.sc.fsu.edu/~jburkardt/cpp_src/tiler_3d/tiler_3d.html）
		int* iend,//face:[umax,vmax,1],Block:[jmax,kmax,lmax]
		unsigned dim1, //face:umax  block:jmax
		unsigned dim2, //face:vmax  block:kmax
		unsigned dim3, //face:1     block:lmax
		unsigned dimax);//face:max(umax,vmax),block:max(jmax,kmax,lmax)


    //! 生成Coons曲面（可以用于参数化）
	/*!
	  \param  udim u方向维数
	  \param  vdim v方向维数
	  \param  us,ue u方向边界点数
	  \param  vs,ve v方向边界点数
	*/
	void MeshCoonsSurface(
		int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve,
		Vertex3D * dompts);


	//! 生成四边形网格
	/*!
	  \param  udim u方向维数
	  \param  vdim v方向维数
	  \param  us,ue u方向边界点数
	  \param  vs,ve v方向边界点数
	  \return
	*/
	void MeshPolyMesh(Polygon_mesh& msh,
		int udim, int vdim,
		std::vector<Vertex3D> &us,
		std::vector<Vertex3D> &ue,
		std::vector<Vertex3D> &vs,
		std::vector<Vertex3D> &ve, bool isappran = true);

	//! 生成三角形网格
	/*!
	  \param  udim u方向维数
	  \param  vdim v方向维数
	  \param  us,ue u方向边界点数
	  \param  vs,ve v方向边界点数
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