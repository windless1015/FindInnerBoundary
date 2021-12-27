#pragma once
#include <OpenMeshWarpper.h>
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/SparseLU>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <time.h>
#include <iostream>
#include <string>

/*
 * ����������ԭ��
 * 1.������ˮ�ܵ����ף��������ά����������û���Խ�����
 * 2.�ܹ�֧������׶�(������״�ͳߴ�,�������κͷ�����)
 * 3.������ʱ���ܸı�ԭʼ����
 * 4.�����ܹ���Ӧ��ͬ�ĳ��ϣ����ʺͷ����ʻָ���
 * 5.����Ҫ������Ӧԭʼ����ߴ�
 * 6.������Ҫ����߽�Լ������
 */

namespace MeshFix
{
	class HoleFiller
	{
	public:
		using Vertex3D = OpenMesh::Vec3f;
		using Point3f = Vertex3D;
		using ftype = float;
		using vvectorf = std::vector<std::vector<ftype>>;
		using vvectori = std::vector<std::vector<int>>;
	public:
		/*
		 * @brief ���ڿ׶�����Ԥ����
		 * @detail ���ڱ߽���С��nb�Ŀ׶�ֱ�ӽ��в���
		 * @param[in|out] mesh,��ҪԤ���������
		 * @param[in] nb Ԥ�����������߽���
		 * @date 2020.3.10
		 */
		static void pretreatmentHole(Triangle_mesh& mesh, int nb = 10);


	public:
		/*
		 * @brief c0�������񲹶�
		 * @detail �����ó�������������,ǰn�����Ǳ߽��
		 * @param [in] mesh �ṩ�����߽�����,������������
		 * @param [in] hole ����׶����ṩ������Ϊ�˽���
		 * @param [out] hmesh ������������
		 */
		static bool fillHoleLaplaceC0(Triangle_mesh& mesh,
			const std::vector<int>& hole, Triangle_mesh& hmesh);

		/*
		 * @brief c0�������񲹶�
		 * @detail �����ó�������������,ǰn�����Ǳ߽��
		 * @param [in] mesh �ṩ�����߽�����,������������
		 * @param [in] hole ����׶����ṩ������Ϊ�˽���
		 * @param [out] hmesh ������������
		 */
		static bool fillHoleLaplaceC0(Triangle_mesh& mesh,
			const std::vector<int>& hole);

		/*
		 * @brief c0�������񲹶�
		 * @detail �����ó�������������,ǰn�����Ǳ߽��
		 * @param [in] mesh �ṩ�����߽�����,������������
		 * @param [in] hole ����׶����ṩ������Ϊ�˽���
		 * @param [out] hmesh ������������
		 */
		static bool fillHoleLaplaceC0(Triangle_mesh& mesh,
			const std::vector<Vertex3D>& hole);
		/*
		 * @brief c1�������񲹶�
		 * @detail �����ó�������������,ǰn�����Ǳ߽��
		 * @param [in] mesh �ṩ�����߽�����,������������
		 * @param [in] hole ����׶����Լ����������
		 */
		static bool fillHoleLaplaceC1(Triangle_mesh& mesh,
			const std::vector<int>& hole);

		/*
		 * @brief ��������c1��������
		 * @detail �����ó�������������,ǰn�����Ǳ߽��
		 * @param [in] mesh �ṩ�����߽�����,������������
		 * @param [in] hole ����׶����ṩ������Ϊ�˽���
		 * @param [out] hmesh ������������
		 */
		static bool fillHoleLaplaceKFC1(Triangle_mesh& mesh,
			const std::vector<int>& hole, Triangle_mesh& hmesh);

	public:
		//��ȡ���߽�
		static void extraBoundary(Triangle_mesh* mesh, std::vector<Point3f>& boundary, std::vector<int>& ids);
		//ȥ�������α߽���ȡ�㷨(��Է����б߽绷���н���)
		static void extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& aids);
		//��ȡ���б߽磬aids�����б߽��id��ids�����߽磬ik�����߽���Ӧ��id
		static void extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& aids, std::vector<int>& ids, int& ik);
	public:
		template<class T>
		static void debug(std::string filaname, const std::vector<T>& pts) {
			std::ofstream fcout(filaname);
			for (auto ip : pts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << std::endl;
			}
			fcout.close();
		}
	public://�������
		/*
		 * @brief �������(�����ʷ��)
		 * @datail ���ε����������
		 * @param [in] loops ����ϵ�Ȧ
		 * @param [out] tris ��������
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinArea(
			const std::vector<Vertex3D>& loops,
			std::vector<int>& tris);
		/*
		 * @brief �������(��С����Ǽ�Ȩ��Ͽ׶�)
		 * @datail ���ε����������
		 * @param [in] loops ����ϵ�Ȧ
		 * @param [out] tris ��Ϻ����������
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinDihedral(
			const std::vector<Vertex3D>& loops,
			std::vector<int>& tris);
#ifdef _USE_MATH_DEFINES
		/*
		 * @brief �������(��С�������С����Ǽ�Ȩ��Ͽ׶�)
		 * @datail ���ε����������
		 * @param [in] loops ����ϵ�Ȧ
		 * @param [out] tris ��Ϻ����������
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinDihedralArea(
			Triangle_mesh& mesh, const std::vector<OpenMesh::VertexHandle>& loops,
			std::vector<int>& tris);
#endif
		/*
		 * @brief �������(�����ʷ������)
		 * @datail ���ε����������
		 * @param [in] loops ����ϵ�Ȧ
		 * @param [in] normals 3d�����ϵ�Ȧ
		 * @param [out] tris ��Ϻ����������
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinNormal(
			const std::vector<Vertex3D>& loops,
			const std::vector<Vertex3D>& normals,
			std::vector<int>& tris);

	public:
		/*
		 * @brief ���߽�remesh
		 * @param [in|out] mesh
		 * @param [in] avelen remesh��Ŀ��ߴ�
		 * @param [in] angle remesh �ĽǶ�
		 * @author lyc
		 * @date 2020.1.3
		 */
		static bool remeshTriangle(Triangle_mesh& mesh, ftype avelen, ftype angle = 2.0f);
		/*
		 * @brief �������ѡ����
		 * @param [in] mesh ��ѡ������
		 * @param [in] that ��ѡ�ĳ�ʼ�߽�
		 * @param [in] n nȦ����
		 * @param [out] spts ѡ��ĵ��id
		 * @author lyc
		 * @date 2020.1.3
		 */
		static bool selectVerts(Triangle_mesh& mesh, const std::vector<int>& that, 
			int n, std::vector<int>&spts, bool issave = false);
		static bool selectVerts(Triangle_mesh& mesh, const VertexHandle& that, 
			int n, std::vector<VertexHandle>&spts, bool issave = false);
		static bool selectVerts(Triangle_mesh& mesh, const std::vector<VertexHandle>& that, 
			int n, std::vector<VertexHandle>&spts, bool issave = false);
	public:
	protected:
		static Point3f rbf_Intelpoation(std::vector<Point3f>& pts,Point3f& udir,Point3f& vdir,Point3f& wdir, Point3f& p);
		static bool nlocalFairMesh(Triangle_mesh& mesh, std::vector<VertexHandle>& vhs, int k);
	};
}