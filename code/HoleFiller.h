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
 * 补洞的六大原则：
 * 1.补洞能水密的填充孔，且满足二维流形条件，没有自交现象
 * 2.能够支持任意孔洞(任意形状和尺寸,包括流形和非流形)
 * 3.补洞的时候不能改变原始数据
 * 4.补洞能够适应不同的场合（曲率和非曲率恢复）
 * 5.补洞要尽量适应原始网格尺寸
 * 6.补洞需要满足边界约束条件
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
		 * @brief 对于孔洞进行预处理
		 * @detail 对于边界数小于nb的孔洞直接进行补洞
		 * @param[in|out] mesh,需要预处理的网格
		 * @param[in] nb 预处理补洞的最大边界数
		 * @date 2020.3.10
		 */
		static void pretreatmentHole(Triangle_mesh& mesh, int nb = 10);


	public:
		/*
		 * @brief c0连续网格补洞
		 * @detail 单独拿出来方便做并行,前n个点是边界点
		 * @param [in] mesh 提供连续边界条件,待补洞的网格
		 * @param [in] hole 网格孔洞，提供索引是为了解耦
		 * @param [out] hmesh 补洞的新网格
		 */
		static bool fillHoleLaplaceC0(Triangle_mesh& mesh,
			const std::vector<int>& hole, Triangle_mesh& hmesh);

		/*
		 * @brief c0连续网格补洞
		 * @detail 单独拿出来方便做并行,前n个点是边界点
		 * @param [in] mesh 提供连续边界条件,待补洞的网格
		 * @param [in] hole 网格孔洞，提供索引是为了解耦
		 * @param [out] hmesh 补洞的新网格
		 */
		static bool fillHoleLaplaceC0(Triangle_mesh& mesh,
			const std::vector<int>& hole);

		/*
		 * @brief c0连续网格补洞
		 * @detail 单独拿出来方便做并行,前n个点是边界点
		 * @param [in] mesh 提供连续边界条件,待补洞的网格
		 * @param [in] hole 网格孔洞，提供索引是为了解耦
		 * @param [out] hmesh 补洞的新网格
		 */
		static bool fillHoleLaplaceC0(Triangle_mesh& mesh,
			const std::vector<Vertex3D>& hole);
		/*
		 * @brief c1连续网格补洞
		 * @detail 单独拿出来方便做并行,前n个点是边界点
		 * @param [in] mesh 提供连续边界条件,待补洞的网格
		 * @param [in] hole 网格孔洞，自己计算出单环
		 */
		static bool fillHoleLaplaceC1(Triangle_mesh& mesh,
			const std::vector<int>& hole);

		/*
		 * @brief 保持特征c1连续补洞
		 * @detail 单独拿出来方便做并行,前n个点是边界点
		 * @param [in] mesh 提供连续边界条件,待补洞的网格
		 * @param [in] hole 网格孔洞，提供索引是为了解耦
		 * @param [out] hmesh 补洞的新网格
		 */
		static bool fillHoleLaplaceKFC1(Triangle_mesh& mesh,
			const std::vector<int>& hole, Triangle_mesh& hmesh);

	public:
		//提取最大边界
		static void extraBoundary(Triangle_mesh* mesh, std::vector<Point3f>& boundary, std::vector<int>& ids);
		//去除非流形边界提取算法(会对非流行边界环进行解耦)
		static void extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& aids);
		//提取所有边界，aids是所有边界的id，ids是最大边界，ik是最大边界相应的id
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
	public://缝合网格
		/*
		 * @brief 缝合网格(无曲率缝合)
		 * @datail 流形单环，不打结
		 * @param [in] loops 待缝合的圈
		 * @param [out] tris 三角网格
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinArea(
			const std::vector<Vertex3D>& loops,
			std::vector<int>& tris);
		/*
		 * @brief 缝合网格(最小二面角加权缝合孔洞)
		 * @datail 流形单环，不打结
		 * @param [in] loops 待缝合的圈
		 * @param [out] tris 缝合后的三角网格
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinDihedral(
			const std::vector<Vertex3D>& loops,
			std::vector<int>& tris);
#ifdef _USE_MATH_DEFINES
		/*
		 * @brief 缝合网格(最小面积与最小二面角加权缝合孔洞)
		 * @datail 流形单环，不打结
		 * @param [in] loops 待缝合的圈
		 * @param [out] tris 缝合后的三角网格
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinDihedralArea(
			Triangle_mesh& mesh, const std::vector<OpenMesh::VertexHandle>& loops,
			std::vector<int>& tris);
#endif
		/*
		 * @brief 缝合网格(保曲率缝合网格)
		 * @datail 流形单环，不打结
		 * @param [in] loops 待缝合的圈
		 * @param [in] normals 3d网格缝合的圈
		 * @param [out] tris 缝合后的三角网格
		 * @date 2019.12.31
		 * @author lyc
		 */
		static bool sewUpLoopMinNormal(
			const std::vector<Vertex3D>& loops,
			const std::vector<Vertex3D>& normals,
			std::vector<int>& tris);

	public:
		/*
		 * @brief 保边界remesh
		 * @param [in|out] mesh
		 * @param [in] avelen remesh的目标尺寸
		 * @param [in] angle remesh 的角度
		 * @author lyc
		 * @date 2020.1.3
		 */
		static bool remeshTriangle(Triangle_mesh& mesh, ftype avelen, ftype angle = 2.0f);
		/*
		 * @brief 网格快速选择器
		 * @param [in] mesh 待选的网格
		 * @param [in] that 待选的初始边界
		 * @param [in] n n圈领域
		 * @param [out] spts 选择的点的id
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