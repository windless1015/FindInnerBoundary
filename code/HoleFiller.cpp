#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include "RBF_interpolation.hpp"
#include "HoleFiller.h"
#include "MeshDoctor.h"
#include <iostream>
#include <map>
#include <stack>
using namespace std;
namespace MeshFix
{
	static OpenMesh::VPropHandleT<int> DEGREE;
	static OpenMesh::VPropHandleT<bool> isFix;
	using SparseMatrix = Eigen::SparseMatrix<HoleFiller::ftype>;
	using Triplet = Eigen::Triplet<HoleFiller::ftype>;
	struct weight
	{
		using ftype = HoleFiller::ftype;
		weight(ftype _angle = FLT_MAX, ftype _area = FLT_MAX)
			: s_angle(_angle), s_area(_area) {}
		weight operator+(const weight& _rhs) const {
			return weight(std::max(s_angle, _rhs.s_angle), s_area + _rhs.s_area);
		}
		bool operator<(const weight& _rhs) const {
			return (s_angle < _rhs.s_angle ||
				(s_angle == _rhs.s_angle && s_area < _rhs.s_area));
		}
		ftype s_angle;
		ftype s_area;
	};

	/*
	* @brief 提供数值修正与检测算法
	* @date 2020.3.16
	*/
	template <class Real>
	struct utils {
		static Real realThreshold() { return Real(FLT_EPSILON); }
		static Real realPrecision() { return Real(FLT_EPSILON); }
		static Real pi() { return Real(3.14159265358979323846264338327950288); }
		static Real tau() { return Real(2) * pi<Real>(); }
		//正向修正
		static Real fixPos(Real rl) {
			if (rl >= 0 && rl < realPrecision())
				return realPrecision();
			else if (rl<0 && rl>-realPrecision())
				return -realPrecision();
			else
				return rl;
		}
		//负向修正
		static Real fixNag(Real rl) {
			if (rl > 0 && rl < realPrecision())
				return realPrecision();
			else if (rl <= 0 && rl > -realPrecision())
				return -realPrecision();
			else
				return rl;
		}
		//检查数值有效性
		static bool checkEval(Real rl) {
			return !std::isnan<Real>(rl) && !std::isinf<Real>(rl);
		}
		//检查数值的零值
		static bool isZero(Real rl) {
			return rl >= -realPrecision() && rl <= realPrecision();
		}
		static bool isEqual(Real r1, Real r2) {
			return isZero(r1 - r2);
		}
	};
	using util = utils<HoleFiller::ftype>;

	//预处理网格
	void HoleFiller::pretreatmentHole(Triangle_mesh& mesh, int nb)
	{
		//1.删除孤立点
		int ncount = 0;
		MeshDoctor mdc(mesh);
		mdc.deleteIrrVerts();
		mdc.repeatMesh();
		stack<int> nboundary;
		vector<vector<int>> ids;
		extraBoundary(&mesh, ids);
		nboundary.push((int)ids.size());
		//小于nb的边界自动填充
		while (true && ncount <= 10)
		{
			for (auto& ip : ids)
				if (ip.size() <= nb)
					fillHoleLaplaceC0(mesh, ip);
			extraBoundary(&mesh, ids);
			bool isok = true;
			for (auto& icp : ids) {
				if (icp.size() < nb) {
					isok = false;
				}
			}
			if (!isok) {
				if (nboundary.top() == ids.size())
					++ncount;
				nboundary.push((int)ids.size());
			}
			if (isok)
				break;
		}
	}

	/*
	* @brief c0连续网格补洞
	* @detail 单独拿出来方便做并行,前n个点是边界点
	* @param [in] mesh 提供连续边界条件
	* @param [in] hole 网格孔洞，提供索引是为了解耦
	* @param [out] ntriangle 补洞的点
	* @param [out] tris 三角形下表
	*/
	bool HoleFiller::fillHoleLaplaceC0(Triangle_mesh& mesh,
		const std::vector<int>& hole, Triangle_mesh& hmesh)
	{
		hmesh.clear();
		vector<VertexHandle> cloops;
		int nhole = static_cast<int>(hole.size());

		auto chole = hole;
		sort(chole.begin(), chole.end());
		unique(chole.begin(), chole.end());
		if (hole.size() != chole.size())
			return false;
		cloops.reserve(nhole);
		auto points = mesh.points();
		ftype clen = static_cast<ftype>(0.0);
		for (int i = 0; i < nhole; ++i) {
			auto p = points[hole[i]];
			auto cp = points[hole[(i + 1) % nhole]];
			clen += static_cast<ftype>((p - cp).length());
		}
		ftype pclen = static_cast<ftype>(0.0);
		for (int i = 0; i < nhole; ++i) {
			cloops.emplace_back(VertexHandle(hole[i]));
		}
		vector<int> ctris;
		sewUpLoopMinDihedralArea(mesh, cloops, ctris);
		hmesh.reserve(cloops.size(), ctris.size() * 3, ctris.size() / 3);
		vector<VertexHandle> vhs;
		for (auto ip : cloops)
			vhs.push_back(hmesh.add_vertex(mesh.point(ip)));
		int ntris = static_cast<int>(ctris.size());
		for (int i = 0; i < ntris; i += 3) {
			hmesh.add_face(vhs[ctris[i]], vhs[ctris[i + 1]], vhs[ctris[i + 2]]);
		}
		remeshTriangle(hmesh, static_cast<ftype>(clen / nhole));
		if (!hmesh.has_vertex_status())
			hmesh.request_vertex_status();
		for (int i = 0; i < nhole; ++i) {
			hmesh.status(VertexHandle(i)).set_locked(true);
		}
		for (int i = nhole; i < hmesh.n_vertices(); ++i) {
			hmesh.status(VertexHandle(i)).set_locked(false);
		}
		//最小二乘原始网格
		OpenMesh::Smoother::JacobiLaplaceSmootherT< Triangle_mesh > smoother(hmesh);
		smoother.initialize(OpenMesh::Smoother::SmootherT< Triangle_mesh >::
			Tangential_and_Normal, OpenMesh::Smoother::SmootherT< Triangle_mesh >::C1);
		smoother.smooth(20);
		hmesh.release_vertex_status();
		return true;
	}

	/*
	 * @brief c0连续网格补洞
	 * @detail 单独拿出来方便做并行,前n个点是边界点
	 * @param [in] mesh 提供连续边界条件
	 * @param [in] hole 网格孔洞，提供索引是为了解耦
	 * @param [out] ntriangle 补洞的点
	 * @param [out] tris 三角形下表
	 */
	bool HoleFiller::fillHoleLaplaceC0(Triangle_mesh& mesh,
		const std::vector<int>& hole)
	{
		Triangle_mesh hmesh;
		hmesh.clear();
		vector<Vertex3D> cloops;
		int nhole = static_cast<int>(hole.size());
		if (hole.size() < 3)
			return false;
		auto chole = hole;
		sort(chole.begin(), chole.end());
		unique(chole.begin(), chole.end());
		if (hole.size() != chole.size())
			return false;
		cloops.reserve(nhole);
		auto points = mesh.points();
		ftype clen = static_cast<ftype>(0.0);
		for (int i = 0; i < nhole; ++i) {
			auto p = points[hole[i]];
			auto cp = points[hole[(i + 1) % nhole]];
			clen += static_cast<ftype>((p - cp).length());
		}
		ftype pclen = static_cast<ftype>(0.0);
		for (int i = 0; i < nhole; ++i) {
			auto p = points[hole[i]];
			cloops.emplace_back(p);
		}
		vector<int> ctris;
		vector<VertexHandle> handles;
		for (auto ip : hole) {
			handles.push_back(VertexHandle(ip));
		}
		sewUpLoopMinDihedralArea(mesh, handles, ctris);
		hmesh.reserve(nhole, ctris.size() * 3, ctris.size() / 3);
		for (auto ip : cloops)
			hmesh.add_vertex(ip);
		int ntris = static_cast<int>(ctris.size());
		for (int i = 0; i < ntris; i += 3) {
			hmesh.add_face(VertexHandle(ctris[i]), VertexHandle(ctris[i + 1]), VertexHandle(ctris[i + 2]));
		}
		vector<VertexHandle> chv;
		if (remeshTriangle(hmesh, static_cast<ftype>(clen / nhole)))
		{
			for (int i = 0; i < nhole; ++i) {
				chv.push_back(VertexHandle(hole[i]));
			}
			
			for (int i = nhole; i < hmesh.n_vertices(); ++i) {
				chv.push_back(mesh.add_vertex(hmesh.point(VertexHandle(i))));
			}
			for (auto ff = hmesh.faces_begin(); ff != hmesh.faces_end(); ++ff)
			{
				auto fv = hmesh.fv_begin(*ff);
				auto p1 = fv->idx(); ++fv;
				auto p2 = fv->idx(); ++fv;
				auto p3 = fv->idx();
				if (!mesh.add_face(chv[p1], chv[p2], chv[p3]).is_valid())
					mesh.add_face(chv[p1], chv[p3], chv[p2]);
			}
			vector<VertexHandle> smoothverts;
			selectVerts(mesh, chv, 3, smoothverts, true);
			if (!mesh.has_vertex_status())
				mesh.request_vertex_status();
			vector<bool> ismooth(mesh.n_vertices(), true);
			for (auto& iv : smoothverts)
				ismooth[iv.idx()] = false;
			for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv)
				mesh.status(*iv).set_locked(ismooth[iv->idx()]);
			//最小二乘原始网格
			OpenMesh::Smoother::JacobiLaplaceSmootherT< Triangle_mesh > smoother(mesh);
			smoother.initialize(OpenMesh::Smoother::SmootherT< Triangle_mesh >::
				Tangential_and_Normal, OpenMesh::Smoother::SmootherT< Triangle_mesh >::C1);
			smoother.smooth(20);
			mesh.release_vertex_status();
		}
		return true;
	}


	/*
  * @brief c1连续网格补洞
  * @detail 单独拿出来方便做并行,前n个点是边界点
  * @param [in] mesh 提供连续边界条件
  * @param [in] hole 网格孔洞，提供索引是为了解耦
  * @param [out] ntriangle 补洞的点
  * @param [out] tris 三角形下表
  */
	bool HoleFiller::fillHoleLaplaceC0(Triangle_mesh& mesh,
		const std::vector<Vertex3D>& hole)
	{
		Triangle_mesh hmesh;
		hmesh.clear();
		vector<Vertex3D> cloops;
		int nhole = static_cast<int>(hole.size());
		if (hole.size() < 3)
			return false;
		auto chole = hole;
		sort(chole.begin(), chole.end());
		unique(chole.begin(), chole.end());
		if (hole.size() != chole.size())
			return false;
		cloops.reserve(nhole);
		auto points = mesh.points();
		ftype clen = static_cast<ftype>(0.0);
		for (int i = 0; i < nhole; ++i) {
			auto p = hole[i];
			auto cp = hole[(i + 1) % nhole];
			clen += static_cast<ftype>((p - cp).length());
		}
		ftype pclen = static_cast<ftype>(0.0);
		for (int i = 0; i < nhole; ++i) {
			auto p = hole[i];
			cloops.emplace_back(p);
		}
		vector<int> ctris;
		sewUpLoopMinArea(cloops, ctris);
		hmesh.reserve(cloops.size(), ctris.size() * 3, ctris.size() / 3);
		for (auto ip : cloops)
			hmesh.add_vertex(ip);
		int ntris = static_cast<int>(ctris.size());
		for (int i = 0; i < ntris; i += 3) {
			hmesh.add_face(VertexHandle(ctris[i]), VertexHandle(ctris[i + 1]), VertexHandle(ctris[i + 2]));
		}
		
		vector<VertexHandle> chv;
		if (remeshTriangle(hmesh, static_cast<ftype>(clen / nhole)))
		{
			if (!hmesh.has_vertex_status())
				hmesh.request_vertex_status();
			for (auto iv = hmesh.vertices_begin(); iv != hmesh.vertices_end(); ++iv) {
				if (hmesh.is_boundary(*iv))
					hmesh.status(*iv).set_locked(true);
				else {
					hmesh.status(*iv).set_locked(false);
				}
			}
			//最小二乘原始网格
			OpenMesh::Smoother::JacobiLaplaceSmootherT< Triangle_mesh > smoother(hmesh);
			smoother.initialize(OpenMesh::Smoother::SmootherT< Triangle_mesh >::
				Tangential_and_Normal, OpenMesh::Smoother::SmootherT< Triangle_mesh >::C1);
			smoother.smooth(20);
			mesh.release_vertex_status();
			for (int i = 0; i < hmesh.n_vertices(); ++i) {
				chv.push_back(mesh.add_vertex(hmesh.point(VertexHandle(i))));
			}
			for (auto ff = hmesh.faces_begin(); ff != hmesh.faces_end(); ++ff)
			{
				auto fv = hmesh.fv_begin(*ff);
				auto p1 = fv->idx(); ++fv;
				auto p2 = fv->idx(); ++fv;
				auto p3 = fv->idx();
				if (!mesh.add_face(chv[p1], chv[p2], chv[p3]).is_valid())
					mesh.add_face(chv[p1], chv[p3], chv[p2]);
			}
		}
		hmesh.release_vertex_status();
		return true;
	}

	void DebugPts(std::string filename, std::vector<OpenMesh::Vec3f> &boundary)
	{
		ofstream fcout(filename);
		for (auto ip : boundary)
		{
			fcout << ip[0] << " " << ip[1] << " " << ip[2] << endl;
		}
		fcout.close();
	}

	/*
	 * @brief c1连续网格补洞
	 * @detail 单独拿出来方便做并行,前n个点是边界点
	 * @param [in] mesh 提供连续边界条件,待补洞的网格
	 * @param [in] hole 网格孔洞，提供索引是为了解耦
	 * @param [out] hmesh 补洞的新网格
	 */
	bool HoleFiller::fillHoleLaplaceC1(Triangle_mesh& mesh, const std::vector<int>& hole)
	{
		Triangle_mesh hmesh;
		fillHoleLaplaceC0(mesh, hole, hmesh);
		hmesh.request_vertex_status();
		int nhole = static_cast<int>(hole.size());
		int nverts = static_cast<int>(hmesh.n_vertices());
		vector<VertexHandle> handles;
		handles.reserve(nverts);
		for (auto ih : hole) {
			handles.push_back(VertexHandle(ih));
		}
		for (int i = nhole; i < nverts; ++i) {
			auto vh = mesh.add_vertex(hmesh.point(VertexHandle(i)));
			handles.emplace_back(vh);
		}
		if (!mesh.has_vertex_status())
			mesh.request_vertex_status();

		for (auto f = hmesh.faces_begin(); f != hmesh.faces_end(); ++f) {
			auto fv = hmesh.fv_begin(*f);
			auto o = fv->idx(); ++fv;
			auto t = fv->idx(); ++fv;
			auto h = fv->idx();
			mesh.add_face(handles[o], handles[t], handles[h]);
		}
		//最小二乘化网格
		vector<bool> isok(mesh.n_vertices(), true);
		vector<VertexHandle> smoothverts;
		selectVerts(mesh, handles, 4, smoothverts, true);
		for (auto& iv : smoothverts)
			isok[iv.idx()] = false;
		for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
			mesh.status(*iv).set_locked(isok[iv->idx()]);
		}
		OpenMesh::Smoother::JacobiLaplaceSmootherT< Triangle_mesh > smoother(mesh);
		smoother.initialize(OpenMesh::Smoother::SmootherT< Triangle_mesh >::
			Tangential_and_Normal, OpenMesh::Smoother::SmootherT< Triangle_mesh >::C1);
		smoother.smooth(20);
		mesh.release_vertex_status();
		if (smoothverts.size() > handles.size()) {
			//对handle区域进行曲率恢复
			selectVerts(mesh, handles, 2, smoothverts, true);
			if (nlocalFairMesh(mesh, smoothverts, 2))
				cout << "hole fill is ok!" << endl;
			else {
				if(nlocalFairMesh(mesh, smoothverts, 1))
					cout << "hole fill is ok!" << endl;
				else {
					cout << "hole fill is ok,but hole curvature recovery is fail!" << endl;
				}
			}
		}
		return true;
	}

	HoleFiller::Point3f HoleFiller::rbf_Intelpoation(std::vector<Point3f>& pts,
		Point3f& udir, Point3f& vdir, Point3f& wdir, Point3f& p)
	{
		//参数化pts
		assert(pts.size() > 0 && "the size of pts must be max 0!");
		using Rbf = FussenAlgo::RBFCore<ftype>;
		Point3f ccp = pts[0];
		udir.normalize();
		vdir.normalize();
		wdir.normalize();
		Rbf::Matrix x, y;
		int npt = static_cast<int>(pts.size());
		x.resize(npt, 2);
		y.resize(npt, 1);
		for (int i = 0; i < npt; ++i) {
			Point3f cp = pts[i] - ccp;
			x(i, 0) = cp | udir;
			y(i, 0) = cp | vdir;
			y(i, 0) = cp | wdir;
		}
		Rbf::Matrix rbfCoff = Rbf::rbfcreate(x, y, Rbf::rbfphi_cubic, 1.0f, 0.0f);
		Rbf::Matrix xx, yy;
		xx.resize(1, 2);
		Point3f dp = p - ccp;
		xx(0, 0) = dp | udir;
		xx(0, 1) = dp | vdir;
		yy = Rbf::rbfinterp(x, rbfCoff, xx, Rbf::rbfphi_cubic, 1.0f);
		Point3f aimp = ccp + xx(0, 0)*udir + xx(0, 1)*vdir + yy(0, 0)*wdir;
		return aimp;
	}

	/*
	 * @brief 保持特征c1连续补洞
	 * @detail 单独拿出来方便做并行,前n个点是边界点
	 * @param [in] mesh 提供连续边界条件,待补洞的网格
	 * @param [in] hole 网格孔洞，提供索引是为了解耦
	 * @param [out] hmesh 补洞的新网格
	 */
	bool HoleFiller::fillHoleLaplaceKFC1(Triangle_mesh& mesh,
		const std::vector<int>& hole, Triangle_mesh& hmesh)
	{



		return true;
	}

	/*
	 * @brief 缝合网格
	 * @param [in] loops 待缝合的圈
	 * @param [out] tris 三角网格
	 * @date 2019.12.31
	 * @author lyc
	 */
	bool HoleFiller::sewUpLoopMinArea(const std::vector<Vertex3D>& loops, std::vector<int>& tris)
	{
		int size = static_cast<int>(loops.size());
		if (size < 3)
			return false;
		auto area = [&loops](int i, int j, int k) {
			Vertex3D a = loops[i];
			Vertex3D b = loops[j];
			Vertex3D c = loops[k];
			Vertex3D n((b - a) % (c - b));
			return static_cast<ftype>(0.5 * n.norm());
		};
		std::function<void(vector<int>& tris,
			const vector<int>& hole_id,
			HoleFiller::vvectori &minimum_weight_index,
			int begin, int end)>AddHoleToMesh;
		AddHoleToMesh = [&AddHoleToMesh](
			vector<int>& tris,
			const vector<int>& hole_id,
			vvectori &minimum_weight_index,
			int begin, int end) {
			if (end - begin > 1) {
				int cu = minimum_weight_index[begin][end];
				tris.push_back(hole_id[begin]);
				tris.push_back(hole_id[cu]);
				tris.push_back(hole_id[end]);
				AddHoleToMesh(tris, hole_id, minimum_weight_index, begin, cu);
				AddHoleToMesh(tris, hole_id, minimum_weight_index, cu, end);
			}
		};

		vvectorf minimum_weight(size, vector<ftype>(size, 0));
		vvectori minimum_weight_index(size, vector<int>(size, -1));
		vector<int> ids;
		ids.reserve(size);
		for (int ic = 0; ic < size; ++ic)
			ids.push_back(ic);
		tris.clear();
		tris.reserve(size * 3);
		for (int j = 2; j < size; ++j) {
			for (int i = 0; i < size - j; ++i) {
				ftype min = FLT_MAX;
				int index = -1;
				int k = i + j;
				for (int m = i + 1; m < k; m++) {
					ftype farea = area(i, m, k);
					ftype val = minimum_weight[i][m] + minimum_weight[m][k] + farea;
					if (val < min) {
						min = val; index = m;
					}
				}
				minimum_weight[i][k] = min;
				minimum_weight_index[i][k] = index;
			}
		}
		AddHoleToMesh(tris, ids, minimum_weight_index, 0, size - 1);
		return true;
	}

#ifdef _USE_MATH_DEFINES
	/*
	* @brief 缝合网格(最小面积与最小二面角加权缝合孔洞)
	* @datail 流形单环，不打结
	* @param [in] loops 待缝合的圈
	* @param [out] tris 缝合后的三角网格
	* @date 2019.12.31
	* @author lyc
	*/
	bool HoleFiller::sewUpLoopMinDihedralArea(
		Triangle_mesh& mesh,
		const std::vector<OpenMesh::VertexHandle>& cloops,
		std::vector<int>& tris)
	{
		tris.clear();
		auto points = mesh.points();
		int nhole = static_cast<int>(cloops.size());
		vector<Vertex3D> loops(nhole, Vertex3D());
		for (int i = 0; i < nhole; ++i)
			loops[i] = points[cloops[i].idx()];
		vector<vector<weight>> cweight(nhole, vector<weight>(nhole, weight()));
		vector<vector<int>> cindex(nhole, vector<int>(nhole, 0));
		for (int i = 0; i < nhole - 1; ++i) {
			cweight[i][i + 1] = weight(0, 0);
			cindex[i][i + 1] = 0;
		}
		std::function<weight(int i, int j, int k)>compute_weight;
		std::function<ftype(int i, int j, int k)> compute_area;
		compute_area = [&loops](int i, int j, int k) {
			Vertex3D a = loops[i];
			Vertex3D b = loops[j];
			Vertex3D c = loops[k];
			Vertex3D n((b - a) % (c - b));
			return static_cast<ftype>(0.5 * n.norm());
		};
		auto opposite_vertex = [&mesh, &points](OpenMesh::VertexHandle vh) {
			auto cvh = mesh.to_vertex_handle
			(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(mesh.halfedge_handle(vh))));
			return points[cvh.idx()];
		};
		auto exists_edge = [&mesh](const VertexHandle& u, const VertexHandle& w) {
			for (auto vohi = mesh.voh_iter(u); vohi.is_valid(); ++vohi)
				if (!mesh.is_boundary(mesh.edge_handle(*vohi)))
					if (mesh.to_vertex_handle(*vohi) == w)
						return true;
			return false;
		};
		auto compute_angle = [](const Vertex3D& u, const Vertex3D& v, const Vertex3D& a, const Vertex3D& b) {
			Vertex3D n0((v - u) % (a - v));
			Vertex3D n1((u - v) % (b - u));
			n0.normalize();
			n1.normalize();
			return static_cast<ftype>(acos(n0 | n1) * 180.0 / M_PI);
		};
		compute_weight =
			[&cloops, &loops, &compute_area, &cindex,
			&exists_edge, &compute_angle, &opposite_vertex](int i, int j, int k) {
			if (exists_edge(cloops[i], cloops[j]) ||
				exists_edge(cloops[j], cloops[k]) ||
				exists_edge(cloops[k], cloops[i]))
				return weight();
			if (cindex[i][j] == -1) return weight();
			if (cindex[j][k] == -1) return weight();

			auto a = loops[i];
			auto b = loops[j];
			auto c = loops[k];
			ftype area = compute_area(i, j, k);
			ftype angle = static_cast<ftype>(0);
			if (i + 1 == j)
				angle = std::max(angle, compute_angle(
					loops[i], loops[j], loops[k], opposite_vertex(cloops[i])));
			else
				angle = std::max(angle, compute_angle(
					loops[i], loops[j], loops[k], loops[cindex[i][j]]));
			if (j + 1 == k)
				angle = std::max(angle, compute_angle(
					loops[j], loops[k], loops[i], opposite_vertex(cloops[j])));
			else
				angle = std::max(angle, compute_angle(
					loops[j], loops[k], loops[i], loops[cindex[j][k]]));
			if (i == 0 && k == (int)loops.size() - 1)
				angle = std::max(angle, compute_angle(loops[k],
					loops[i], loops[j], opposite_vertex(cloops[k])));
			return weight(angle, area);
		};
		for (int j = 2; j < nhole; ++j) {
#pragma omp parallel for shared(j, nhole)
			for (int i = 0; i < nhole - j; ++i) {
				int k = i + j;
				weight wmin = weight();
				int imin = -1;
				for (int m = i + 1; m < k; ++m) {
					weight w = cweight[i][m] + compute_weight(i, m, k) + cweight[m][k];
					if (w < wmin) {
						wmin = w;
						imin = m;
					}
				}
				cweight[i][k] = wmin;
				cindex[i][k] = imin;
			}
		}
		std::function<bool(int i, int j)>fill;
		fill = [&cindex, &tris, &fill](int i, int j) {
			if (i + 1 == j)
				return true;
			tris.push_back(i);
			tris.push_back(cindex[i][j]);
			tris.push_back(j);
			if (i == j || i == cindex[i][j] || j == cindex[i][j])
				return false;
			if (!fill(i, cindex[i][j]) || !fill(cindex[i][j], j))
				return false;
			else
				return true;
		};
		fill(0, nhole - 1);
		return true;
	}
#endif
	/*
	 * @brief 缝合网格(最小面积加权缝合孔洞)
	 * @datail 流形单环，不打结
	 * @param [in] loops 待缝合的圈
	 * @param [out] tris 缝合后的三角网格
	 * @date 2019.12.31
	 * @author lyc
	 */
	bool HoleFiller::sewUpLoopMinDihedral(const std::vector<Vertex3D>& loops, std::vector<int>& tris)
	{
		int nhole = static_cast<int>(loops.size());
		if (nhole < 3)
			return false;
		using vvectorw = vector<vector<weight>>;
		vvectorw cw;
		cw.resize(nhole, vector<weight>(nhole, weight()));
		vvectori cl;
		cl.resize(nhole, vector<int>(nhole, 0));
		for (int i = 0; i < nhole - 1; ++i)
			cw[i][i + 1] = weight(0, 0);
		auto dihedral = [&loops](int i, int j, int k, int l) {
			Vertex3D u = loops[i];
			Vertex3D v = loops[j];
			Vertex3D a = loops[k];
			Vertex3D b = loops[l];
			Vertex3D n0((v - u) % (a - v));
			Vertex3D n1((u - v) % (b - u));
			n0.normalize();
			n1.normalize();
			return static_cast<ftype>(acos(n0 | n1) * 180.0 / M_PI);
		};
		auto area = [&loops](int i, int j, int k) {
			Vertex3D a = loops[i];
			Vertex3D b = loops[j];
			Vertex3D c = loops[k];
			Vertex3D n((b - a) % (c - b));
			return static_cast<ftype>(0.5 * n.norm());
		};
		auto cweight = [&cl, &dihedral, &nhole, &area](int i, int j, int k) {
			if (cl[i][j] == -1)return weight();
			if (cl[j][k] == -1)return weight();
			ftype angle = static_cast<ftype>(0.0);
			if (i + 1 == j)
				angle = std::max(angle, dihedral(i, j, k, (i - 1 + nhole) % nhole));
			else
				angle = std::max(angle, dihedral(i, j, k, cl[i][j]));
			if (j + 1 == k)
				angle = std::max(angle, dihedral(j, k, i, (j - 1 + nhole) % nhole));
			else
				angle = std::max(angle, dihedral(j, k, i, cl[j][k]));
			if (i == 0 && k == nhole - 1)
				angle = std::max(angle, dihedral(k, i, j, (k - 1 + nhole) % nhole));
			return weight(angle, area(i, j, k));
		};
		for (int j = 2; j < nhole; ++j) {
#pragma omp parallel for shared(j, nhole)
			for (int i = 0; i < nhole - j; ++i) {
				weight valmin;
				int   argmin = -1;
				for (int m = i + 1; m < i + j; ++m) {
					weight newval = cw[i][m] + cw[m][i + j] + cweight(i, m, i + j);
					if (newval < valmin) {
						valmin = newval;
						argmin = m;
					}
				}
				cw[i][i + j] = valmin;
				cl[i][i + j] = argmin;
			}
		}
		std::function<bool(int i, int j)>fill;
		fill = [&cl, &tris, &fill](int i, int j) {
			if (i + 1 == j)
				return true;
			tris.push_back(i);
			tris.push_back(cl[i][j]);
			tris.push_back(j);
			if (i == j || i == cl[i][j] || j == cl[i][j])
				return false;
			if (!fill(i, cl[i][j]) || !fill(cl[i][j], j))
				return false;
			else
				return true;
		};
		tris.clear();
		tris.reserve(nhole * 3);
		return fill(0, nhole - 1);
	}


	/*
	* @brief 缝合网格(保曲率缝合网格)
	* @datail 流形单环，不打结
	* @param [in] loops 待缝合的圈
	* @param [out] tris 缝合后的三角网格
	* @date 2019.12.31
	* @author lyc
	*/
	bool HoleFiller::sewUpLoopMinNormal(
		const std::vector<Vertex3D>& loops,
		const std::vector<Vertex3D>& normals,
		std::vector<int>& tris)
	{
		int nhole = static_cast<int>(loops.size());
		if (nhole < 3)
			return false;
		std::function<ftype(int i, int j, int k)>DAngle;
		DAngle = [&loops, &normals](int i, int j, int k) {
			auto v1 = loops[j] - loops[i];
			auto v2 = loops[k] - loops[i];
			auto norm = v1 % v2;
			norm.normalize();
			auto ci = norm | normals[i];
			auto cj = norm | normals[j];
			auto ck = norm | normals[k];
			ftype ijk = static_cast<ftype>((ci + cj + ck) / 3.0);
			ftype dijk = static_cast<ftype>((ci - ijk)*(ci - ijk) + (cj - ijk)*(cj - ijk) + (ck - ijk)*(ck - ijk));
			return dijk;
		};
		std::function<ftype(int i, int j, int k, int l)>dihedral;
		dihedral = [&loops](int i, int j, int k, int l) {
			Vertex3D u = loops[i];
			Vertex3D v = loops[j];
			Vertex3D a = loops[k];
			Vertex3D b = loops[l];
			Vertex3D n0((v - u) % (a - v));
			Vertex3D n1((u - v) % (b - u));
			n0.normalize();
			n1.normalize();
			return static_cast<ftype>(acos(n0 | n1) * 180.0 / M_PI);
		};
		vvectori minimum_weight_index(nhole, vector<int>(nhole, -1));
		auto CAngle = [&minimum_weight_index, &dihedral, &nhole](int i, int j, int k) {
			if (minimum_weight_index[i][j] == -1)return static_cast<ftype>(0.0);
			if (minimum_weight_index[j][k] == -1)return static_cast<ftype>(0.0);
			ftype angle = static_cast<ftype>(0.0);
			if (i + 1 == j)
				angle = std::max(angle, dihedral(i, j, k, (i - 1 + nhole) % nhole));
			else
				angle = std::max(angle, dihedral(i, j, k, minimum_weight_index[i][j]));
			if (j + 1 == k)
				angle = std::max(angle, dihedral(j, k, i, (j - 1 + nhole) % nhole));
			else
				angle = std::max(angle, dihedral(j, k, i, minimum_weight_index[j][k]));
			if (i == 0 && k == nhole - 1)
				angle = std::max(angle, dihedral(k, i, j, (k - 1 + nhole) % nhole));
			return angle;
		};
		auto area = [&loops](int i, int j, int k) {
			Vertex3D a = loops[i];
			Vertex3D b = loops[j];
			Vertex3D c = loops[k];
			Vertex3D n((b - a) % (c - b));
			return static_cast<ftype>(0.5 * n.norm());
		};
		auto length = [&loops](int i, int j, int k) {
			ftype len1 = static_cast<ftype>((loops[i] - loops[j]).length());
			ftype len2 = static_cast<ftype>((loops[i] - loops[k]).length());
			ftype len3 = static_cast<ftype>((loops[k] - loops[j]).length());
			return len1 + len2 + len3;
		};
		vvectorf minimum_weight(nhole, vector<ftype>(nhole, 0));
		vector<int> ids;
		ids.reserve(nhole);
		for (int ic = 0; ic < nhole; ++ic)
			ids.push_back(ic);
		tris.clear();
		tris.reserve(nhole * 3);
		for (int j = 2; j < nhole; ++j) {
			for (int i = 0; i < nhole - j; ++i) {
				ftype min = FLT_MAX;
				int index = -1;
				int k = i + j;
				for (int m = i + 1; m < k; m++) {
					ftype val = minimum_weight[i][m] + minimum_weight[m][k] +
						area(i, m, k) + DAngle(i, m, k) + CAngle(i, m, k);
					if (val < min) {
						min = val; index = m;
					}
				}
				minimum_weight[i][k] = min;
				minimum_weight_index[i][k] = index;
			}
		}
		std::function<void(vector<int>& tris,
			const vector<int>& hole_id,
			HoleFiller::vvectori &minimum_weight_index,
			int begin, int end)>AddHoleToMesh;
		AddHoleToMesh = [&AddHoleToMesh](
			vector<int>& tris,
			const vector<int>& hole_id,
			vvectori &minimum_weight_index,
			int begin, int end) {
			if (end - begin > 1) {
				int cu = minimum_weight_index[begin][end];
				tris.push_back(hole_id[begin]);
				tris.push_back(hole_id[cu]);
				tris.push_back(hole_id[end]);
				AddHoleToMesh(tris, hole_id, minimum_weight_index, begin, cu);
				AddHoleToMesh(tris, hole_id, minimum_weight_index, cu, end);
			}
		};
		AddHoleToMesh(tris, ids, minimum_weight_index, 0, nhole - 1);
		return true;
	}


	/*
	 * @brief 保边界remesh
	 * @author lyc
	 * @date 2020.1.3
	 */
	bool HoleFiller::remeshTriangle(Triangle_mesh& mesh, ftype avelen, ftype angle) {
		mesh.request_edge_status();
		mesh.request_face_status();
		mesh.request_vertex_status();
		mesh.request_face_normals();
		mesh.request_vertex_normals();
		mesh.update_normals();
		OpenMesh::VPropHandleT<double>  mVertScale;
		mesh.add_property(mVertScale, "scale");
		Triangle_mesh::VertexOHalfedgeIter voh_iter;
		auto v_end = mesh.vertices_end();
		ftype msize = static_cast<ftype>(avelen);
		for (auto v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
			if (mesh.is_boundary(*v_it)) {
				mesh.status(*v_it).set_locked(true);
				double total = 0.0;
				int neiborsize = 0;
				voh_iter = mesh.voh_iter(*v_it);
				const auto& v = mesh.point(*v_it);
				for (; voh_iter.is_valid(); ++voh_iter)
				{
					const auto& p = mesh.point(mesh.to_vertex_handle(*voh_iter));
					total += (p - v).length();
					neiborsize++;
				}
				mesh.property(mVertScale, *v_it) = msize;
			}
		}
		vector<EdgeHandle> mBoundaryEH;
		vector<FaceHandle> mNewFillFace;
		auto e_end = mesh.edges_end();
		for (auto e_it = mesh.edges_begin(); e_it != e_end; ++e_it)
			if (mesh.is_boundary(*e_it))
				mesh.status(*e_it).set_locked(true);
		auto f_end = mesh.faces_end();
		mNewFillFace.reserve(mesh.n_faces());
		for (auto f_it = mesh.faces_begin(); f_it != f_end; ++f_it) {
			mNewFillFace.push_back(*f_it);
		}
		ftype  alpha = static_cast<ftype>(sqrt(angle));
		auto InCircumsphere = [](
			const Vertex3D & x,
			const Vertex3D & a,
			const Vertex3D & b,
			const Vertex3D & c) {
			auto ab = b - a;
			auto ac = c - a;
			double a00 = -2.0f * (ab | a);
			double a01 = -2.0f * (ab | b);
			double a02 = -2.0f * (ab | c);
			double b0 = a.sqrnorm() - b.sqrnorm();
			double a10 = -2.0f * (ac | a);
			double a11 = -2.0f * (ac | b);
			double a12 = -2.0f * (ac | c);
			double b1 = a.sqrnorm() - c.sqrnorm();
			double alpha = -(-a11 * a02 + a01 * a12 - a12 * b0 + b1 * a02 + a11 * b0 - a01 * b1)
				/ (-a11 * a00 + a11 * a02 - a10 * a02 + a00 * a12 + a01 * a10 - a01 * a12);
			double beta = (a10*b0 - a10 * a02 - a12 * b0 + a00 * a12 + b1 * a02 - a00 * b1)
				/ (-a11 * a00 + a11 * a02 - a10 * a02 + a00 * a12 + a01 * a10 - a01 * a12);
			double gamma = (-a11 * a00 - a10 * b0 + a00 * b1 + a11 * b0 + a01 * a10 - a01 * b1)
				/ (-a11 * a00 + a11 * a02 - a10 * a02 + a00 * a12 + a01 * a10 - a01 * a12);
			auto center = alpha * a + beta * b + gamma * c;
			return (x - center).sqrnorm() < (a - center).sqrnorm();
		};
		auto Relax = [&InCircumsphere](EdgeHandle eh, Triangle_mesh& mHoleMesh) {
			if (mHoleMesh.status(eh).locked())
				return false;
			HalfedgeHandle h0 = mHoleMesh.halfedge_handle(eh, 0);
			HalfedgeHandle h1 = mHoleMesh.halfedge_handle(eh, 1);
			auto u(mHoleMesh.point(mHoleMesh.to_vertex_handle(h0)));
			auto v(mHoleMesh.point(mHoleMesh.to_vertex_handle(h1)));
			auto a(mHoleMesh.point(mHoleMesh.to_vertex_handle(mHoleMesh.next_halfedge_handle(h0))));
			auto b(mHoleMesh.point(mHoleMesh.to_vertex_handle(mHoleMesh.next_halfedge_handle(h1))));
			if (InCircumsphere(a, u, v, b) || InCircumsphere(b, u, v, a)) {
				if (mHoleMesh.is_flip_ok(eh)) {
					mHoleMesh.flip(eh);
					return true;
				}
				else
					mHoleMesh.status(eh).set_selected(true);
			}
			return false;
		};
		auto Subdivide = [&alpha, &mesh, &mNewFillFace, &mVertScale, &mBoundaryEH, &Relax]() {
			bool status = false;
			size_t facenum = mNewFillFace.size();
			for (size_t i = 0; i < facenum; ++i) {
				HalfedgeHandle hh = mesh.halfedge_handle(mNewFillFace[i]);
				VertexHandle vi = mesh.to_vertex_handle(hh);
				VertexHandle vj = mesh.to_vertex_handle(mesh.prev_halfedge_handle(hh));
				VertexHandle vk = mesh.to_vertex_handle(mesh.next_halfedge_handle(hh));
				const auto& vip = mesh.point(vi);
				const auto& vjp = mesh.point(vj);
				const auto& vkp = mesh.point(vk);
				auto c = (vip + vjp + vkp) / 3.0f;
				const double vis = mesh.property(mVertScale, vi);
				const double vjs = mesh.property(mVertScale, vj);
				const double vks = mesh.property(mVertScale, vk);
				double sac = (vis + vjs + vks) / 3.0f;
				double dist_c_vi = (c - vip).length();
				double dist_c_vj = (c - vjp).length();
				double dist_c_vk = (c - vkp).length();
				if ((dist_c_vi + dist_c_vj + dist_c_vk) / 3.0f < sac)
					continue;
				if ((alpha * dist_c_vi > sac) &&
					(alpha * dist_c_vj > sac) &&
					(alpha * dist_c_vk > sac) &&
					(alpha * dist_c_vi > vis) &&
					(alpha * dist_c_vj > vjs) &&
					(alpha * dist_c_vk > vks)) {
					VertexHandle ch = mesh.add_vertex(c);
					mesh.split(mNewFillFace[i], ch);
					for (auto vfi = mesh.vf_iter(ch); vfi.is_valid(); ++vfi)
						if (*vfi != mNewFillFace[i])
							mNewFillFace.push_back(*vfi);
					for (auto vei = mesh.ve_iter(ch); vei.is_valid(); ++vei)
						mBoundaryEH.push_back(*vei);
					auto fei = mesh.fe_iter(mNewFillFace[i]);
					EdgeHandle  e0 = *fei; ++fei;
					EdgeHandle  e1 = *fei; ++fei;
					EdgeHandle  e2 = *fei; ++fei;
					Relax(e0, mesh);
					Relax(e1, mesh);
					Relax(e2, mesh);
					mesh.property(mVertScale, ch) = sac;
					status = true;
				}
			}
			return status;
		};
		auto CRelax = [&mBoundaryEH, &Relax](Triangle_mesh& mHoleMesh) {
			bool status = false;
			for (size_t i = 0; i < mBoundaryEH.size(); ++i)
				if (Relax(mBoundaryEH[i], mHoleMesh))
					status = true;
			return status;
		};
		for (int i = 0; i < 10; ++i) {
			bool is_subdivided = Subdivide();
			if (!is_subdivided)
				break;
			bool is_relaxed = CRelax(mesh);
			if (!is_relaxed)
				break;
		}
		for (auto e_it = mesh.edges_begin(); e_it != e_end; ++e_it)
			if (mesh.is_boundary(*e_it))
				mesh.status(*e_it).set_locked(false);
		mesh.remove_property(mVertScale);
		mesh.release_edge_status();
		mesh.release_face_status();
		mesh.release_vertex_status();
		mesh.release_face_normals();
		mesh.release_vertex_normals();
		return true;
	}
	/*
	* @brief 网格选择器
	* @param [in] mesh 待选的网格
	* @param [in] that 待选的初始边界
	* @param [in] n n圈领域
	* @param [out] spts 选择的点的id
	* @author lyc
	* @date 2020.1.3
	*/
	bool HoleFiller::selectVerts(Triangle_mesh& mesh, const std::vector<int>& that, int n, std::vector<int>&spts, bool issave)
	{
		int nvert = static_cast<int>(mesh.n_vertices());
		vector<bool> iscall(nvert, false);
		for (auto ip : that)
			iscall[ip] = true;
		vector<int> loops;
		loops.reserve(nvert);
		loops.insert(loops.end(), that.begin(), that.end());
		int i = 0;
		spts.clear();
		spts.reserve(nvert);
		if (issave)
			spts.insert(spts.end(), loops.begin(), loops.end());
		while ((i++) < n) {
			vector<int> curloop;
			curloop.reserve(nvert);
			for (auto ip : loops) {
				auto vh = VertexHandle(ip);
				auto eh = mesh.vv_end(vh);
				for (auto h = mesh.vv_begin(vh); h != eh; ++h) {
					if (iscall[h->idx()] == false)
					{
						curloop.push_back(h->idx());
						iscall[h->idx()] = true;
					}
				}
			}
			spts.insert(spts.end(), curloop.begin(), curloop.end());
			loops.swap(curloop);
		}
		return true;
	}

	bool HoleFiller::selectVerts(Triangle_mesh& mesh, const VertexHandle& that, int n, std::vector<VertexHandle>&spts, bool issave)
	{
		int nvert = static_cast<int>(mesh.n_vertices());
		vector<bool> iscall(nvert, false);
		iscall[that.idx()] = true;
		vector<VertexHandle> loops;
		loops.reserve(nvert);
		loops.push_back(that);
		int i = 0;
		spts.clear();
		spts.reserve(nvert);
		if (issave)
			spts.insert(spts.end(), loops.begin(), loops.end());
		while ((i++) < n) {
			vector<VertexHandle> curloop;
			curloop.reserve(nvert);
			for (auto ip : loops) {
				auto vh = VertexHandle(ip);
				auto eh = mesh.vv_end(vh);
				for (auto h = mesh.vv_begin(vh); h != eh; ++h) {
					if (iscall[h->idx()] == false)
					{
						curloop.push_back(*h);
						iscall[h->idx()] = true;
					}
				}
			}
			spts.insert(spts.end(), curloop.begin(), curloop.end());
			loops.swap(curloop);
		}
		return true;
	}

	bool HoleFiller::selectVerts(Triangle_mesh& mesh, const std::vector<VertexHandle>& thats,
		int n, std::vector<VertexHandle>&spts, bool issave)
	{
		int nvert = static_cast<int>(mesh.n_vertices());
		vector<bool> iscall(nvert, false);
		for(auto& that:thats)
			iscall[that.idx()] = true;
		vector<VertexHandle> loops;
		loops.reserve(nvert);
		loops.insert(loops.end(),thats.begin(),thats.end());
		int i = 0;
		spts.clear();
		spts.reserve(nvert);
		if (issave)
			spts.insert(spts.end(), loops.begin(), loops.end());
		while ((i++) < n) {
			vector<VertexHandle> curloop;
			curloop.reserve(nvert);
			for (auto ip : loops) {
				auto vh = VertexHandle(ip);
				auto eh = mesh.vv_end(vh);
				for (auto h = mesh.vv_begin(vh); h != eh; ++h) {
					if (iscall[h->idx()] == false)
					{
						curloop.push_back(*h);
						iscall[h->idx()] = true;
					}
				}
			}
			spts.insert(spts.end(), curloop.begin(), curloop.end());
			loops.swap(curloop);
		}
		return true;
	}

	void HoleFiller::extraBoundary(Triangle_mesh* mesh, std::vector<Point3f>& boundary, std::vector<int>& ids)
	{
		vector<vector<int>> aids;
		vector<bool> iscall(mesh->n_vertices(), false);
		for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
		{
			if (!mesh->is_boundary(*iv) || iscall[iv->idx()] != false) {
				continue;
			}
			auto vh = *iv;
			iscall[vh.idx()] = true;
			auto h = mesh->halfedge_handle(vh);
			if (!mesh->is_boundary(h))
				h = mesh->opposite_halfedge_handle(h);
			auto nh = h;
			ids.clear();
			do {
				auto vh = mesh->from_vertex_handle(nh);
				ids.push_back(vh.idx());
				nh = mesh->next_halfedge_handle(nh);
				iscall[vh.idx()] = true;
			} while (nh != h);
			aids.push_back(ids);
		}
		ids.clear();
		int ik = 0;
		int size = 0;
		for (int i = 0; i < aids.size(); ++i)
		{
			if (size < aids[i].size())
			{
				size = (int)aids[i].size();
				ik = i;
			}
		}
		ids = aids[ik];
		boundary.clear();
		for (auto i : ids) {
			boundary.push_back(mesh->point(VertexHandle(i)));
		}
	}

	void HoleFiller::extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& aids, std::vector<int>& ids, int& ik)
	{
		aids.clear();
		vector<bool> iscall(mesh->n_vertices(), false);
		for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
		{
			if (!mesh->is_boundary(*iv) || iscall[iv->idx()] != false) {
				continue;
			}
			auto vh = *iv;
			iscall[vh.idx()] = true;
			auto h = mesh->halfedge_handle(vh);
			if (!mesh->is_boundary(h))
				h = mesh->opposite_halfedge_handle(h);
			auto nh = h;
			ids.clear();
			do {
				auto vh = mesh->from_vertex_handle(nh);
				ids.push_back(vh.idx());
				nh = mesh->next_halfedge_handle(nh);
				iscall[vh.idx()] = true;
			} while (nh != h);
			aids.push_back(ids);
		}
		ids.clear();
		ik = -1;
		int size = 0;
		for (int i = 0; i < aids.size(); ++i)
		{
			if (size < aids[i].size())
			{
				size = (int)aids[i].size();
				ik = i;
			}
		}
		if (aids.size() > 0 && ik < aids.size())
		{
			ids = aids[ik];
		}
	}

	void HoleFiller::extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& caids)
	{
		auto GetReatEle = [](
			std::vector<int> &pcu,
			std::set<int> &sct,
			std::vector<int> &repeat)
		{
			vector<int> cur = pcu;
			sort(cur.begin(), cur.end());
			for (int ic = 0; ic < cur.size(); ++ic)
			{
				for (int ik = ic + 1; ik < cur.size(); ++ik)
				{
					if (cur[ic] == cur[ik])
					{
						if (sct.count(cur[ic]) == 0)
						{
							sct.insert(cur[ic]);
							repeat.push_back(cur[ic]);
						}
					}
					else
					{
						break;
					}
				}
			}
		};
		auto decoupHole = [&GetReatEle](
			const std::vector<std::vector<int>>&vchole,
			std::vector<std::vector<int>>& vdecode) {
			int nhole = static_cast<int>(vchole.size());
			vdecode.clear();
			vdecode.reserve(nhole);
			for (int ic = 0; ic < nhole; ++ic)
			{
				vector<int> pcu = vchole[ic];
				vector<int> repeat;
				set<int> sct;
				GetReatEle(pcu, sct, repeat);
				if (repeat.size() == 0)
				{
					vdecode.push_back(pcu);
					continue;
				}
				stack<int> dcoup;
				stack<int> mark;
				for (int ik = 0; ik < pcu.size(); ++ik)
				{
					dcoup.push(pcu[ik]);
					if (sct.count(pcu[ik]) != 0)
					{
						if (mark.size() == 0 || mark.top() != pcu[ik])
						{
							mark.push(pcu[ik]);
						}
						else
						{
							vector<int> cir;
							do
							{
								cir.push_back(dcoup.top());
								dcoup.pop();
							} while (dcoup.top() != pcu[ik]);
							reverse(cir.begin(), cir.end());
							vdecode.push_back(cir);
						}
					}
				}
				vector<int> cir;
				while (dcoup.size() != 0)
				{
					cir.push_back(dcoup.top());
					dcoup.pop();
				}
				std::reverse(cir.begin(), cir.end());
				vector<int> crepeat;
				set<int> csct;
				GetReatEle(cir, csct, crepeat);
				if (csct.size() == 0)
					vdecode.push_back(cir);
			}
		};
		std::vector<std::vector<int>> aids;
		aids.clear();
		vector<bool> iscall(mesh->n_vertices(), false);
		vector<int> ids;
		for (auto iv = mesh->vertices_begin(); iv != mesh->vertices_end(); ++iv)
		{
			if (!mesh->is_boundary(*iv) || iscall[iv->idx()] != false) {
				continue;
			}
			auto vh = *iv;
			iscall[vh.idx()] = true;
			auto h = mesh->halfedge_handle(vh);
			if (!mesh->is_boundary(h))
				h = mesh->opposite_halfedge_handle(h);
			auto nh = h;
			ids.clear();
			ids.reserve(mesh->n_vertices());
			do {
				auto vh = mesh->from_vertex_handle(nh);
				ids.push_back(vh.idx());
				nh = mesh->next_halfedge_handle(nh);
				iscall[vh.idx()] = true;
			} while (nh != h);
			aids.push_back(ids);
		}
		decoupHole(aids, caids);
	}

	struct Triple
	{
		Triple() {};
		Triple(VertexHandle v, HoleFiller::ftype weight, int degree)
			: vertex_(v), weight_(weight), degree_(degree)
		{}
		VertexHandle vertex_;
		HoleFiller::ftype weight_;
		int degree_;
	};

	static HoleFiller::ftype getAverageEdgeLength(Triangle_mesh &mesh)
	{
		MeshDoctor::ftype average_edge_length = 0.0f;
		for (Triangle_mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
			average_edge_length += mesh.calc_edge_length(*e_it);
		MeshDoctor::ftype edgeNum = (MeshDoctor::ftype)mesh.n_edges();
		average_edge_length /= edgeNum;
		return average_edge_length;
	}


	//! \fn计算网格的cotan laplace矩阵
		/*!
		  \param[in]  mesh 网格数据
		  \param[out] L    laplace矩阵
		  \param[out] M    质量矩阵，一般M*L表示拉普拉斯矩阵
		  \return
		  - true  成功
		  - false 失败
		*/
	static bool CotLaplaceMatrix(
		const Triangle_mesh& mesh,
		Eigen::SparseMatrix<MeshDoctor::ftype>& L,
		Eigen::SparseMatrix<MeshDoctor::ftype>& M)
	{
		using Point = OpenMesh::Vec3f;
		int vsize = static_cast<int>(mesh.n_vertices());
		if (vsize < 3)
			return false;
		if (!mesh.is_trimesh())
			return false;
		L.resize(vsize, vsize);
		M.resize(vsize, vsize);
		std::vector<Eigen::Triplet<MeshDoctor::ftype, int> > IJV;
		std::vector<Eigen::Triplet<MeshDoctor::ftype, int> > MTrip;
		IJV.reserve(8 * mesh.n_vertices());
		MTrip.reserve(mesh.n_vertices());
		auto v_end = mesh.vertices_end();
		auto  TrigCut = [](MeshDoctor::ftype value)
		{
			if (value < -1)
				value = -1.0f;
			else if (value > 1)
				value = 1.0f;
			return value;
		};
		auto IsObtuse = [](const Point& p0, const Point& p1, const Point& p2)
		{
			MeshDoctor::ftype a0 = ((p1 - p0) | (p2 - p0));
			if (a0 < 0.0f)
				return 1;
			else {
				MeshDoctor::ftype a1 = ((p0 - p1) | (p2 - p1));
				if (a1 < 0.0f)
					return 2;
				else {
					MeshDoctor::ftype a2 = ((p0 - p2) | (p1 - p2));
					if (a2 < 0.0f)
						return 3;
					else return 0;
				}
			}
			return 0;
		};
		for (auto v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		{
			MeshDoctor::ftype w_ii = 0;
			MeshDoctor::ftype area = 0.0f;
			auto voh_it = mesh.cvoh_iter(*v_it);
			if (!voh_it->is_valid())
				return false;
			for (; voh_it.is_valid(); ++voh_it)
			{
				if (mesh.is_boundary(mesh.edge_handle(*voh_it)))
					continue;
				HalfedgeHandle lefthh = mesh.prev_halfedge_handle(*voh_it);
				HalfedgeHandle righthh = mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(*voh_it));
				const Point& p0 = mesh.point(*v_it);
				const Point& p1 = mesh.point(mesh.to_vertex_handle(*voh_it));
				const Point& p2 = mesh.point(mesh.from_vertex_handle(lefthh));
				const Point& p3 = mesh.point(mesh.to_vertex_handle(righthh));
				MeshDoctor::ftype alpha = std::acos(TrigCut((p0 - p2).normalize() | (p1 - p2).normalize()));
				MeshDoctor::ftype beta = std::acos(TrigCut((p0 - p3).normalize() | (p1 - p3).normalize()));
				MeshDoctor::ftype cotw = 0.0f;
				if (!util::isEqual(alpha, (MeshDoctor::ftype)M_PI_2))
					cotw += 1.0f / std::tan(alpha);
				if (!util::isEqual(beta, (MeshDoctor::ftype)M_PI_2))
					cotw += 1.0f / std::tan(beta);
				if (util::isZero(cotw) || std::isinf(cotw) || std::isnan(cotw))
					continue;
				const int obt = IsObtuse(p0, p1, p2);
				if (obt == 0) {
					MeshDoctor::ftype gamma = acos(TrigCut((p0 - p1).normalize() | (p2 - p1).normalize()));
					MeshDoctor::ftype tmp = 0.0f;
					if (!util::isEqual(alpha, (MeshDoctor::ftype)M_PI_2))
						tmp += (p0 - p1).sqrnorm()*1.0f / tan(alpha);
					if (!util::isEqual(gamma, (MeshDoctor::ftype)M_PI_2))
						tmp += (p0 - p2).sqrnorm()*1.0f / tan(gamma);
					if (util::isZero(tmp) || std::isinf(tmp) || std::isnan(tmp))
						continue;
					area += 0.125f*(tmp);
				}
				else {
					if (obt == 1)
						area += ((p1 - p0) % (p2 - p0)).norm() * 0.5f * 0.5f;
					else
						area += ((p1 - p0) % (p2 - p0)).norm() * 0.5f * 0.25f;
				}
				w_ii -= cotw;
				VertexHandle vh = mesh.to_vertex_handle(*voh_it);
				IJV.push_back(Eigen::Triplet<MeshDoctor::ftype, int>(v_it->idx(), vh.idx(), cotw));
			}
			IJV.push_back(Eigen::Triplet<MeshDoctor::ftype, int>(v_it->idx(), v_it->idx(), w_ii));
			MTrip.push_back(Eigen::Triplet<MeshDoctor::ftype, int>(v_it->idx(), v_it->idx(), 0.5f / area));
		}
		L.setFromTriplets(IJV.begin(), IJV.end());
		M.setFromTriplets(MTrip.begin(), MTrip.end());
		return true;
	}


	bool HoleFiller::nlocalFairMesh(Triangle_mesh& mesh, vector<VertexHandle>& handle, int c)
	{
		using namespace std;
		vector<VertexHandle> freevertx;
		freevertx.reserve(mesh.n_vertices());
		for (auto iv : handle) {
			if (iv.idx() < mesh.n_vertices() && !mesh.is_boundary(iv) && mesh.is_manifold(iv)) {
				freevertx.push_back(iv);
			}
		}
		Eigen::Matrix<ftype, Eigen::Dynamic, 3> b(freevertx.size(), 3);
		Eigen::Matrix<ftype, Eigen::Dynamic, 3> x(freevertx.size(), 3);
		b.setZero();
		Eigen::SparseMatrix<ftype>  L, Mat, M, LU;
		CotLaplaceMatrix(mesh, L, M);
		if (0 == c)
			Mat = M * L;
		else if (1 == c)
			Mat = M * L* M *L;
		else if (2 == c) {
			Mat = M * L*M*L*M*L;
		}
		else {
			return false;
		}
		int id = 0;
		std::map<VertexHandle, int> idmap;
		for (auto vh : freevertx)
			idmap[vh] = id++;

		std::vector<Eigen::Triplet<ftype, int> > IJV;
		for (int k = 0; k < Mat.outerSize(); ++k) {
			for (Eigen::SparseMatrix<ftype>::InnerIterator it(Mat, k); it; ++it) {
				ftype v = it.value();
				VertexHandle rowvh = VertexHandle((int)it.row());
				VertexHandle colvh = VertexHandle((int)it.col());
				if (idmap.find(rowvh) == idmap.end())
					continue;
				if (idmap.find(colvh) != idmap.end())
					IJV.push_back(Eigen::Triplet<ftype, int>(idmap[rowvh], idmap[colvh], v));
				else {
					auto p = mesh.point(colvh);
					b(idmap[rowvh], 0) += -1.0f*v*p[0];
					b(idmap[rowvh], 1) += -1.0f*v*p[1];
					b(idmap[rowvh], 2) += -1.0f*v*p[2];
				}
			}
		}
		LU.resize(freevertx.size(), freevertx.size());
		LU.setFromTriplets(IJV.begin(), IJV.end());
		LU.makeCompressed();
		Eigen::SparseLU<Eigen::SparseMatrix<ftype>, Eigen::COLAMDOrdering<int> >   solver;
		solver.compute(LU);
		if (solver.info() != Eigen::Success) {
			std::cout << "faild compute:" << solver.info() << std::endl;
			return false;
		}
		x = solver.solve(b);
		if (solver.info() != Eigen::Success)
		{
			std::cout << "faild solve:" << solver.info() << std::endl;
			return false;
		}
		auto v_end = mesh.vertices_end();
		//检查误差
		ftype ave_edge = 10 * getAverageEdgeLength(mesh);
		for (auto v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		{
			if (idmap.find(*v_it) == idmap.end())
				continue;
			int index = idmap[*v_it];
			Point3f& p = mesh.point(*v_it);
			Point3f ep;
			ep[0] = x(index, 0);
			ep[1] = x(index, 1);
			ep[2] = x(index, 2);
			ftype cfd = (p - ep).length();
			if (cfd >= ave_edge)
				return false;
		}
		for (auto v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		{
			if (idmap.find(*v_it) == idmap.end())
				continue;
			int index = idmap[*v_it];
			Point3f& p = mesh.point(*v_it);
			p[0] = x(index, 0);
			p[1] = x(index, 1);
			p[2] = x(index, 2);
		}
		return true;
	}

}