#include "MeshDoctor.h"
#include <iostream>
#include <Eigen/Eigen>
#include <stack>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
using namespace std;

#ifdef M_PI
#undef M_PI
#define M_PI 3.1415926535897932384626433832795f
#endif
#define AG_10 0.17453292519943295769236907684886f
#define AG_145 2.5307274153917778865393516143085f
#define AG_60 1.0471975511965977461542144610932f

static OpenMesh::VPropHandleT<int> CompIndex;
static const OpenMesh::Vec3uc RED = { 255,0,0 };
static const OpenMesh::Vec3uc BLUE = { 27,171,241 };

template<class TVert>
class crepeator
{
public:
	bool operator()(const TVert& _v0, const TVert& _v1) const
	{
		if (fabs(_v0[0] - _v1[0]) <= eps_)
		{
			if (fabs(_v0[1] - _v1[1]) <= eps_)
			{
				return (_v0[2] < _v1[2] - eps_);
			}
			else return (_v0[1] < _v1[1] - eps_);
		}
		else return (_v0[0] < _v1[0] - eps_);
	}
public:
	static float eps_;
};
template<class TVert>
float crepeator<TVert>::eps_;
namespace GeniusCore
{
	using ftype = float;
	class  CMeanValueParameterizer
	{	
	public:
		typedef Eigen::Matrix<ftype, Eigen::Dynamic, 1> EigenVector;
		typedef Eigen::SparseMatrix<ftype>              EigenMatrix;
		typedef Eigen::Triplet<ftype, int>              Triplet;
		typedef Eigen::BiCGSTAB<EigenMatrix, Eigen::IncompleteLUT<ftype>>  Solver;
		using Point = OpenMesh::Vec3f;
		using PointXYd = OpenMesh::Vec2f;
		using Vertex3D = Point;
	public:
		CMeanValueParameterizer(Triangle_mesh& mesh);
		~CMeanValueParameterizer();
		bool Parameterize();
	private:
		inline ftype FixSine(ftype sine)
		{
			if (sine >= 1.0)
				return 1.0;
			else if (sine <= -1.0)
				return -1.0;
			else
				return sine;
		}
		//
		bool BorderParameterize(const HalfedgeHandle& hh, const ftype borderlen);
		void InitializeSystemBorder(EigenMatrix& A, std::vector<Triplet>& Atriplet,
			EigenVector& Bu, EigenVector& Bv, HalfedgeHandle& hh);
		bool SetInnerVertexRelations(EigenMatrix& A, std::vector<Triplet>& Atriplet,
			EigenVector& Bu, EigenVector& Bv);
		ftype ComputeWeight(const VertexHandle& vert, const VertexHandle& around);
		ftype ComputeAngleRad(const Point& P, const Point& Q, const Point& R);
	private:
		Triangle_mesh&                   mTriMesh;
		OpenMesh::VPropHandleT<PointXYd> mUVProperty;
		OpenMesh::VPropHandleT<bool>     mVertStatus;
		std::string                      mPropertyName;
	};
}


OpenMesh::FPropHandleT<MeshDoctor::ftype> MeshDoctor::MinAngle; //最小角
OpenMesh::FPropHandleT<MeshDoctor::ftype> MeshDoctor::MaxAngle; //最大角
OpenMesh::FPropHandleT<MeshDoctor::ftype> MeshDoctor::AspectRatio; //纵横比
OpenMesh::FPropHandleT<MeshDoctor::ftype> MeshDoctor::SkewnessEquiangle; //偏斜等角
OpenMesh::VPropHandleT<MeshDoctor::Point3f> MeshDoctor::SmoothVertex;//特征度
OpenMesh::FPropHandleT<MeshDoctor::ftype> MeshDoctor::SkewnessEquiarea; //偏斜比
OpenMesh::EPropHandleT<MeshDoctor::ftype> MeshDoctor::AreaRatio; //偏斜比
OpenMesh::VPropHandleT<bool> MeshDoctor::Genius; //偏斜比


using Matrix = Eigen::Matrix<MeshDoctor::ftype, Eigen::Dynamic, Eigen::Dynamic>;
static  MeshDoctor::ftype fix(MeshDoctor::ftype val) {
	val = val > 1 ? 1 : val;
	val = val < -1 ? -1 : val;
	return val;
}

static MeshDoctor::ftype fixCot(MeshDoctor::ftype v) {
	MeshDoctor::ftype bound = 5.6713f; // 10 degrees
	return (v < -bound ? -bound : (v > bound ? bound : v));
}

static MeshDoctor::ftype fixCos(MeshDoctor::ftype v) {
	MeshDoctor::ftype bound = 0.985f; // 10 degrees
	return (v < -bound ? -bound : (v > bound ? bound : v));
}

static MeshDoctor::ftype angle(const MeshDoctor::Point3f& v0,const MeshDoctor::Point3f& v1) {
	return atan2(OpenMesh::norm(OpenMesh::cross(v0, v1)), OpenMesh::dot(v0, v1));
}

static  MeshDoctor::ftype sin(const MeshDoctor::Point3f& v0, const MeshDoctor::Point3f& v1)
{
	using namespace OpenMesh;
	return norm(cross(v0, v1)) / (norm(v0) * norm(v1));
}

static  MeshDoctor::ftype cos(const MeshDoctor::Point3f& v0, const MeshDoctor::Point3f& v1)
{
	using namespace OpenMesh;
	return dot(v0, v1) / (norm(v0) * norm(v1));
}

static  MeshDoctor::ftype cotan(const MeshDoctor::Point3f& v0, const MeshDoctor::Point3f& v1)
{
	using namespace OpenMesh;
	return fixCot(dot(v0, v1) / norm(cross(v0, v1)));
}

MeshDoctor::MeshDoctor(TriMesh& msh):
	mesh(msh)
{
	if (!mesh.has_face_colors())
		mesh.request_face_colors();
	if (!mesh.has_vertex_colors())
		mesh.request_vertex_colors();
}

MeshDoctor::~MeshDoctor()
{
	if (MinAngle.is_valid())
		mesh.remove_property(MinAngle);
	if (MaxAngle.is_valid())
		mesh.remove_property(MaxAngle);
	if (AspectRatio.is_valid())
		mesh.remove_property(AspectRatio);
	if (SkewnessEquiangle.is_valid())
		mesh.remove_property(SkewnessEquiangle);
	if (SmoothVertex.is_valid())
		mesh.remove_property(SmoothVertex);
	if (SkewnessEquiarea.is_valid())
		mesh.remove_property(SkewnessEquiarea);
	if (AreaRatio.is_valid())
		mesh.remove_property(AreaRatio);
	if (mesh.has_face_normals())
		mesh.release_face_normals();
	if (Genius.is_valid())
		mesh.remove_property(Genius);
	if (mesh.has_face_colors())
		mesh.release_face_colors();
	if (mesh.has_vertex_colors())
		mesh.release_vertex_colors();
}
void MeshDoctor::updateAngle()
{
	if (!MinAngle.is_valid())
		mesh.add_property(MinAngle, "minAngle");
	if (!MaxAngle.is_valid())
		mesh.add_property(MaxAngle, "maxAngle");
	int nface =(int) mesh.n_faces();
	auto points = mesh.points();
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh = FaceHandle(i);
		auto fv = mesh.fv_begin(fh);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx(); 
		ftype a = (points[o]-points[t]).length();
		ftype b = (points[o]-points[h]).length();
		ftype c = (points[h]-points[t]).length();
		ftype ag = acos(fix((b*b + c * c - a * a) / (2 * b*c)));
		ftype ab = acos(fix((a*a + c * c - b * b) / (2 * a*c)));
		ftype ac = acos(fix((a*a + b * b - c * c) / (2 * a*b)));
		ftype minAgle = ag < ab ? ag : ab;
		minAgle = minAgle < ac ? minAgle : ac;
		ftype maxAngle = ag > ab ? ag : ab;
		maxAngle = maxAngle > ac ? maxAngle : ac;
		mesh.property(MinAngle,fh)= minAgle;
		mesh.property(MaxAngle,fh)= maxAngle;
	}
}

void MeshDoctor::calcAngle(const FaceHandle& fh, ftype& minAgle, ftype& maxAngle)
{
	auto points = mesh.points();
	auto fv = mesh.fv_begin(fh);
	int o = fv->idx(); ++fv;
	int t = fv->idx(); ++fv;
	int h = fv->idx();
	ftype a = (points[o] - points[t]).length();
	ftype b = (points[o] - points[h]).length();
	ftype c = (points[h] - points[t]).length();
	ftype ag = acos(fix((b*b + c * c - a * a) / (2 * b*c)));
	ftype ab = acos(fix((a*a + c * c - b * b) / (2 * a*c)));
	ftype ac = acos(fix((a*a + b * b - c * c) / (2 * a*b)));
	minAgle = ag < ab ? ag : ab;
	minAgle = minAgle < ac ? minAgle : ac;
	maxAngle = ag > ab ? ag : ab;
	maxAngle = maxAngle > ac ? maxAngle : ac;
}

void MeshDoctor::updateAspectRatio()
{
	if (!AspectRatio.is_valid())
		mesh.add_property(AspectRatio);
	int nface = (int)mesh.n_faces();
	auto points = mesh.points();
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh = FaceHandle(i);
		auto fv = mesh.fv_begin(fh);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		ftype a = (points[o] - points[t]).length();
		ftype b = (points[o] - points[h]).length();
		ftype c = (points[h] - points[t]).length();
		ftype minlen;
		ftype maxlen;
		minlen = a < b ? a : b;
		minlen = minlen < c ? minlen : c;
		maxlen = a > b ? a : b;
		maxlen = maxlen > c ? maxlen : c;
		if (minlen < FLT_EPSILON)
			minlen = FLT_EPSILON;
		mesh.property(AspectRatio, fh) = maxlen/minlen;
	}
}

MeshDoctor::ftype MeshDoctor::calcAspectRatio(const FaceHandle& fh) {
	auto points = mesh.points();
	auto fv = mesh.fv_begin(fh);
	int o = fv->idx(); ++fv;
	int t = fv->idx(); ++fv;
	int h = fv->idx();
	ftype a = (points[o] - points[t]).length();
	ftype b = (points[o] - points[h]).length();
	ftype c = (points[h] - points[t]).length();
	ftype minlen;
	ftype maxlen;
	minlen = a < b ? a : b;
	minlen = minlen < c ? minlen : c;
	maxlen = a > b ? a : b;
	maxlen = maxlen > c ? maxlen : c;
	if (minlen < FLT_EPSILON)
		minlen = FLT_EPSILON;
	return maxlen / minlen;
}

void MeshDoctor::updateSkewnessEquiarea()
{
	if (mesh.n_faces() == 0)
		return;
	if (!SkewnessEquiarea.is_valid())
		mesh.add_property(SkewnessEquiarea);
	int nface = (int)mesh.n_faces();
	auto points = mesh.points();
	vector<ftype> area;
	area.resize(nface);
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh = FaceHandle(i);
		auto fv = mesh.fv_begin(fh);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		ftype a = (points[o] - points[t]).length();
		ftype b = (points[o] - points[h]).length();
		ftype c = (points[h] - points[t]).length();
		ftype p = (a + b + c) / 2;
		ftype ar = static_cast<ftype>(sqrt(p * (p - a) * (p - b) * (p - c)));
		area[i] = ar;
	}
	ftype avg = ftype(0);
	for (auto ia : area)
		avg += ia;
	avg /= (ftype)nface;
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh = FaceHandle(i);
		mesh.property(SkewnessEquiangle, fh) = (area[i]-avg)/avg;
	}
}

MeshDoctor::ftype MeshDoctor::calcSkewnessEquiarea(const FaceHandle& fh, int cn)
{
	vector<FaceHandle> fhs;
	selectFace(fh, cn, fhs, true);
	int nface = (int)fhs.size();
	vector<ftype> area;
	area.resize(nface);
	auto points = mesh.points();
	for (int i = 0; i < nface; ++i) {
		FaceHandle& fh = fhs[i];
		auto fv = mesh.fv_begin(fh);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		ftype a = (points[o] - points[t]).length();
		ftype b = (points[o] - points[h]).length();
		ftype c = (points[h] - points[t]).length();
		ftype p = (a + b + c) / 2;
		ftype ar = static_cast<ftype>(sqrt(p * (p - a) * (p - b) * (p - c)));
		area[i] = ar;
	}
	ftype avg = ftype(0);
	for (auto ia : area)
		avg += ia;
	avg /= (ftype)nface;
	return (area[0] - avg) / avg;

}

void MeshDoctor::updateSkewnessEquiangle()
{
	if (mesh.n_faces() == 0)
		return;
	if (!SkewnessEquiarea.is_valid())
		mesh.add_property(SkewnessEquiangle);
	int nface = (int)mesh.n_faces();
	updateAngle();
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh(i);
		ftype qmin = mesh.property(MinAngle, fh);
		ftype qmax = mesh.property(MaxAngle, fh);
		ftype qe = (qmax - AG_60) / (M_PI - AG_60);
		ftype cq = (AG_60 - qmin) / AG_60;
		mesh.property(SkewnessEquiangle, fh) = max(qe, cq);
	}
}

MeshDoctor::ftype MeshDoctor::calcSkewnessEquiangle(const FaceHandle& fh)
{
	ftype qmax, qmin;
	calcAngle(fh, qmin, qmax);
	ftype qe = (qmax - AG_60) / (M_PI - AG_60);
	ftype cq = (AG_60 - qmin) / AG_60;
	return max(qe, cq);
}

void MeshDoctor::updateAreaRatio()
{
	if (!AreaRatio.is_valid())
		mesh.add_property(AreaRatio);
	int nface = (int)mesh.n_faces();
	auto points = mesh.points();
	vector<ftype> area;
	area.resize(nface);
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh = FaceHandle(i);
		auto fv = mesh.fv_begin(fh);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		ftype a = (points[o] - points[t]).length();
		ftype b = (points[o] - points[h]).length();
		ftype c = (points[h] - points[t]).length();
		ftype p = (a + b + c) / 2;
		ftype ar = static_cast<ftype>(sqrt(p * (p - a) * (p - b) * (p - c)));
		if (ar < FLT_EPSILON)
			ar = FLT_EPSILON;
		area[i] = ar;
	}
	int nedge = (int)mesh.n_edges();
	for (int i = 0; i < nedge;++i) {
		auto ie = EdgeHandle(i);
		auto h0 = mesh.halfedge_handle(ie, 0);
		auto h1 = mesh.halfedge_handle(ie, 1);
		ftype a0=FLT_EPSILON, a1= FLT_EPSILON;
		if (!mesh.is_boundary(h0)) {
			int f = mesh.face_handle(h0).idx();
			a0 = area[f];
		}
		if (!mesh.is_boundary(h1)) {
			int f = mesh.face_handle(h1).idx();
			a1 = area[f];
		}
		mesh.property(AreaRatio, ie) = min(a0, a1) / max(a0, a1);
	}

}

MeshDoctor::ftype MeshDoctor::calcFaceArea(const FaceHandle& fh)
{
	if (!mesh.is_valid_handle(fh))
		return 0;
	auto points = mesh.points();
	auto fv = mesh.fv_begin(fh);
	int o = fv->idx(); ++fv;
	int t = fv->idx(); ++fv;
	int h = fv->idx();
	ftype a = (points[o] - points[t]).length();
	ftype b = (points[o] - points[h]).length();
	ftype c = (points[h] - points[t]).length();
	ftype p = (a + b + c) / 2;
	ftype ar = static_cast<ftype>(sqrt(p * (p - a) * (p - b) * (p - c)));
	if (ar < FLT_EPSILON)
		ar = FLT_EPSILON;
	return ar;
}

MeshDoctor::ftype MeshDoctor::calcCotanWeight(TriMesh& mesh, const EdgeHandle& eh)
{
	ftype weight = 0.0f;
	if (!mesh.is_valid_handle(eh))
		return weight;
	using namespace OpenMesh;
	const auto h0 = mesh.halfedge_handle(eh, 0);
	const auto h1 = mesh.halfedge_handle(eh, 1);
	const Point3f p0 = mesh.point(mesh.to_vertex_handle(h0));
	const Point3f p1 = mesh.point(mesh.to_vertex_handle(h1));
	if (!mesh.is_boundary(h0)){
		const Point3f p2 =
			(Point3f)mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(h0)));
		const Point3f d0 = p0 - p2;
		const Point3f d1 = p1 - p2;
	    ftype area = norm(cross(d0, d1));
		area = area < FLT_EPSILON ? FLT_EPSILON : area;
		const ftype cot = dot(d0, d1) / area;
		weight += fixCot(cot);
	}
	if (!mesh.is_boundary(h1)) {
		const Point3f p2 =
			(Point3f)mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(h1)));
		const Point3f d0 = p0 - p2;
		const Point3f d1 = p1 - p2;
		ftype area = norm(cross(d0, d1));
		area = area < FLT_EPSILON ? FLT_EPSILON : area;
		const ftype cot = dot(d0, d1) / area;
		weight += fixCot(cot);
	}
	return weight;
}

MeshDoctor::ftype MeshDoctor::calcVoronoiArea(TriMesh& mesh, const VertexHandle& v)
{
	ftype area(0.0f);
	using namespace OpenMesh;
	if (!mesh.is_isolated(v))
	{
		HalfedgeHandle h0, h1, h2;
		Point3f p, q, r, pq, qr, pr;
		ftype dotp, dotq, dotr, triArea;
		ftype cotq, cotr;
		for (auto h = mesh.voh_begin(v); h != mesh.voh_end(v);++h) {
			h0 = *h;
			h1 = mesh.next_halfedge_handle(h0);
			h2 = mesh.next_halfedge_handle(h1);
			if (mesh.is_boundary(h0))
				continue;
			p = mesh.point(mesh.to_vertex_handle(h2));
			q = mesh.point(mesh.to_vertex_handle(h0));
			r = mesh.point(mesh.to_vertex_handle(h1));
			(pq = q) -= p;
			(qr = r) -= q;
			(pr = r) -= p;
			triArea = norm(cross(pq, pr));
			triArea = triArea < FLT_EPSILON ? FLT_EPSILON : triArea;
			dotp = dot(pq, pr);
			dotq = -dot(qr, pq);
			dotr = dot(qr, pr);
			if (dotp < 0.0f)
				area += 0.25f * triArea;
			else if (dotq < 0.0f || dotr < 0.0f)
				area += 0.125f * triArea;
			else
			{
				cotq = dotq / triArea;
				cotr = dotr / triArea;
				area += 0.125f * (sqrnorm(pr) * fixCot(cotq) + sqrnorm(pq) * fixCot(cotr));
			}
		}
	}
	area = area < FLT_EPSILON ? FLT_EPSILON : area;
	return area;
}

MeshDoctor::ftype MeshDoctor::calcVoronoiAreaBarycentric(TriMesh& mesh, const VertexHandle& v)
{
	ftype area(0.0f);
	using namespace OpenMesh;
	if (!mesh.is_isolated(v))
	{
		const Point3f p = mesh.point(v);
		HalfedgeHandle h0, h1;
		Point3f q, r, pq, pr;
		for (auto h = mesh.voh_begin(v); h != mesh.voh_end(v); ++h)
		{
			if (mesh.is_boundary(*h))
				continue;
			h0 = *h;
			h1 = mesh.next_halfedge_handle(h0);
			pq = mesh.point(mesh.to_vertex_handle(h0));
			pq -= p;
			pr = mesh.point(mesh.to_vertex_handle(h1));
			pr -= p;
			area += norm(cross(pq, pr)) / 3.0f;
		}
	}
	return area;
}

MeshDoctor::Point3f MeshDoctor::calcLaplace(TriMesh& mesh, const VertexHandle& v)
{
	Point3f laplace(0.0f, 0.0f, 0.0f);
	if (!mesh.is_isolated(v))
	{
		ftype weight, sumWeights(0.0f);
		for (auto h = mesh.voh_begin(v); h != mesh.voh_end(v); ++h){
			weight = calcCotanWeight(mesh, mesh.edge_handle(*h));
			sumWeights += weight;
			laplace += weight * mesh.point(mesh.to_vertex_handle(*h));
		}
		laplace -= sumWeights * mesh.point(v);
		laplace /= ftype(2.0f) * calcVoronoiArea(mesh, v);
	}
	return laplace;
}

MeshDoctor::ftype MeshDoctor::calcAngleSum(TriMesh& mesh, const VertexHandle& v)
{
	ftype angles(0.0f);
	using namespace OpenMesh;
	if (!mesh.is_boundary(v))
	{
		const Point3f& p0 = mesh.point(v);
		for (auto h = mesh.voh_begin(v); h != mesh.voh_end(v); ++h){
			const Point3f& p1 = mesh.point(mesh.to_vertex_handle(*h));
			const Point3f& p2 =
				mesh.point(mesh.to_vertex_handle(mesh.ccw_rotated_halfedge_handle(*h)));
			const Point3f p01 = (p1 - p0).normalized();
			const Point3f p02 = (p2 - p0).normalized();
			ftype cos_angle = fixCos(dot(p01, p02));
			angles += acos(cos_angle);
		}
	}
	return angles;
}

MeshDoctor::ftype MeshDoctor::calcVertexCurvature(TriMesh& mesh, const VertexHandle& v, CurvatureType ct) {
	using namespace OpenMesh;
	ftype ctval = 0;
	const ftype area = calcVoronoiArea(mesh, v);
	switch (ct)
	{
	case CurvatureType::idMean:
		ctval= ftype(0.5f) * norm(calcLaplace(mesh, v));
		break;
	case CurvatureType::idGauss:
		ctval = (2.0f * M_PI - calcAngleSum(mesh, v)) / area;
		break;
	case CurvatureType::idMax:
	{
		ftype mean = ftype(0.5f) * norm(calcLaplace(mesh, v));
		ftype gauss= (2.0f * M_PI - calcAngleSum(mesh, v)) / area;
	    ftype s = sqrt(std::max(ftype(0.0f), mean * mean - gauss));
		ctval = mean + s;
	}
		break;
	case CurvatureType::idMin:
	{
		ftype mean = ftype(0.5f) * norm(calcLaplace(mesh, v));
		ftype gauss = (2.0f * M_PI - calcAngleSum(mesh, v)) / area;
		ftype s = sqrt(std::max(ftype(0.0f), mean * mean - gauss));
		ctval = mean - s;
	}
		break;
	}
	return ctval;
}

void MeshDoctor::calcVertexCurvature(TriMesh& mesh, const VertexHandle& v,
	ftype& minvc, ftype& maxvc, ftype& gaussvc, ftype& meanvc)
{
	using namespace OpenMesh;
	ftype ctval = 0;
	const ftype area = calcVoronoiArea(mesh, v);
	meanvc = ftype(0.5f) * norm(calcLaplace(mesh, v));
	gaussvc = (2.0f * M_PI - calcAngleSum(mesh, v)) / area;
	ftype s = sqrt(std::max(ftype(0.0f), meanvc * meanvc - gaussvc));
	minvc = meanvc - s;
	maxvc = meanvc + s;
}

MeshDoctor::ftype MeshDoctor::calcAreaRatio(const EdgeHandle& eh)
{
	auto h0 = mesh.halfedge_handle(eh, 0);
	auto h1 = mesh.halfedge_handle(eh, 1);
	ftype a0 = FLT_EPSILON, a1 = FLT_EPSILON;
	if (!mesh.is_boundary(h0)) {
		a0 = calcFaceArea(mesh.face_handle(h0));
	}
	if (!mesh.is_boundary(h1)) {
		a1 = calcFaceArea(mesh.face_handle(h1));
	}
	return min(a0, a1) / max(a0, a1);
}

void MeshDoctor::updateSmoothVertex(int n)
{
	if (!SmoothVertex.is_valid())
		mesh.add_property(SmoothVertex);
	int nvert =(int) mesh.n_vertices();
	if (nvert == 0)
		return;
	int nface = (int)mesh.n_faces();
	auto points = mesh.points();
	vector<ftype> area;
	area.resize(nface);
#pragma omp parallel for
	for (int i = 0; i < nface; ++i) {
		FaceHandle fh = FaceHandle(i);
		auto fv = mesh.fv_begin(fh);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		ftype a = (points[o] - points[t]).length();
		ftype b = (points[o] - points[h]).length();
		ftype c = (points[h] - points[t]).length();
		ftype p = (a + b + c) / 2;
		ftype ar = static_cast<ftype>(sqrt(p * (p - a) * (p - b) * (p - c)));
		area[i] = ar;
	}

	if (!mesh.has_face_normals())
		mesh.request_face_normals();
	mesh.update_face_normals();
#pragma omp parallel for
	for (int i = 0; i < nvert; ++i) {
		VertexHandle vh(i);
		vector<FaceHandle> fhs;
		selectVFace(vh, n, fhs, true);
		Matrix mat(3, 3);
		mat.setZero();
		ftype ar = ftype(0);
		for (auto& fh : fhs) {
			ar += area[fh.idx()];
		}
		ar /= static_cast<ftype>(area.size());
		for (auto& fh : fhs) {
			Point3f norm = mesh.calc_face_normal(fh);
			ftype cu = area[fh.idx()] / ar*exp(-(mesh.calc_face_centroid(fh) - points[i]).norm());
			mat(0, 0) += cu*norm[0] * norm[0];
			mat(0, 1) += cu*norm[0] * norm[1];
			mat(0, 2) += cu*norm[0] * norm[2];
			mat(1, 1) += cu*norm[1] * norm[1];
			mat(1, 2) += cu*norm[1] * norm[2];
			mat(2, 2) += cu*norm[2] * norm[2];
		}
		mat(1, 0) = mat(0, 1);
		mat(2, 1) = mat(1, 2);
		mat(2, 0) = mat(0, 2);
		Eigen::SelfAdjointEigenSolver<Matrix> solver(mat);
		if (solver.info() != Eigen::Success) {
			std::cout << "solver is faiure!" << std::endl;
			mesh.property(SmoothVertex, vh)=Point3f(0,0,0);
			continue;
		}
		auto eigen_val = solver.eigenvalues();
		Point3f p3f;
		ftype maxv = eigen_val(0) > eigen_val(1) ? eigen_val(0) : eigen_val(1);
		maxv = maxv > eigen_val(2) ? maxv : eigen_val(2);
		ftype minv = eigen_val(0) < eigen_val(1) ? eigen_val(0) : eigen_val(1);
		minv = minv < eigen_val(2) ? minv : eigen_val(2);
		ftype midv = eigen_val(0)+ eigen_val(1)+ eigen_val(2)-minv-maxv;
		ftype sum = maxv + midv + minv;
		p3f[0] = maxv;
		p3f[1] = midv;
		p3f[2] = minv;
		mesh.property(SmoothVertex, vh) = p3f;
	}

	ftype maxv=-2E-10F,minv=FLT_MAX;
	for (int i = 0; i < nvert; ++i) {
		auto p = mesh.property(SmoothVertex, VertexHandle(i));
		maxv = maxv > p[0] ? maxv : p[0];
		minv = minv < p[2] ? minv : p[2];
		minv = minv < p[1] ? minv : p[1];
		minv = minv < p[0] ? minv : p[0];
	}
	for (int i = 0; i < nvert; ++i) {
		auto& p = mesh.property(SmoothVertex, VertexHandle(i));
		p[0] = (p[0]-minv)/(maxv-minv);
		p[1] = (p[1]-minv)/(maxv-minv);
		p[2] = (p[2]-minv)/(maxv-minv);
	}
}

MeshDoctor::Point3f MeshDoctor::calcSmoothVertex(const VertexHandle& vh, int n)
{
	vector<FaceHandle> fhs;
	auto points = mesh.points();
	selectVFace(vh, n, fhs, true);
	Matrix mat(3, 3);
	mat.setZero();
	ftype ar = ftype(0);
	vector<ftype> area;
	for (auto& fh : fhs) {
		ftype arr = calcFaceArea(fh);
		ar += arr;
		area.push_back(arr);
	}
	int nfhs = (int)fhs.size();
	for (int i = 0; i < nfhs;++i) {
		auto& fh = fhs[i];
		Point3f norm = mesh.calc_face_normal(fh);
		ftype cu = area[fh.idx()] / (ar*(mesh.calc_face_centroid(fh) - points[vh.idx()]).norm());
		mat(0, 0) += cu * norm[0] * norm[0];
		mat(0, 1) += cu * norm[0] * norm[1];
		mat(0, 2) += cu * norm[0] * norm[2];
		mat(1, 1) += cu * norm[1] * norm[1];
		mat(1, 2) += cu * norm[1] * norm[2];
		mat(2, 2) += cu * norm[2] * norm[2];
	}
	mat(1, 0) = mat(0, 1);
	mat(2, 1) = mat(1, 2);
	mat(2, 0) = mat(0, 2);
	Eigen::SelfAdjointEigenSolver<Matrix> solver(mat);
	if (solver.info() != Eigen::Success) {
		std::cout << "solver is faiure!" << std::endl;
		return Point3f(0,0,0);
	}
	auto eigen_val = solver.eigenvalues();
	Point3f p3f;
	ftype maxv = eigen_val(0) > eigen_val(1) ? eigen_val(0) : eigen_val(1);
	maxv = maxv > eigen_val(1) ? maxv : eigen_val(2);
	ftype minv = eigen_val(0) < eigen_val(1) ? eigen_val(0) : eigen_val(1);
	minv = minv < eigen_val(1) ? minv : eigen_val(2);
	ftype midv = eigen_val(0) + eigen_val(1) + eigen_val(2) - minv - maxv;
	p3f[0] = maxv - midv;
	p3f[1] = midv - minv;
	p3f[2] = minv;
	return p3f;
}

void MeshDoctor::selectVerts(Triangle_mesh& mesh,
	const std::vector<int>& that, int n, std::vector<int>&spts, bool issave)
{
	using namespace std;
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
		if (issave == false)
			spts.insert(spts.end(), loops.begin(), loops.end());
		else
			spts.insert(spts.end(), curloop.begin(), curloop.end());
		loops.swap(curloop);
	}
	spts.insert(spts.end(), loops.begin(), loops.end());
}
void MeshDoctor::selectFace(const FaceHandle& fh, int n, std::vector<FaceHandle>& fhs, bool issave)
{
	fhs.clear();
	if (!mesh.is_valid_handle(fh))
		return;
	int nface = static_cast<int>(mesh.n_faces());
	vector<bool> iscall(nface, false);
	vector<FaceHandle> loop;
	loop.push_back(fh);
	fhs.reserve(nface);
	if (issave)
		fhs.push_back(fh);
	int c = 0;
	iscall[fh.idx()] = true;
	while (c < n) {
		vector<FaceHandle> curloop;
		curloop.reserve(50);
		for (auto& vh : loop) {
			for (auto ff = mesh.ff_begin(vh); ff != mesh.ff_end(vh); ++ff) {
				if (iscall[ff->idx()]==false) {
					curloop.push_back(*ff);
					iscall[ff->idx()] = true;
				}
			}
		}
		fhs.insert(fhs.end(),curloop.begin(), curloop.end());
		loop.swap(curloop);
		++c;
	}
}


void MeshDoctor::deleteIrrVerts( ftype esmallcomp)
{
	if (!mesh.has_edge_status())
		mesh.request_edge_status();
	if (!mesh.has_vertex_status())
		mesh.request_vertex_status();
	if (!mesh.has_face_status())
		mesh.request_face_status();
	mesh.delete_isolated_vertices(); //删除固定点
	while (true)
	{
		bool isok = false;
		for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv)
			if (!mesh.is_manifold(*iv))
				isok = true, mesh.delete_vertex(*iv);
		if (!isok)
			break;
	}
	mesh.garbage_collection();
	//标记小组件
	string compent = "compent";
	mesh.add_property(CompIndex, compent);
	vector<bool> iscall(mesh.n_vertices(), false);
	int id = 0;
	for (auto cm = mesh.vertices_begin(); cm != mesh.vertices_end(); ++cm) {
		if (iscall[cm->idx()])
			continue;
		iscall[cm->idx()] = true;
		mesh.property(CompIndex, *cm) = id;
		stack<int> stck;
		stck.push(cm->idx());
		do
		{
			int tp = stck.top();
			stck.pop();
			auto cv = mesh.vertex_handle(tp);
			for (auto iv = mesh.vv_begin(cv); iv != mesh.vv_end(cv); ++iv)
			{
				if (iscall[iv->idx()] == false)
				{
					stck.push(iv->idx());
					mesh.property(CompIndex, *iv) = id;
					iscall[iv->idx()] = true;
				}
			}
		} while (stck.size() > 0);
		++id;
	}
	ftype nvert = (ftype)mesh.n_vertices();
	vector<int> ids(id, 0);
	for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
		++ids[mesh.property(CompIndex, *iv)];
	}
	vector<int> delid;
	for (int i = 0; i < id; ++i) {
		if ((ftype)ids[i] / nvert < esmallcomp)
			delid.push_back(i);
	}
	for (auto iddel : delid) {
		for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
			if (mesh.property(CompIndex, *iv) == iddel)
				mesh.delete_vertex(*iv);
		}
	}
	mesh.garbage_collection();
	mesh.release_vertex_status();
	mesh.release_face_status();
	mesh.release_edge_status();
	mesh.remove_property(CompIndex);
}


void MeshDoctor::deleteIrrVerts()
{
	if (!mesh.has_edge_status())
		mesh.request_edge_status();
	if (!mesh.has_vertex_status())
		mesh.request_vertex_status();
	if (!mesh.has_face_status())
		mesh.request_face_status();
	mesh.delete_isolated_vertices(); //删除固定点
	while (true)
	{
		bool isok = false;
		for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv)
			if (!mesh.is_manifold(*iv))
				isok = true, mesh.delete_vertex(*iv);
		if (!isok)
			break;
	}
	mesh.garbage_collection();
	//标记小组件
	string compent = "compent";
	mesh.add_property(CompIndex, compent);
	vector<bool> iscall(mesh.n_vertices(), false);
	int id = 0;
	for (auto cm = mesh.vertices_begin(); cm != mesh.vertices_end(); ++cm) {
		if (iscall[cm->idx()])
			continue;
		iscall[cm->idx()] = true;
		mesh.property(CompIndex, *cm) = id;
		stack<int> stck;
		stck.push(cm->idx());
		do
		{
			int tp = stck.top();
			stck.pop();
			auto cv = mesh.vertex_handle(tp);
			for (auto iv = mesh.vv_begin(cv); iv != mesh.vv_end(cv); ++iv)
			{
				if (iscall[iv->idx()] == false)
				{
					stck.push(iv->idx());
					mesh.property(CompIndex, *iv) = id;
					iscall[iv->idx()] = true;
				}
			}
		} while (stck.size() > 0);
		++id;
	}
	ftype nvert = (ftype)mesh.n_vertices();
	vector<int> ids(id, 0);
	for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
		++ids[mesh.property(CompIndex, *iv)];
	}
	int maxid = 0;
	size_t coucc = 0;
	for (int i = 0; i < id; ++i) {
		if (coucc < ids[i]) {
			coucc = ids[i];
			maxid = i;
		}
	}
	vector<int> delid;
	for (int i = 0; i < id; ++i) {
		if (i!=maxid)
			delid.push_back(i);
	}
	for (auto iddel : delid) {
		for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
			if (mesh.property(CompIndex, *iv) == iddel)
				mesh.delete_vertex(*iv);
		}
	}
	mesh.garbage_collection();
	mesh.release_vertex_status();
	mesh.release_face_status();
	mesh.release_edge_status();
	mesh.remove_property(CompIndex);
}


/*
* @brief 网格去重
* @detail
* @date 2020.3.24
*/
void MeshDoctor::repeatMesh( ftype _eps)
{
	Triangle_mesh rmesh;
	using namespace std;
	crepeator<Point3f>::eps_ = _eps;
	map<Point3f, int, crepeator<Point3f>> repoter;
	rmesh.clean();
	rmesh.reserve(mesh.n_vertices(), mesh.n_edges(), mesh.n_faces());
	int nvert = (int)mesh.n_vertices();
	auto points = mesh.points();
	int ncount = -1;
	vector<VertexHandle> handles;
	handles.reserve(nvert);
	for (int i = 0; i < nvert; ++i) {
		if (repoter.count(points[i]) == 0) {
			repoter.insert(make_pair(points[i], ++ncount));
			handles.push_back(rmesh.add_vertex(points[i]));
		}
	}
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		auto fv = mesh.fv_begin(*f);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		auto oh = handles[repoter[points[o]]];
		auto th = handles[repoter[points[t]]];
		auto hh = handles[repoter[points[h]]];
		rmesh.add_face(oh, th, hh);
	}
	mesh = rmesh;
}

/*
* @brief 网格去重(保存纹理坐标)
* @detail
* @date 2020.3.24
*/
void MeshDoctor::repeatMeshTextCoord( ftype _eps)
{
	Triangle_mesh rmesh;
	using namespace std;
	crepeator<Point3f>::eps_ = _eps;
	map<Point3f, int, crepeator<Point3f>> repoter;
	rmesh.clean();
	rmesh.reserve(mesh.n_vertices(), mesh.n_edges(), mesh.n_faces());
	int nvert = (int)mesh.n_vertices();
	auto points = mesh.points();
	int ncount = -1;
	vector<VertexHandle> handles;
	handles.reserve(nvert);
	for (int i = 0; i < nvert; ++i) {
		if (repoter.count(points[i]) == 0) {
			repoter.insert(make_pair(points[i], ++ncount));
			handles.push_back(rmesh.add_vertex(points[i]));
		}
	}
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		auto fv = mesh.fv_begin(*f);
		int o = fv->idx(); ++fv;
		int t = fv->idx(); ++fv;
		int h = fv->idx();
		auto oh = handles[repoter[points[o]]];
		auto th = handles[repoter[points[t]]];
		auto hh = handles[repoter[points[h]]];
		rmesh.add_face(oh, th, hh);
	}
	if (mesh.has_vertex_texcoords2D()) {
		if (!rmesh.has_vertex_texcoords2D())
			rmesh.request_vertex_texcoords2D();
		for (int i = 0; i < nvert; ++i) {
			rmesh.set_texcoord2D(handles[repoter[points[i]]], mesh.texcoord2D(VertexHandle(i)));
		}
	}
	mesh = rmesh;
}

void MeshDoctor::nLocalRemoveAspectRatio(ftype fc)
{
	vector<bool> iscallf(mesh.n_faces(), false);
	vector<bool> iscallv(mesh.n_vertices(), true);
	vector<FaceHandle> fhs;
	vector<VertexHandle>remeshfhs;
	fhs.reserve(mesh.n_faces());
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		if (mesh.property(AspectRatio, *f) > fc) {
			iscallf[f->idx()] = true;
		}
	}
	selectFVace(fhs, 3, remeshfhs, true);
	if(!mesh.has_vertex_status())
		mesh.request_vertex_status();
	for (auto& v : remeshfhs)
		iscallv[v.idx()] = false;
	//for (auto v = mesh.vertices_begin(); v != mesh.vertices_end(); ++v) {
	//	mesh.status(*v).set_locked(iscallv[v->idx()]);
	//}
	OpenMesh::Smoother::JacobiLaplaceSmootherT< Triangle_mesh > smoother(mesh);
	smoother.initialize(OpenMesh::Smoother::SmootherT< Triangle_mesh >::
		Tangential_and_Normal,
		OpenMesh::Smoother::SmootherT< Triangle_mesh >::C1);
	smoother.smooth(40);
	mesh.release_vertex_status();
}

void MeshDoctor::checkGenius()
{
	using namespace GeniusCore;
	if (!Genius.is_valid())
		mesh.add_property(Genius,"Genius");
	CMeanValueParameterizer cparam(mesh);
	cparam.Parameterize();
	Triangle_mesh uvmesh;
	uvmesh.reserve(mesh.n_vertices(), mesh.n_edges(), mesh.n_faces());
	for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
		auto uv = mesh.texcoord2D(*iv);
		uvmesh.add_vertex(Point3f(uv[0], uv[1], 0));
	}
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		auto fv = mesh.fv_begin(*f);
		auto o = *fv; ++fv;
		auto t = *fv; ++fv;
		auto h = *fv; 
		uvmesh.add_face(o, t, h);
	}
	mesh.release_vertex_texcoords2D();
	uvmesh.request_vertex_normals();
	uvmesh.request_face_normals();
	uvmesh.update_face_normals();
	uvmesh.update_vertex_normals();
	Point3f p3f = { 0,0,0 };
	for (auto iv = uvmesh.vertices_begin(); iv != uvmesh.vertices_end(); ++iv) {
		p3f += uvmesh.calc_vertex_normal(*iv);
	}
	p3f /= static_cast<ftype>(uvmesh.n_vertices());
	for (auto iv = uvmesh.vertices_begin(); iv != uvmesh.vertices_end(); ++iv) {
		auto p=uvmesh.calc_vertex_normal(*iv);
		if ((p | p3f) > 0)
			mesh.property(Genius, *iv) = false;
		else
			mesh.property(Genius, *iv) = true;
	}
	vector<vector<int>> ids;
	extraBoundary(&mesh, ids);
	for (auto& id : ids) {
		for (auto& icd : id) {
			mesh.property(Genius, VertexHandle(icd)) = false;
		}
	}
	uvmesh.release_face_normals();
	uvmesh.release_vertex_status();
}

void MeshDoctor::selectVFace(const VertexHandle& fh, int n, std::vector<FaceHandle>& fhs, bool issave)
{
	fhs.clear();
	if (!mesh.is_valid_handle(fh))
		return;
	int nface = static_cast<int>(mesh.n_faces());
	vector<bool> iscall(nface, false);
	vector<FaceHandle> loop;
	for (auto vf = mesh.vf_begin(fh); vf != mesh.vf_end(fh); ++vf)
		loop.push_back(*vf);
	fhs.reserve(nface);
	if (issave)
		fhs.insert(fhs.end(), loop.begin(), loop.end());
	int c = 0;
	iscall[fh.idx()] = true;
	while (c < n) {
		vector<FaceHandle> curloop;
		curloop.reserve(50);
		for (auto& vh : loop) {
			for (auto ff = mesh.ff_begin(vh); ff != mesh.ff_end(vh); ++ff) {
				if (iscall[ff->idx()] == false) {
					curloop.push_back(*ff);
					iscall[ff->idx()] = true;
				}
			}
		}
		fhs.insert(fhs.end(), curloop.begin(), curloop.end());
		loop.swap(curloop);
		++c;
	}
}

void MeshDoctor::selectFVace(const std::vector<FaceHandle>& fh, int n, std::vector<VertexHandle>& fhs, bool issave)
{
	fhs.clear();
	int nface = static_cast<int>(mesh.n_vertices());
	vector<bool> iscall(nface, false);
	vector<VertexHandle> loop;
	for (auto&f : fh) {
		if (mesh.is_valid_handle(f)) {
			for (auto fv = mesh.fv_begin(f); fv != mesh.fv_end(f); ++fv) {
				if (iscall[fv->idx()] == false) {
					iscall[fv->idx()] = true;
					loop.push_back(*fv);
				}
			}
		}
	}
	fhs.reserve(nface);
	if (issave)
		fhs.insert(fhs.end(),loop.begin(),loop.end());
	int c = 0;
	while (c < n) {
		vector<VertexHandle> curloop;
		curloop.reserve(50);
		for (auto& vh : loop) {
			for (auto ff = mesh.vv_begin(vh); ff != mesh.vv_end(vh); ++ff) {
				if (iscall[ff->idx()] == false) {
					curloop.push_back(*ff);
					iscall[ff->idx()] = true;
				}
			}
		}
		fhs.insert(fhs.end(), curloop.begin(), curloop.end());
		loop.swap(curloop);
		++c;
	}
}

void MeshDoctor::selectVerts(const VertexHandle& fh, int n, std::vector<VertexHandle>& fhs, bool issave)
{
	fhs.clear();
	if (!mesh.is_valid_handle(fh))
		return;
	int nface = static_cast<int>(mesh.n_faces());
	vector<bool> iscall(nface, false);
	vector<VertexHandle> loop;
	loop.push_back(fh);
	fhs.reserve(nface);
	if (issave)
		fhs.push_back(fh);
	int c = 0;
	iscall[fh.idx()] = true;
	while (c < n) {
		vector<VertexHandle> curloop;
		curloop.reserve(50);
		for (auto& vh : loop) {
			for (auto ff = mesh.vv_begin(vh); ff != mesh.vv_end(vh); ++ff) {
				if (iscall[ff->idx()] == false) {
					curloop.push_back(*ff);
					iscall[ff->idx()] = true;
				}
			}
		}
		fhs.insert(fhs.end(), curloop.begin(), curloop.end());
		loop.swap(curloop);
		++c;
	}
}

void MeshDoctor::selectFace(const std::vector<FaceHandle>& fh, 
	int n, std::vector<FaceHandle>& fhs, bool issave)
{
	fhs.clear();
	int nface = static_cast<int>(mesh.n_faces());
	vector<bool> iscall(nface, false);
	vector<FaceHandle> loop;
	loop=fh;
	fhs.reserve(nface);
	if (issave)
		fhs.insert(fhs.end(),fh.begin(),fh.end());
	int c = 0;
	for(auto& f:fh)
		iscall[f.idx()] = true;
	while (c < n) {
		vector<FaceHandle> curloop;
		curloop.reserve(50);
		for (auto& vh : loop) {
			for (auto ff = mesh.ff_begin(vh); ff != mesh.ff_end(vh); ++ff) {
				if (iscall[ff->idx()] == false) {
					curloop.push_back(*ff);
					iscall[ff->idx()] = true;
				}
			}
		}
		fhs.insert(fhs.end(), curloop.begin(), curloop.end());
		loop.swap(curloop);
		++c;
	}
}


void MeshDoctor::nRemeshKeepBoundry(Triangle_mesh& mesh,ftype avedge)
{
	
}


void MeshDoctor::extraBoundary(Triangle_mesh* mesh, std::vector<std::vector<int>>& caids)
{
	using namespace std;
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



bool MeshDoctor::smoothBoundaryDelVertex()
{
	using Point = OpenMesh::Vec3f;
	auto clearHoleFunc = [](Triangle_mesh& mesh)
	{
		bool isdegnerate = true;
		while (isdegnerate)
		{
			isdegnerate = false;
			auto v_end = mesh.vertices_end();
			Point point;
			for (auto v_it = mesh.vertices_sbegin(); v_it != v_end; ++v_it)
			{
				VertexHandle vh = *v_it;
				if (!mesh.is_boundary(*v_it))
					continue;
				point.vectorize(0.0f);
				int ivalence = 0;
				Point vhp = mesh.point(vh);
				std::vector<VertexHandle>  vcircle;
				for (auto vv_it = mesh.cvv_iter(vh); vv_it.is_valid(); ++vv_it)
				{
					if (mesh.is_boundary(*vv_it))
					{
						point += 0.5f*(vhp + mesh.point(*vv_it));
						ivalence++;
						vcircle.push_back(*vv_it);
					}
				}
				if (ivalence > 2)
				{
					mesh.delete_vertex(vh);
					isdegnerate = true;
				}
				else if (2 == ivalence)
				{
					HalfedgeHandle starthh0 = mesh.find_halfedge(vh, vcircle[0]);
					HalfedgeHandle starthh1 = mesh.opposite_halfedge_handle(starthh0);
					FaceHandle fh0;
					if (mesh.is_boundary(starthh0))
					{
						fh0 = mesh.face_handle(starthh1);
					}
					else
					{
						fh0 = mesh.face_handle(starthh0);
					}
					//
					HalfedgeHandle endh0 = mesh.find_halfedge(vh, vcircle[1]);
					HalfedgeHandle endh1 = mesh.opposite_halfedge_handle(endh0);
					FaceHandle fh1;
					if (mesh.is_boundary(endh0))
					{
						fh1 = mesh.face_handle(endh1);
					}
					else
					{
						fh1 = mesh.face_handle(endh0);
					}
					const Point& vp0 = mesh.point(vcircle[0]);
					const Point& vp1 = mesh.point(vcircle[1]);
					Point e0 = (vhp - vp0).normalized();
					Point e1 = (vhp - vp1).normalized();
					float dote01 = (e0 | e1);
					if (dote01 > -0.3090f)
					{
						point *= (1.0f / (1.0f*ivalence));
						mesh.set_point(vh, point);
					}
				}
				else
				{
					isdegnerate = true;
					mesh.delete_vertex(vh);
				}
			}
		}
		mesh.garbage_collection();
	};
	//
	if (!mesh.has_edge_status())
		mesh.request_edge_status();
	if (!mesh.has_face_status())
		mesh.request_face_status();
	if (!mesh.has_vertex_status())
		mesh.request_vertex_status();
	//删除孤立顶点
	mesh.delete_isolated_vertices();
	mesh.garbage_collection();
	HalfedgeHandle bhh;
	auto fh_end = mesh.halfedges_end();
	for (auto fh_it = mesh.halfedges_begin(); fh_it != fh_end; ++fh_it)
	{
		if (mesh.is_boundary(*fh_it))
		{
			bhh = *fh_it;
			break;
		}
	}
	float maxangle = 0.523599f;//30
	if (!mesh.is_boundary(bhh))
	{
		return false;
	}
	HalfedgeHandle current = bhh;
	HalfedgeHandle breakhh = bhh;

	while (true)
	{
		bool newface = false;
		if (!mesh.is_boundary(current))
		{
			break;
		}
		HalfedgeHandle nexthh = mesh.next_halfedge_handle(current);
		HalfedgeHandle prehh = mesh.prev_halfedge_handle(current);
		VertexHandle vh0 = mesh.to_vertex_handle(current);
		VertexHandle vh1 = mesh.to_vertex_handle(prehh);
		VertexHandle vh2 = mesh.to_vertex_handle(nexthh);
		const Point& p0 = mesh.point(vh0);
		const Point& p1 = mesh.point(vh1);
		const Point& p2 = mesh.point(vh2);
		Point e1 = p0 - p1;
		Point e2 = p2 - p1;
		float v0angle = (e1 | e2) / (e1.length()*e2.length());
		if (v0angle < -1.0f)
			v0angle = -1.0f;
		if (v0angle > 1.0f)
			v0angle = 1.0f;
		v0angle = acos(v0angle);
		if (v0angle < 0.3f)
		{
			current = nexthh;
			if ((current == breakhh) && (!newface))
				break;
			continue;
		}
		FaceHandle currentfh = mesh.face_handle(mesh.opposite_halfedge_handle(current));
		FaceHandle nextfh = mesh.face_handle(mesh.opposite_halfedge_handle(nexthh));
		Point nor0 = (e1%e2).normalized();
		Point nor1 = mesh.calc_face_normal(currentfh);
		Point nor2 = mesh.calc_face_normal(nextfh);
		float v1dot = (nor0 | nor1) / (nor0.length()*nor1.length());
		float v2dot = (nor0 | nor2) / (nor0.length()*nor2.length());
		if (v1dot < -1.0f)
			v1dot = -1.0f;
		if (v2dot > 1.0f)
			v2dot = 1.0f;
		float v1angle = acos(v1dot);
		float v2angle = acos(v2dot);
		if ((v1angle > maxangle) || (v2angle > maxangle))
		{
			current = nexthh;
			if ((current == breakhh) && (!newface))
				break;
			continue;
		}
		//
		FaceHandle fh = mesh.add_face(vh1, vh0, vh2);
		if (!fh.is_valid())
			break;
		current = prehh;
		breakhh = prehh;
		newface = true;
	}
	clearHoleFunc(mesh);
	return true;
}


void MeshDoctor::debugAngle(const std::string& filename)
{
	updateAngle();
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		mesh.set_color(*f, BLUE);
		if (mesh.property(MinAngle, *f) < AG_10 || mesh.property(MaxAngle, *f) > AG_145)
			mesh.set_color(*f, RED);
		else
			mesh.set_color(*f, BLUE);
	}
	OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::FaceColor);
}

void MeshDoctor::debugAspectRatio(const std::string& filename) 
{
	updateAspectRatio();
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		mesh.set_color(*f, BLUE);
		if (mesh.property(AspectRatio, *f) >=4)
			mesh.set_color(*f, RED);
		else
			mesh.set_color(*f, BLUE);
	}
	OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::FaceColor);
}

void MeshDoctor::debugSkewnessEquiangle(const std::string& filename) {
	updateSkewnessEquiangle();
	for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) {
		mesh.set_color(*f, BLUE);
		if (mesh.property(SkewnessEquiangle, *f) >=0.8f)
			mesh.set_color(*f, RED);
		else
			mesh.set_color(*f, BLUE);
	}
	OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::FaceColor);
}

void MeshDoctor::debugGenus(const std::string& filename)
{
	//repeatMesh();
	deleteIrrVerts();
	checkGenius();
	for (auto iv = mesh.vertices_begin(); iv != mesh.vertices_end(); ++iv) {
		mesh.set_color(*iv, BLUE);
		if (mesh.property(MeshDoctor::Genius, *iv)) {
			mesh.set_color(*iv, RED);
		}
		else {
			mesh.set_color(*iv, BLUE);
		}
	}
	OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::VertexColor);
}

void MeshDoctor::debugFaceColor(const std::string& filename, int faceIdx)
{
	if (!mesh.has_face_colors())
		mesh.request_face_colors();
	/*for (auto f = mesh.faces_begin(); f != mesh.faces_end(); ++f) 
	{
		int faceIdx = f->idx();
		mesh.set_color(*f, { 255,0,0 });
	}*/
	mesh.set_color(FaceHandle(faceIdx), { 0,0,255 });
	//输出ply格式
	OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::FaceColor);
}

void MeshDoctor::debugFaceColor(const std::string& filename, std::map<int, int>& map)
{
	//if (!mesh.has_face_colors())
	//	mesh.request_face_colors();
	//std::map<int, int>::const_iterator it = map.begin();
	//for (; it != map.end(); ++it)
	//{
	//	if (it->second == 15) {
	//		mesh.set_color(FaceHandle(it->first), { 255,0,0 });
	//	}
	//}
	//mesh.set_color(FaceHandle(2782), { 0,0,255 });
	////输出ply格式
	//OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::FaceColor);

	if (!mesh.has_face_colors())
		mesh.request_face_colors();
	std::map<int, int>::const_iterator it = map.begin();
	for (; it != map.end(); ++it)
	{
		if (it->second == 15) {
			mesh.set_color(FaceHandle(it->first), { 255,0,0 });
		}
	}
	mesh.set_color(FaceHandle(2879), { 0,0,255 });
	//输出ply格式
	OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::FaceColor);



}


namespace GeniusCore
{


	CMeanValueParameterizer::CMeanValueParameterizer(Triangle_mesh& mesh)
		:mTriMesh(mesh), mPropertyName("UVCoord")
	{
		if (!mTriMesh.has_vertex_texcoords2D())
			mTriMesh.request_vertex_texcoords2D();
		mesh.add_property(mUVProperty, mPropertyName);
	}

	CMeanValueParameterizer::~CMeanValueParameterizer()
	{
		mTriMesh.remove_property(mUVProperty);
	}

	std::pair<HalfedgeHandle, ftype> FindLongestBoundary(const Triangle_mesh& mesh)
	{
		std::vector<bool>  halfedgestatus(mesh.n_halfedges(), false);
		auto h_end = mesh.halfedges_end();
		std::vector<ftype>          eathboundlength;
		std::vector<HalfedgeHandle>  eathboundhalfhandle;
		int currentboundsize = 0;
		for (auto h_it = mesh.halfedges_begin(); h_it != h_end; ++h_it)
		{
			if (halfedgestatus[h_it->idx()])
				continue;
			if (mesh.is_boundary(*h_it))
			{
				ftype hhlen = mesh.calc_edge_length(*h_it);
				halfedgestatus[h_it->idx()] = true;
				HalfedgeHandle currenthh = *h_it;
				eathboundlength.push_back(0.0);
				eathboundhalfhandle.push_back(*h_it);
				do
				{
					eathboundlength[currentboundsize] += hhlen;
					currenthh = mesh.next_halfedge_handle(currenthh);
					hhlen = mesh.calc_edge_length(currenthh);
					halfedgestatus[currenthh.idx()] = true;
				} while (currenthh != (*h_it));
				++currentboundsize;
			}
		}
		if (0 == eathboundlength.size())
			return std::make_pair(HalfedgeHandle(), 0.0);
		ftype firstsize = eathboundlength[0];
		size_t maxnumber = 0;
		for (size_t i = 0; i < eathboundlength.size(); ++i)
		{
			if (firstsize < eathboundlength[i])
			{
				firstsize = eathboundlength[i];
				maxnumber = i;
			}
		}

		return std::make_pair(eathboundhalfhandle[maxnumber], eathboundlength[maxnumber]);
	}

	bool CMeanValueParameterizer::Parameterize()
	{
		bool haveproperty = mTriMesh.get_property_handle(mUVProperty, mPropertyName);
		if (!haveproperty)
			return false;
		mTriMesh.add_property(mVertStatus, "ParaStatus");
		auto v_end = mTriMesh.vertices_end();
		for (auto v_it = mTriMesh.vertices_begin(); v_it != v_end; ++v_it)
			mTriMesh.property(mVertStatus, *v_it) = false;
		std::pair<HalfedgeHandle, ftype> boundaryhh = FindLongestBoundary(mTriMesh);
		if (!BorderParameterize(boundaryhh.first, boundaryhh.second))
		{
			mTriMesh.remove_property(mVertStatus);
			return false;
		}
		int nbVertices = static_cast<int>(mTriMesh.n_vertices());
		EigenMatrix A(nbVertices, nbVertices);
		EigenVector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);
		Xu.setZero();
		Xv.setZero();
		Bu.setZero();
		Bv.setZero();
		std::vector<Triplet> Atriplets;
		Atriplets.reserve(nbVertices);
		InitializeSystemBorder(A, Atriplets, Bu, Bv, boundaryhh.first);
		if (!SetInnerVertexRelations(A, Atriplets, Bu, Bv))
			return false;
		A.setFromTriplets(Atriplets.begin(), Atriplets.end());
		A.makeCompressed();
		Solver linesolver;
		linesolver.compute(A);
		Xu = linesolver.solve(Bu);
		Xv = linesolver.solve(Bv);
		for (auto v_it = mTriMesh.vertices_begin(); v_it != v_end; ++v_it){
			bool vstatus = mTriMesh.property(mVertStatus, *v_it);
			if (vstatus)
				continue;
			int index = v_it->idx();
			ftype u = Xu[index];
			ftype v = Xv[index];
			mTriMesh.property(mUVProperty, *v_it) = PointXYd(u, v);
			mTriMesh.set_texcoord2D(*v_it, PointXYd(u, v));
			mTriMesh.property(mVertStatus, *v_it) = true;
		}
		mTriMesh.remove_property(mVertStatus);
		return false;
	}

	bool CMeanValueParameterizer::BorderParameterize(
		const HalfedgeHandle& hh, 
		const ftype borderlen)
	{
		if (!hh.is_valid())
			return false;
		const ftype tmp = 2 * M_PI / borderlen;
		ftype len = 0.0;
		HalfedgeHandle currenthh = hh;
		do
		{
			VertexHandle vh = mTriMesh.to_vertex_handle(currenthh);
			ftype angle = len * tmp;
			PointXYd uv(1.0*std::cos(-angle), 1.0*std::sin(-angle));
			mTriMesh.property(mUVProperty, vh) = uv;
			mTriMesh.property(mVertStatus, vh) = true;
			len += mTriMesh.calc_edge_length(currenthh);
			currenthh = mTriMesh.next_halfedge_handle(currenthh);
		} while (currenthh != hh);
		return true;

	}

	void CMeanValueParameterizer::InitializeSystemBorder(EigenMatrix& A, std::vector<Triplet>& Atriplet,
		EigenVector& Bu, EigenVector& Bv, HalfedgeHandle& hh)
	{
		HalfedgeHandle currenthh = hh;
		do
		{
			VertexHandle vh = mTriMesh.to_vertex_handle(currenthh);
			int index = vh.idx();
			Atriplet.push_back(Triplet(index, index, 1.0));
			const PointXYd uv = mTriMesh.property(mUVProperty, vh);
			Bu[index] = uv[0];
			Bv[index] = uv[1];
			//
			currenthh = mTriMesh.next_halfedge_handle(currenthh);
		} while (currenthh != hh);
	}

	bool CMeanValueParameterizer::SetInnerVertexRelations(EigenMatrix& A, std::vector<Triplet>& Atriplet,
		EigenVector& Bu, EigenVector& Bv)
	{
		auto v_end = mTriMesh.vertices_end();
		for (auto v_it = mTriMesh.vertices_begin(); v_it != v_end; ++v_it)
		{
			bool vstatus = mTriMesh.property(mVertStatus, *v_it);
			if (vstatus)
				continue;
			int i = v_it->idx();
			ftype w_ii = 0;
			int vertexIndex = 0;
			auto vv_circle = mTriMesh.vv_ccwiter(*v_it);
			for (; vv_circle.is_valid(); ++vv_circle)
			{
				ftype w_ij = ComputeWeight(*v_it, *vv_circle);
				w_ii -= w_ij;
				int j = vv_circle->idx();
				Atriplet.push_back(Triplet(i, j, w_ij));
				vertexIndex++;
			}
			if (vertexIndex < 2)
				return false;
			Atriplet.push_back(Triplet(i, i, w_ii));
		}

		return true;
	}

	ftype CMeanValueParameterizer::ComputeWeight(const VertexHandle& vert, const VertexHandle& around)
	{
		HalfedgeHandle hh = mTriMesh.find_halfedge(vert, around);
		HalfedgeHandle opphh = mTriMesh.opposite_halfedge_handle(hh);
		VertexHandle next = mTriMesh.to_vertex_handle(mTriMesh.next_halfedge_handle(hh));
		VertexHandle pre = mTriMesh.to_vertex_handle(mTriMesh.next_halfedge_handle(opphh));
		//
		const Point& position_v_i = mTriMesh.point(vert);
		const Point& position_v_j = mTriMesh.point(around);
		ftype len = (position_v_i - position_v_j).length();
		const Point& position_v_k = mTriMesh.point(pre);
		const Point& position_v_l = mTriMesh.point(next);

		ftype gamma_ij = ComputeAngleRad(position_v_j, position_v_i, position_v_k);
		ftype delta_ij = ComputeAngleRad(position_v_l, position_v_i, position_v_j);

		ftype weight = 0.0f;
		if (len < 1e-7f)
			return -1.0f;
		weight = static_cast<ftype>((std::tan(0.5*gamma_ij) + std::tan(0.5f*delta_ij)) / len);
		if (weight < 1e-7f)
			return -1.0f;

		return weight;
	}

	ftype CMeanValueParameterizer::ComputeAngleRad(const Point& P, const Point& Q, const Point& R)
	{
		Point u = P - Q;
		Point v = R - Q;
		ftype product = u.length() * v.length();
		if (product == 0)
			return 0.0f;
		ftype dot = u | v;
		ftype cosine = dot / product;
		Point w = u % v;
		ftype abs_sine = w.length() / product;
		if (cosine >= 0)
			return std::asin(FixSine(abs_sine));
		else
			return M_PI - std::asin(FixSine(abs_sine));
		return 0.0f;
	}
}