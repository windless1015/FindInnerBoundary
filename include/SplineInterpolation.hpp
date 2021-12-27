#pragma once
#include <vector>
#include <cmath>
#include <assert.h>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/unsupported/Eigen/Splines>
#include <array>
#include <numeric>
#include <set>
#include <nanoflann.hpp>
#include <fstream>
#pragma warning(push)
#pragma warning(disable:4267) 
#ifndef _MATH_DEFINES_DEFINED
#define M_E        2.71828182845904523536   // e
#define M_PI       3.14159265358979323846   // pi
#define M_PI_2     1.57079632679489661923   // pi/2
#define M_PI_4     0.785398163397448309616  // pi/4
#define M_1_PI     0.318309886183790671538  // 1/pi
#define M_2_PI     0.636619772367581343076  // 2/pi
#endif
/*
 * @brief 曲线的插值和逼近算法
 * @author lyc
 * @date 2020.3.3
 * @detail
 * @modify
 *  --add--2020.2.29 HermiteCubic,Interpolation(Point)
 *  --add--2020.2.30 CubeicCurve,Intelpolation(Points)
 *  --add--2020.3.01 PolynomialCurve,Intelpolation(Points)
 *  --add--2020.3.02 SegmentCurve,Intelpolation(Points,Points+Normal)
 *  --add--2020.3.03 BSpline,Intelpolation(Points,Points+Normal)
 *  --add--2020.3.07 BSpline,Approximation(Points)
 *  --add--2020.3.09 distribution points
			--(Average（均匀分布）,Linear（线性）,exp(指数),Parabola(抛物线),copy(拷贝),Gompertz(S曲线))
 *  --add--2020.3.09 BSpline,Intelpolation-Approximation(闭合曲线功能)
 *  --add--2020.3.11 BSpline,parameterization methods(曲线参数化方法)
 *       --see also Parameterizatio Method on B-Spline Curve [H] H.Haron A.Rehman D.I.S.Adi S.P.Lim T.Saba
 *  --add--2020.3.12 Nurbs,添加nurbs内核的功能
 *  --add--2020.3.12 IO,统一IO接口，方便调试
 *  --fix--2020.3.13 修改B样条计算切向量的bug
 *  --add--2020.3.13 添加B样条计算切向量的调试接口
 *  --add--2020.3.16 添加多段线的反算接口
 *  --add--2020.3.25 添加封闭曲线拟合接口
 */

namespace FussenAlgo
{

	namespace BasicTypeCore
	{
		const static double epslion = (1e-10);
#ifndef M_PI
		const static double  pi = double(3.14159265358979323846264338327950288);
#else
		const static double  pi = double(M_PI);
#endif
	}
	enum class ParamType {
		idAverage, //均匀参数
		idChordLength,//累计弦长参数化
		idCentripetal, //向心参数化
		idFoley, //福利参数化
		idHybrid, //混合参数化
		idExp, //指数参数化
	};
	enum class DerivativeType {
		idSub, //前插值法
		idFmill, //弗米尔求导
		idBessel, //贝赛尔求导
	};

	enum class CurveType {

		idHermiteCubic, //加权hermite曲线，可分段
		idCubeicCurve,  //标准Hermite曲线，可分段

		idPolynomialCurve, //多项式曲线，不可分段
		idSegmentCurve, //分段曲线，不可分段
		idBSpline,    //带约束B样条曲线

		//有待实现
		idNurbsCurve, //非均匀B样条曲线
		idBSplineRational, //有理B样条曲线
		idNurbsCurveRational, //非均匀有理B样条曲线
	};

	template<class TVert, class ftype, int dim = 3>
	class kdtreenode
	{
	public:
		virtual ~kdtreenode() {}
	public:
		size_t kdtree_get_point_count() const { return mpts.size(); }
		ftype kdtree_get_pt(const size_t idx, const size_t idm) const {
			return mpts[idx][idm];
		}
		template <class BBOX>
		bool kdtree_get_bbox(BBOX&) const { return false; }
	public:
		std::vector<TVert> mpts;
		std::vector<ftype> param;
	};



	template <typename T> class array2
	{
	public:
		array2() = default;
		array2(const array2<T> &arr) = default;
		array2 &operator=(const array2 &arr) = default;
		array2(size_t rows, size_t cols, T default_value = T()) { resize(rows, cols, default_value); }
		array2(size_t rows, size_t cols, const std::vector<T> &arr)
			: rows_(rows), cols_(cols), data_(arr)
		{
			if (arr.size() != rows * cols)
			{
				throw std::runtime_error("Dimensions do not match with size of vector");
			}
		}
		void resize(size_t rows, size_t cols, T val = T())
		{
			data_.resize(rows * cols, val);
			rows_ = rows;
			cols_ = cols;
		}
		void clear()
		{
			rows_ = cols_ = 0;
			data_.clear();
		}
		T operator()(size_t row, size_t col) const
		{
			assert(row < rows_ && col < cols_);
			return data_[row * cols_ + col];
		}
		T &operator()(size_t row, size_t col)
		{
			assert(row < rows_ && col < cols_);
			return data_[row * cols_ + col];
		}
		T operator[](size_t idx) const
		{
			assert(idx < data_.size());
			return data_[idx];
		}
		T &operator[](size_t idx)
		{
			assert(idx < data_.size());
			return data_[idx];
		}
		size_t rows() const { return rows_; }
		size_t cols() const { return cols_; }
		size_t size() const { return data_.size(); }

	private:
		size_t rows_, cols_;
		std::vector<T> data_;
	};

	template<class TVert, class ftype, int dim = 3>
	class AlgoUtil
	{
		//using ftype = float;
		//const static int dim = 3;
		//using TVert = std::array<ftype, dim>;
	public: //基础算法
		static ftype length(const TVert& p, const TVert& q) {
			ftype len = static_cast<ftype>(0);
			for (int i = 0; i < dim; ++i)
				len += (p[i] - q[i])*(p[i] - q[i]);
			return static_cast<ftype>(std::sqrt(len));
		}
		static ftype length(const TVert& p) {
			ftype len = static_cast<ftype>(0);
			for (int i = 0; i < dim; ++i)
				len += p[i] * p[i];
			return static_cast<ftype>(std::sqrt(len));
		}
		static ftype length2(const TVert& p, const TVert& q) {
			ftype len = static_cast<ftype>(0);
			for (int i = 0; i < dim; ++i)
				len += (p[i] - q[i])*(p[i] - q[i]);
			return len;
		}
		static ftype length2(const TVert& p) {
			ftype len = static_cast<ftype>(0);
			for (int i = 0; i < dim; ++i)
				len += p[i] * p[i];
			return len;
		}
		static ftype dot(const TVert& p, const TVert& q)
		{
			ftype pq = static_cast<ftype>(0);
			for (int i = 0; i < dim; ++i)
				pq += p[i] * q[i];
			return pq;
		}
		static TVert cross(const TVert& p, const TVert& q)
		{
			static_assert(dim == 3, "only for dim==3");
			TVert pq;
			pq[0] = p[1] * q[2] - p[2] * q[1];
			pq[1] = p[2] * q[0] - p[0] * q[2];
			pq[2] = p[0] * q[1] - p[1] * q[2];
			return pq;
		}
		static TVert sub(const TVert& p, const TVert& q)
		{
			TVert pq;
			for (int i = 0; i < dim; ++i)
				pq[i] = p[i] - q[i];
			return pq;
		}
		static TVert add(const TVert& p, const TVert& q)
		{
			TVert pq;
			for (int i = 0; i < dim; ++i)
				pq[i] = p[i] + q[i];
			return pq;
		}
		static TVert muilt(const ftype& p, const TVert& q)
		{
			TVert pq;
			for (int i = 0; i < dim; ++i)
				pq[i] = p * q[i];
			return pq;
		}
		static TVert divd(const TVert& p, const ftype& q)
		{
			assert(q != 0 && "q can not be zero!");
			TVert pq;
			for (int i = 0; i < dim; ++i)
				pq[i] = p[i] / q;
			return pq;
		}
		static TVert normal(const TVert&p) {
			using namespace BasicTypeCore;
			ftype len = length(p);
			if (len < static_cast<ftype>(epslion))
				len = static_cast<ftype>(epslion);
			return divd(p, len);
		}
		static void zero(TVert& tp) {	
			for (int i = 0; i < dim; ++i)
				tp[i] = ftype();
		}
	public:
		struct RPoint :public std::array<ftype, dim>{
			RPoint(const TVert& p) {
				for (int i = 0; i < dim; ++i)
					(*this)[i] = p[i];
			}
			bool operator<(const RPoint& that) {
				for (int i = 0; i < dim; ++i)
					if ((*this)[i] < that[i])
						return true;
				return false;
			}
		};

		static void repeat(const std::vector<TVert>& pts, std::vector<TVert>& dpt)
		{
			int npt = static_cast<int>(pts.size());
			dpt.clear();
			dpt.reserve(npt);
			if (npt > 0)
				dpt.push_back(pts[0]);
			for (int i = 1; i < npt; ++i) {
				int nd = static_cast<int>(dpt.size());
				if (length2(dpt[nd - 1], pts[i]) > static_cast<ftype>(BasicTypeCore::epslion)) {
					dpt.push_back(pts[i]);
				}
			}
		}

	public: //曲线参数化算法

		//均匀参数化
		static void getAverage(int npts, std::vector<ftype>& u)
		{
			assert(npts >= 2 && "npts must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			for (int i = 0; i < npts; ++i) {
				u.push_back(static_cast<ftype>(i / (ftype)(npts - 1)));
			}
		}
	public:
		//积累弦长参数化
		static void getChordLength(const std::vector<TVert>& pts, std::vector<ftype>& u)
		{
			using namespace BasicTypeCore;
			assert(pts.size() >= 2 && "pts size must be max or equal 2");
			size_t npt = pts.size();
			u.clear();
			u.reserve(npt);
			u.push_back(0);
			ftype clen = static_cast<ftype>(0);
			for (int i = 1; i < npt; ++i) {
				clen += length(pts[i - 1], pts[i]);
				u.push_back(clen);
			}
			for (int i = 0; i < npt; ++i)
				u[i] /= clen;
		}

		//向心参数化
		static void getCentripetal(const std::vector<TVert>& pts, std::vector<ftype>& u)
		{
			using namespace BasicTypeCore;
			assert(pts.size() >= 2 && "pts size must be max or equal 2");
			size_t npt = pts.size();
			u.clear();
			u.reserve(npt);
			u.push_back(0);
			ftype clen = static_cast<ftype>(0);
			for (int i = 1; i < npt; ++i) {
				clen += static_cast<ftype>(sqrt(length(pts[i - 1], pts[i])));
				u.push_back(clen);
			}
			for (int i = 0; i < npt; ++i)
				u[i] /= clen;
		}

		//指数参数化
		static void getExp(const std::vector<TVert>& pts, std::vector<ftype>& u,
			int degree = 3, ftype alpha = static_cast<ftype>(0.8))
		{
			using namespace BasicTypeCore;
			assert(pts.size() >= 2 && "pts size must be max or equal 2");
			size_t npt = pts.size();
			u.clear();
			u.reserve(npt);
			u.push_back(0);
			ftype clen = static_cast<ftype>(0);
			for (int i = 1; i < npt; ++i) {
				clen += static_cast<ftype>(std::pow(length(pts[i - 1], pts[i]), alpha));
				u.push_back(clen);
			}
			for (int i = 0; i < npt; ++i)
				u[i] /= clen;
		}

		//混合参数化
		static void getHybrid(const std::vector<TVert>& pts,
			std::vector<ftype>& u, int degree = 3, ftype alpha = static_cast<ftype>(0.8))
		{
			using namespace BasicTypeCore;
			assert(pts.size() >= 2 && "pts size must be max or equal 2");
			assert(alpha > -1 && "alpha must be max -1!");
			size_t npt = pts.size(); u.clear();
			u.reserve(npt); u.push_back(0);
			ftype clen = static_cast<ftype>(0);
			for (int i = 1; i < npt; ++i) {
				clen += static_cast<ftype>(std::pow(length(pts[i - 1], pts[i]), alpha));
				u.push_back(clen);
			}
			for (int i = 0; i < npt; ++i)
				u[i] /= clen;
			for (int i = 0; i < npt; ++i) {
				u[i] = (1 + alpha)*u[i] / (1 + alpha * u[i]);
			}

		}

		//福利参数化（目前最好的参数化）
		static void getFoley(const std::vector<TVert>& pts, std::vector<ftype>& u)
		{
			using namespace BasicTypeCore;
			assert(pts.size() >= 2 && "pts size must be max or equal 2");
			auto FoleyK = [&pts](int i) {
				if (i == 1) {
					TVert p0 = sub(pts[0], pts[1]);
					TVert p1 = sub(pts[2], pts[1]);
					ftype t0 = length(p0);
					ftype t1 = length(p1);
					if (t0 < static_cast<ftype>(epslion))
						t0 = static_cast<ftype>(epslion);
					if (t1 < static_cast<ftype>(epslion))
						t1 = static_cast<ftype>(epslion);
					ftype ff = dot(p0, p1) / (t0*t1);
					if (ff < -1)
						ff = -1;
					if (ff > 1)
						ff = 1;
					ftype angle = static_cast<ftype>(pi) - std::acos(ff);
					angle = std::fmin(angle, static_cast<ftype>(pi) / 2);
					return static_cast<ftype>(1 + 1.5*(angle*t0 / (t0 + t1)));
				}
				else if (i < pts.size() - 1) {
					TVert pf2 = sub(pts[i - 2], pts[i - 1]);
					TVert pf1 = sub(pts[i], pts[i - 1]);
					TVert p1 = muilt(-1, pf1);
					TVert p2 = sub(pts[i + 1], pts[i]);
					ftype lenf2 = length(pf2);
					ftype lenf1 = length(pf1);
					ftype len1 = length(p1);
					ftype len2 = length(p2);
					if (lenf2 < static_cast<ftype>(epslion))
						lenf2 = static_cast<ftype>(epslion);
					if (lenf1 < static_cast<ftype>(epslion))
						lenf1 = static_cast<ftype>(epslion);
					if (len2 < static_cast<ftype>(epslion))
						len2 = static_cast<ftype>(epslion);
					if (len1 < static_cast<ftype>(epslion))
						len1 = static_cast<ftype>(epslion);
					ftype anglef1 = dot(pf2, pf1) / (lenf1*lenf2);
					ftype angle = dot(p1, p2) / (len1*len2);
					if (anglef1 < -1)anglef1 = -1;
					if (anglef1 > 1)anglef1 = 1;
					if (angle < -1)angle = -1;
					if (angle > 1)angle = 1;
					anglef1 = static_cast<ftype>(pi - std::acos(anglef1));
					anglef1 = static_cast<ftype>(std::fmin(anglef1, pi / 2));
					angle = static_cast<ftype>(pi - std::acos(angle));
					angle = static_cast<ftype>(std::fmin(angle, pi / 2));
					return static_cast<ftype>(1 + 1.5*(lenf2*anglef1 / (lenf1 + lenf2) + len2 * angle / (len1 + len2)));
				}
				else
				{
					TVert pf2 = sub(pts[i - 2], pts[i - 1]);
					TVert pf1 = sub(pts[i], pts[i - 1]);
					ftype lenf2 = length(pf2);
					ftype lenf1 = length(pf1);
					if (lenf2 < static_cast<ftype>(epslion))
						lenf2 = static_cast<ftype>(epslion);
					if (lenf1 < static_cast<ftype>(epslion))
						lenf1 = static_cast<ftype>(epslion);
					ftype anglef1 = dot(pf2, pf1) / (lenf1*lenf2);
					if (anglef1 < -1)anglef1 = -1;
					if (anglef1 > 1)anglef1 = 1;
					anglef1 = static_cast<ftype>(pi - std::acos(anglef1));
					anglef1 = static_cast<ftype>(std::fmin(anglef1, pi / 2));
					return static_cast<ftype>(1 + 1.5*(lenf2*anglef1 / (lenf1 + lenf2)));
				}
			};

			size_t npt = pts.size();
			u.clear();
			u.reserve(npt);
			u.push_back(0);
			ftype clen = static_cast<ftype>(0);
			for (int i = 1; i < npt; ++i) {
				clen += length(pts[i - 1], pts[i])*FoleyK(i);
				u.push_back(clen);
			}
			for (int i = 0; i < npt; ++i)
				u[i] /= clen;
		}
	public://样条函数内核函数
		static int findSpan(size_t degree, const std::vector<ftype> &knots, ftype u) {
			int n = static_cast<int>(knots.size()) - degree - 2;
			if (u > (knots[n + 1] - std::numeric_limits<ftype>::epsilon())) {
				return n;
			}
			if (u < (knots[degree] + std::numeric_limits<ftype>::epsilon())) {
				return degree;
			}
			int low = degree;
			int high = n + 1;
			int mid = (int)std::floor((low + high) / 2.0);
			while (u < knots[mid] || u >= knots[mid + 1]) {
				if (u < knots[mid]) {
					high = mid;
				}
				else {
					low = mid;
				}
				mid = (int)std::floor((low + high) / 2.0);
			}
			return mid;
		}

		/**
		 * Check if two numbers are close enough within eps
		 * @param[in] a First number
		 * @param[in] b Second number
		 * @param[in] eps Tolerance for checking closeness
		 * @return Whether the numbers are close w.r.t. the tolerance
		 */
		static bool close(ftype a, ftype b, double eps = std::numeric_limits<ftype>::epsilon()) {
			return (std::abs(a - b) < eps) ? true : false;
		}
		/**
		 * Compute a single B-spline basis function
		 * @param[in] i The ith basis function to compute.
		 * @param[in] deg Degree of the basis function.
		 * @param[in] knots Knot vector corresponding to the basis functions.
		 * @param[in] u Parameter to evaluate the basis functions at.
		 * @return The value of the ith basis function at u.
		 */
		static ftype  bsplineOneBasis(int i, size_t deg, const std::vector<ftype> &U, ftype u)
		{
			int m = static_cast<int>(U.size()) - 1;
			if ((i == 0 && close(u, U[0])) || (i == m - deg - 1 && close(u, U[m]))) {
				return static_cast<ftype>(1.0);
			}
			if (u < U[i] || u >= U[i + deg + 1]) {
				return static_cast<ftype>(0.0);
			}
			std::vector<ftype> N;
			N.resize(deg + 1);
			for (int j = 0; j <= deg; j++) {
				N[j] = static_cast<ftype>((u >= U[i + j]
					&& u < U[i + j + 1]) ? static_cast<ftype>(1.0) : static_cast<ftype>(0.0));
			}
			for (int k = 1; k <= deg; k++) {
				ftype saved = (close(N[0], static_cast<ftype>(0.0))) ?
					static_cast<ftype>(0.0) : ((u - U[i]) * N[0]) / (U[i + k] - U[i]);
				for (int j = 0; j < deg - k + 1; j++) {
					ftype Uleft = U[i + j + 1];
					ftype Uright = U[i + j + k + 1];
					if (close(N[j + 1], static_cast<ftype>(0.0))) {
						N[j] = saved;
						saved = static_cast<ftype>(0.0);
					}
					else {
						ftype temp = N[j + 1] / (Uright - Uleft);
						N[j] = saved + (Uright - u) * temp;
						saved = (u - Uleft) * temp;
					}
				}
			}
			return N[0];
		}

		/**
		 * Compute all non-zero B-spline basis functions
		 * @param[in] deg Degree of the basis function.
		 * @param[in] span Index obtained from findSpan() corresponding the u and knots.
		 * @param[in] knots Knot vector corresponding to the basis functions.
		 * @param[in] u Parameter to evaluate the basis functions at.
		 * @return N Values of (deg+1) non-zero basis functions.
		 */
		static std::vector<ftype> bsplineBasis(size_t deg,
			int span, const std::vector<ftype> &knots, ftype u)
		{
			std::vector<ftype> N;
			N.resize(deg + 1, ftype(0));
			std::vector<ftype> left, right;
			left.resize(deg + 1, static_cast<ftype>(0.0));
			right.resize(deg + 1, static_cast<ftype>(0.0));
			ftype saved = ftype(), temp = ftype();
			N[0] = static_cast<ftype>(1.0);
			for (int j = 1; j <= deg; j++) {
				left[j] = (u - knots[span + 1 - j]);
				right[j] = knots[span + j] - u;
				saved = 0.0;
				for (int r = 0; r < j; r++) {
					temp = N[r] / (right[r + 1] + left[j - r]);
					N[r] = saved + right[r + 1] * temp;
					saved = left[j - r] * temp;
				}
				N[j] = saved;
			}
			return N;
		}

		/**
		 * Compute all non-zero derivatives of B-spline basis functions
		 * @param[in] deg Degree of the basis function.
		 * @param[in] span Index obtained from findSpan() corresponding the u and knots.
		 * @param[in] knots Knot vector corresponding to the basis functions.
		 * @param[in] u Parameter to evaluate the basis functions at.
		 * @param[in] num_ders Number of derivatives to compute (num_ders <= deg)
		 * @return ders Values of non-zero derivatives of basis functions.
		 */
		static array2<ftype> bsplineDerBasis(size_t deg,
			int span, const std::vector<ftype> &knots, ftype u, int num_ders)
		{
			std::vector<ftype> left, right;
			left.resize(deg + 1, ftype(0));
			right.resize(deg + 1, ftype(0));
			ftype saved = ftype(0), temp = ftype(0);
			array2<ftype> ndu(deg + 1, deg + 1);
			ndu(0, 0) = static_cast<ftype>(1.0);
			for (int j = 1; j <= deg; ++j) {
				left[j] = u - knots[span + 1 - j];
				right[j] = knots[span + j] - u;
				saved = ftype(0);
				for (int r = 0; r < j; ++r) {
					ndu(j, r) = right[r + 1] + left[j - r];
					temp = ndu(r, j - 1) / ndu(j, r);
					ndu(r, j) = saved + right[r + 1] * temp;
					saved = left[j - r] * temp;
				}
				ndu(j, j) = saved;
			}
			array2<ftype> ders(num_ders + 1, deg + 1, ftype(0));
			for (int j = 0; j <= deg; ++j) {
				ders(0, j) = ndu(j, deg);
			}
			array2<ftype> a(2, deg + 1);
			for (int r = 0; r <= deg; ++r) {
				int s1 = 0;
				int s2 = 1;
				a(0, 0) = static_cast<ftype>(1.0);
				for (int k = 1; k <= num_ders; ++k) {
					ftype d = static_cast<ftype>(0.0);
					int rk = r - k;
					int pk = deg - k;
					int j1 = 0;
					int j2 = 0;
					if (r >= k) {
						a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
						d = a(s2, 0) * ndu(rk, pk);
					}
					if (rk >= -1)
						j1 = 1;
					else
						j1 = -rk;
					if (r - 1 <= pk)
						j2 = k - 1;
					else
						j2 = deg - r;
					for (int j = j1; j <= j2; ++j) {
						a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
						d += a(s2, j) * ndu(rk + j, pk);
					}
					if (r <= pk) {
						a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
						d += a(s2, k) * ndu(r, pk);
					}
					ders(k, r) = d;
					int temp = s1;
					s1 = s2;
					s2 = temp;
				}
			}
			ftype fac = static_cast<ftype>(deg);
			for (int k = 1; k <= num_ders; k++) {
				for (int j = 0; j <= deg; j++) {
					ders(k, j) *= fac;
				}
				fac *= static_cast<ftype>(deg - k);
			}
			return ders;
		}
	public://Spline Eval
		/**
		 * Evaluate point on a nonrational NURBS curve
		 * @param[in] degree Degree of the given curve.
		 * @param[in] knots Knot vector of the curve.
		 * @param[in] control_points Control points of the curve.
		 * @param[in] u Parameter to evaluate the curve at.
		 * @return point Resulting point on the curve at parameter u.
		 */
		static TVert curvePoint(size_t degree, const std::vector<ftype> &knots,
			const std::vector<TVert> &control_points, ftype u)
		{
			TVert p;
			for (int i = 0; i < dim; ++i)
				p[i] = 0;
			int span = findSpan(degree, knots, u);
			std::vector<ftype> N = bsplineBasis(degree, span, knots, u);
			for (size_t j = 0; j <= degree; ++j) {
				TVert cp = muilt(N[j], control_points[span - degree + j]);
				p = add(p, cp);
			}
			return p;
		}

		/**
		 * Evaluate derivatives of a non-rational NURBS curve
		 * @param[in] degree Degree of the curve
		 * @param[in] knots Knot vector of the curve.
		 * @param[in] control_points Control points of the curve.
		 * @param[in] num_ders Number of times to derivate.
		 * @param[in] u Parameter to evaluate the derivatives at.
		 * @return curve_ders Derivatives of the curve at u.
		 * E.g. curve_ders[n] is the nth derivative at u, where 0 <= n <= num_ders.
		 */
		static std::vector<TVert> curveDerivatives(size_t degree, const std::vector<ftype> &knots,
			const std::vector<TVert> &control_points, int num_ders, ftype u)
		{
			using std::vector;
			std::vector<TVert> curve_ders;
			curve_ders.resize(num_ders + 1);
			for (int k = degree + 1; k <= num_ders; ++k) {
				for (int d = 0; d < dim; ++d)
					curve_ders[k][d] = static_cast<ftype>(0);
			}
			int span = findSpan(degree, knots, u);
			array2<ftype> ders = bsplineDerBasis(degree, span, knots, u, num_ders);
			int du = num_ders < degree ? num_ders : degree;
			for (int k = 0; k <= du; ++k) {
				for (int j = 0; j <= degree; ++j) {
					TVert cp = muilt(ders(k, j), control_points[span - degree + j]);
					curve_ders[k] = add(curve_ders[k], cp);
				}
			}
			return curve_ders;
		}

		/**
		 * Evaluate the tangent of a B-spline curve
		 * @param[in] crv Curve object
		 * @return Unit tangent of the curve at u.
		 */
		static TVert curveTangent(size_t degree, const std::vector<ftype> &knots,
			const std::vector<TVert> &control_points, ftype u)
		{
			std::vector<TVert> ders = curveDerivatives(degree, knots, control_points, 1, u);
			TVert du = ders[1];
			ftype du_len = length(du);
			if (!close(du_len, ftype(0))) {
				du /= du_len;
			}
			return du;
		}

	public://nurbs Eval
		/**
		 * Evaluate point on a rational NURBS curve
		 * @param[in] knots RationalCurve knotes
		 * @param[in] control_points RationalCurve control_points
		 * @param[in] weights RationalCurve weights
		 * @param[in] u Parameter to evaluate the curve at.
		 * @return point Resulting point on the curve.
		 */
		static TVert curvePoint(size_t degree, const std::vector<ftype> &knots,
			const std::vector<TVert> &control_points, const std::vector<ftype> &weights, ftype u)
		{
			using tvecnp1 = std::array<ftype, dim + 1>;
			std::vector<tvecnp1> Cw;
			int ncp = static_cast<int>(control_points.size());
			Cw.reserve(ncp);
			for (int i = 0; i < ncp; ++i) {
				tvecnp1 tp;
				for (int j = 0; j < dim; ++j)
					tp[j] = weights[i] * control_points[i][j];
				tp[dim] = weights[i];
				Cw.emplace_back(tp);
			}
			tvecnp1 pointw = AlgoUtil<tvecnp1, ftype, dim + 1>::curvePoint(degree, knots, Cw, u);
			TVert tp;
			for (int i = 0; i < dim; ++i) {
				tp[i] = pointw[i] / pointw[dim];
			}
			return tp;
		}

		/**
		 * Evaluate derivatives of a rational NURBS curve
		 * @param[in] u Parameter to evaluate the derivatives at.
		 * @param[in] knots Knot vector of the curve.
		 * @param[in] control_points Control points of the curve.
		 * @param[in] weights Weights corresponding to each control point.
		 * @param[in] num_ders Number of times to differentiate.
		 * @param[inout] curve_ders Derivatives of the curve at u.
		 * E.g. curve_ders[n] is the nth derivative at u, where n is between 0 and
		 * num_ders-1.
		 */
		static std::vector<TVert> curveDerivatives(
			size_t degree,
			const std::vector<ftype> &knots,
			const std::vector<TVert> &control_points,
			const std::vector<ftype> &weights, int num_ders, ftype u)
		{
			using tvecnp1 = std::array<ftype, dim + 1>;
			std::vector<TVert> curve_ders;
			curve_ders.reserve(num_ders + 1);
			std::vector<tvecnp1> Cw;
			int ncp = static_cast<int>(control_points.size());
			Cw.reserve(ncp);
			for (int i = 0; i < ncp; ++i) {
				tvecnp1 tp;
				for (int j = 0; j < dim; ++j)
					tp[j] = weights[i] * control_points[i][j];
				tp[dim] = weights[i];
				Cw.emplace_back(tp[dim]);
			}
			std::vector<tvecnp1> Cwders = AlgoUtil<tvecnp1, ftype, dim + 1>::curveDerivatives(degree, knots, Cw, num_ders, u);
			std::vector<TVert> Aders;
			std::vector<ftype> wders;
			for (const auto &val : Cwders) {
				TVert p;
				for (int i = 0; i < dim; ++i)
					p[i] = val[i];
				Aders.emplace_back(p);
				wders.emplace_back(val[dim]);
			}
			for (int k = 0; k <= num_ders; ++k) {
				TVert v = Aders[k];
				for (int i = 1; i <= k; ++i) {
					TVert tt = muilt(binomial(k, i)*wders[i], curve_ders[k - i]);
					v = sub(v, tt);
				}
				curve_ders.emplace_back(divd(v, wders[0]));
			}

		}

		/**
		 * Evaluate the tangent of a rational B-spline curve
		 * @param[in] crv RationalCurve object
		 * @return Unit tangent of the curve at u.
		 */
		static TVert curveTangent(size_t degree,
			const std::vector<ftype> &knots,
			const std::vector<TVert> &control_points,
			const std::vector<ftype> &weights, ftype u)
		{
			std::vector<TVert> ders = curveDerivatives(degree, knots, control_points, weights, 1, u);
			TVert du = ders[1];
			ftype du_len = length(du);
			if (!close(du_len, ftype(0))) {
				du /= du_len;
			}
			return du;
		}

	private:
		/**
		 * Compute the binomial coefficient (nCk) using the formula
		 * \product_{i=0}^k (n + 1 - i) / i
		 */
		static size_t binomial(size_t n, size_t k)
		{
			size_t result = 1;
			if (k > n) {
				return 0;
			}
			for (unsigned int i = 1; i <= k; ++i) {
				result *= (n + 1 - i);
				result /= i;
			}
			return result;
		}

	public://曲线求导算法
		//前差值法
		static TVert getSub(const std::vector<TVert>& pts, const std::vector<ftype>& uparam, int i)
		{
			assert((pts.size() > 2 && i >= 0 && i < pts.size()) && "the size of pts must be [2,+OO),i must be [0,pts.size())!");
			using namespace BasicTypeCore;
			if (i < pts.size() - 1) {
				TVert tp = sub(pts[i + 1], pts[i]);
				tp = divd(tp, length(tp));
				return tp;
			}
			else {
				TVert tp = sub(pts[i], pts[i - 1]);
				tp = divd(tp, length(tp));
				return tp;
			}
		}

		//弗米尔求导
		static TVert getFmill(const std::vector<TVert>& pts, const std::vector<ftype>& uparam, int i) {
			assert((pts.size() > 2 && i >= 0 && i < pts.size()) && "the size of pts must be [2,+OO),i must be [0,pts.size())!");
			using namespace BasicTypeCore;
			if (i == 0) {
				TVert dp = divd(sub(pts[1], pts[0]), uparam[1] - uparam[0]);
				ftype len = length(dp);
				if (len < static_cast<ftype>(epslion))
					len = static_cast<ftype>(epslion);
				dp = divd(dp, length(dp));
				return dp;
			}
			else if (i > 0 && i < pts.size() - 1) {
				TVert p = sub(pts[i + 1], pts[i - 1]);
				p = divd(p, uparam[i + 1] - uparam[i - 1]);
				return p;
			}
			else {
				TVert dp = divd(sub(pts[pts.size() - 1], pts[pts.size() - 2]), uparam[pts.size() - 1] - uparam[pts.size() - 2]);
				return dp;
			}

		}
		//贝塞尔求导
		static TVert getBessel(const std::vector<TVert>& pts, const std::vector<ftype>& uparam, int i) {
			assert((pts.size() > 2 && i >= 0 && i < pts.size()) && "the size of pts must be [2,+OO),i must be [0,pts.size())!");
			if (i == 0) {
				TVert dp = divd(sub(pts[1], pts[0]), uparam[1] - uparam[0]);
				dp = muilt(static_cast<ftype>(2), dp);
				TVert dp1 = muilt((uparam[2] - uparam[1]) / (uparam[2] - uparam[0]), divd(sub(pts[1], pts[0]), uparam[1] - uparam[0]));
				dp1 = add(dp1, muilt((uparam[1] - uparam[0]) / (uparam[2] - uparam[0]), divd(sub(pts[2], pts[1]), uparam[2] - uparam[1])));
				dp = sub(dp, dp1);
				return dp;
			}
			else if (i > 0 && i < pts.size() - 1) {
				TVert dp1 = muilt((uparam[i + 1] - uparam[i]) / (uparam[i + 1] - uparam[i - 1]), divd(sub(pts[i], pts[i - 1]), uparam[i] - uparam[i - 1]));
				dp1 = add(dp1, muilt((uparam[i] - uparam[i - 1]) / (uparam[i + 1] - uparam[i - 1]), divd(sub(pts[i + 1], pts[i]), uparam[i + 1] - uparam[i])));
				return dp1;
			}
			else {
				int npt = static_cast<int>(pts.size());
				TVert dp = divd(sub(pts[npt - 1], pts[npt - 2]), uparam[npt - 1] - uparam[npt - 2]);
				dp = muilt(static_cast<ftype>(2), dp);
				TVert dp1 = muilt((uparam[i] - uparam[i - 1]) / (uparam[i] - uparam[i - 2]), divd(sub(pts[i - 1], pts[i - 2]), uparam[i - 1] - uparam[i - 2]));
				dp1 = add(dp1, muilt((uparam[i - 1] - uparam[i - 2]) / (uparam[i] - uparam[i - 2]), divd(sub(pts[i], pts[i - 1]), uparam[i] - uparam[i - 1])));
				TVert dp2 = sub(dp, dp1);
				return dp2;
			}
		}
	};



	//曲线参数化的方式
	template<class TVert, class ftype, int dim = 3>
	static void getCurveParam(const std::vector<TVert>& pts, std::vector<ftype>& u,
		ParamType param=ParamType::idChordLength, int degree = 3, ftype alpha = static_cast<ftype>(0.8)) {
		switch (param)
		{
		case FussenAlgo::ParamType::idAverage:
			AlgoUtil<TVert, ftype, dim>::getAverage(static_cast<int>(pts.size()), u);
			break;
		case FussenAlgo::ParamType::idChordLength:
			AlgoUtil<TVert, ftype, dim>::getChordLength(pts, u);
			break;
		case FussenAlgo::ParamType::idCentripetal:
			AlgoUtil<TVert, ftype, dim>::getCentripetal(pts, u);
			break;
		case FussenAlgo::ParamType::idFoley:
			AlgoUtil<TVert, ftype, dim>::getFoley(pts, u);
			break;
		case FussenAlgo::ParamType::idHybrid:
			AlgoUtil<TVert, ftype, dim>::getHybrid(pts, u, degree, alpha);
			break;
		case FussenAlgo::ParamType::idExp:
			AlgoUtil<TVert, ftype, dim>::getExp(pts, u, degree, alpha);
			break;
		}
	}
	//离散点求导的
	template<class TVert, class ftype, int dim = 3>
	static TVert getCurveDeriative(const std::vector<TVert>& pts,
		const std::vector<ftype>& uparam, int i, DerivativeType dt = DerivativeType::idBessel) {
		switch (dt)
		{
		case FussenAlgo::DerivativeType::idFmill:
			return AlgoUtil<TVert, ftype, dim>::getFmill(pts, uparam, i);
		case FussenAlgo::DerivativeType::idBessel:
			return AlgoUtil<TVert, ftype, dim>::getBessel(pts, uparam, i);
		case FussenAlgo::DerivativeType::idSub:
			return AlgoUtil<TVert, ftype, dim>::getSub(pts, uparam, i);
		}
	}
	//点数分布的方式
	template<class ftype, int dim = 3>
	struct distPoints {
		//均匀分布
		static void getAverage(int npts, std::vector<ftype>& u)
		{
			assert(npts >= 2 && "npts must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			for (int i = 0; i < npts; ++i) {
				u.push_back(static_cast<ftype>(i / (ftype)(npts - 1)));
			}
		}
		//线性分布
		static void getLinear(int npts, std::vector<ftype>& u, ftype ek = static_cast<ftype>(1.0))
		{
			assert(npts >= 2 && "npts must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			ftype dx = 1 / (npts - 1);
			ftype dlen = static_cast<ftype>(0);
			u.push_back(dlen);
			for (int i = 1; i < npts; ++i) {
				dlen += static_cast<ftype>(ek*i*dx);
				u.push_back(dlen);
			}
			for (int i = 0; i < npts; ++i)
				u[i] /= dlen;
		}

		//指数分布
		static void getExp(int npts, std::vector<ftype>& u, ftype ep = ftype(1),
			ftype xmin = 0, ftype xmax = static_cast<ftype>(2.0))
		{
			assert(npts >= 2 && "npts must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			ftype dx = (xmax - xmin) / (npts - 1);
			ftype dlen = static_cast<ftype>(0);
			u.push_back(dlen);
			ftype c = std::exp(ep*xmin);
			for (int i = 1; i < npts; ++i) {
				dlen += static_cast<ftype>(std::exp(ep*(xmin + i * dx)) - c);
				u.push_back(dlen);
			}
			for (int i = 0; i < npts; ++i)
				u[i] /= dlen;
		}
		//抛物线分布规律
		static void getParabola(int npts, std::vector<ftype>&u, ftype ea = static_cast<ftype>(1),
			ftype xmin = 0, ftype xmax = static_cast<ftype>(2.0))
		{
			assert(npts >= 2 && "npts must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			ftype dx = (xmax - xmin) / (npts - 1);
			ftype dlen = static_cast<ftype>(0);
			u.push_back(dlen);
			ftype c = ea * xmin*xmin;
			for (int i = 1; i < npts; ++i) {
				ftype cu = xmin + i * dx;
				dlen += static_cast<ftype>(std::abs(ea*cu*cu - c));
				u.push_back(dlen);
			}
			for (int i = 0; i < npts; ++i)
				u[i] /= dlen;
		}
		//拷贝分布
		template<class TVert>
		static void copy(const std::vector<TVert>& ori_pts, std::vector<ftype>&u) {
			int npts = static_cast<int>(ori_pts.size());
			assert(npts >= 2 && "ori_pts size must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			ftype dlen = static_cast<ftype>(0);
			u.push_back(dlen);
			for (int i = 1; i < npts; ++i) {
				dlen += AlgoUtil<TVert, ftype, dim>::length(ori_pts[i - 1], ori_pts[i]);
				u.push_back(dlen);
			}
			for (int i = 0; i < npts; ++i)
				u[i] /= dlen;
		}
		//龚珀兹曲线分布规律
		static void getGompertz(int npts, std::vector<ftype>& u,
			ftype ek = static_cast<ftype>(1), ftype eb = static_cast<ftype>(1),
			ftype xmin = 0, ftype xmax = static_cast<ftype>(2.0))
		{
			assert(npts >= 2 && "npts must be max or equal 2!");
			u.clear();
			u.reserve(npts);
			ftype dx = (xmax - xmin) / (npts - 1);
			ftype dlen = static_cast<ftype>(0);
			u.push_back(dlen);
			ftype c = ek * std::exp(std::pow(eb, xmin));
			for (int i = 1; i < npts; ++i) {
				ftype cu = xmin + i * dx;
				dlen += static_cast<ftype>(ek*std::exp(std::pow(eb, cu)) - c);
				u.push_back(dlen);
			}
			for (int i = 0; i < npts; ++i)
				u[i] /= dlen;
		}

	};


	/*
	 * @struct weightHermiteCubic 带权重的埃米特插值
	 * @detail 只有两个点，区间为[0,1]的插值
	 * @author lyc
	 * @date 2020.2.22
	 */
	template<class TVert, class ftype, int dim = 3>
	struct weightHermiteCubic {
		//using ftype = float;
		//static const int dim = 3;
		//using TVert = std::array<ftype, 3>;
	private:
		TVert mps;
		TVert mpe;
		TVert mpsd;
		TVert mped;
	public:
		void interpolation(const TVert& ps, const TVert& pe, const TVert& psd, const TVert& ped) {
			mps = ps;
			mpe = pe;
			mpsd = psd;
			mped = ped;
		}

		TVert evalunify(ftype u, ftype weights = static_cast<ftype>(1.0), ftype weighte = static_cast<ftype>(1.0)) {
			assert((u >= 0 && u <= 1) && "u must be [0,1]!");
			TVert a0 = mps;
			TVert a1 = AlgoUtil<TVert, ftype, dim>::muilt(weights, mpsd);
			TVert ae = AlgoUtil<TVert, ftype, dim>::muilt(weighte, mped);
			TVert a2 = AlgoUtil<TVert, ftype, dim>::add(
				AlgoUtil<TVert, ftype, dim>::muilt(-3, mps),
				AlgoUtil<TVert, ftype, dim>::muilt(-2, a1));
			a2 = AlgoUtil<TVert, ftype, dim>::add(a2,
				AlgoUtil<TVert, ftype, dim>::muilt(-1, ae));
			a2 = AlgoUtil<TVert, ftype, dim>::add(a2,
				AlgoUtil<TVert, ftype, dim>::muilt(3, mpe));
			TVert a3 = AlgoUtil<TVert, ftype, dim>::add(AlgoUtil<TVert, ftype, dim>::muilt(2, mps), a1);
			a3 = AlgoUtil<TVert, ftype, dim>::add(a3, ae);
			a3 = AlgoUtil<TVert, ftype, dim>::add(a3,
				AlgoUtil<TVert, ftype, dim>::muilt(-2, mpe));
			TVert ap;
			for (int i = 0; i < dim; ++i)
				ap[i] = static_cast<ftype>(0.0);
			ap = AlgoUtil<TVert, ftype, dim>::add(ap, AlgoUtil<TVert, ftype, dim>::muilt(static_cast<ftype>(1.0), a0));
			ap = AlgoUtil<TVert, ftype, dim>::add(ap, AlgoUtil<TVert, ftype, dim>::muilt(u, a1));
			ap = AlgoUtil<TVert, ftype, dim>::add(ap, AlgoUtil<TVert, ftype, dim>::muilt(u*u, a2));
			ap = AlgoUtil<TVert, ftype, dim>::add(ap, AlgoUtil<TVert, ftype, dim>::muilt(u*u*u, a3));
			return ap;
		}

		TVert evalArbit(ftype u, ftype umin, ftype umax, ftype weights = static_cast<ftype>(1.0), ftype weighte = static_cast<ftype>(1.0))
		{
			assert((u >= umin && u <= umax) && "u must be [umin,umax]!");
			return evalunify(((u - umin) / (umax - umin)), weights, weighte);
		}
	};

	/*
	 * @struct CubeicSpline 最优带约束三次B样条插值
	 * @detail 只有两个点，区间为[0,1]的插值
	 * @author lyc
	 * @date 2020.2.22
	 */
	template<class TVert, class ftype, int dim = 3>
	struct CubeicCurve {
		//using ftype = float;
		//static const int dim = 3;
		//using TVert = std::array<ftype, 3>;
	private:
		TVert t0, t1, t2, t3;
	public:
		void interpolation(const TVert& ps, const TVert& pe, const TVert& psd, const TVert& ped) {
			using namespace BasicTypeCore;
			ftype cd0 = AlgoUtil<TVert, ftype, dim>::length(psd);
			ftype cd1 = AlgoUtil<TVert, ftype, dim>::length(ped);
			if (cd0 < static_cast<ftype>(epslion))
				cd0 = static_cast<ftype>(epslion);
			if (cd1 < static_cast<ftype>(epslion))
				cd1 = static_cast<ftype>(epslion);
			ftype cval = AlgoUtil<TVert, ftype, dim>::dot(psd, ped) / (cd0*cd1);
			if (cval > static_cast<ftype>(1.0))
				cval = static_cast<ftype>(1.0);
			if (cval < static_cast<ftype>(-1.0))
				cval = static_cast<ftype>(-1.0);
			ftype alpha = static_cast<ftype>(acos(cval));
			if (alpha < static_cast<ftype>(epslion))
				alpha = static_cast<ftype>(epslion);
			ftype dsend = AlgoUtil<TVert, ftype, dim>::length(ps, pe);
			if (dsend < static_cast<ftype>(epslion))
				dsend = static_cast<ftype>(epslion);
			dsend = static_cast<ftype>(dsend * (((0.5 / sin(alpha / 2.0)) -
				(0.5 / tan(alpha / 2.0)))*4.0) / sin(alpha / 2.0));
			TVert sdir = AlgoUtil<TVert, ftype, dim>::muilt(dsend / cd0, psd);
			TVert edir = AlgoUtil<TVert, ftype, dim>::muilt(dsend / cd1, ped);
			t0 = ps;
			t1 = sdir;
			t2 = AlgoUtil<TVert, ftype, dim>::sub(pe, ps);
			t2 = AlgoUtil<TVert, ftype, dim>::muilt(static_cast<ftype>(3.0), t2);
			t2 = AlgoUtil<TVert, ftype, dim>::sub(t2, edir);
			t2 = AlgoUtil<TVert, ftype, dim>::sub(t2,
				AlgoUtil<TVert, ftype, dim>::muilt(static_cast<ftype>(2.0), sdir));
			t3 = AlgoUtil<TVert, ftype, dim>::add(sdir, edir);
			t3 = AlgoUtil<TVert, ftype, dim>::add(t3,
				AlgoUtil<TVert, ftype, dim>::muilt(static_cast<ftype>(-2.0), pe));
			t3 = AlgoUtil<TVert, ftype, dim>::add(t3,
				AlgoUtil<TVert, ftype, dim>::muilt(static_cast<ftype>(2.0), ps));
		}
		TVert eval(ftype u) {
			using namespace BasicTypeCore;
			assert((u >= 0 && u <= 1) && "u must be [0,1]!");
			TVert p = t0;
			p = AlgoUtil<TVert, ftype, dim>::add(p, AlgoUtil<TVert, ftype, dim>::muilt(u, t1));
			p = AlgoUtil<TVert, ftype, dim>::add(p, AlgoUtil<TVert, ftype, dim>::muilt(u*u, t2));
			p = AlgoUtil<TVert, ftype, dim>::add(p, AlgoUtil<TVert, ftype, dim>::muilt(u*u*u, t3));
			return p;
		}
		TVert eval(ftype u, ftype umin, ftype umax) {
			assert((u >= umin && u <= umax) && "u must be [umin,umax]!");
			return eval((u - umin) / (umax - umin));
		}

		//因为标准借口不可靠
		TVert curveTanght(ftype u) {
			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			TVert p;
			if (u <= static_cast<ftype>(0.999)) {
				p = AlgoUtil<TVert, ftype, dim>::sub(eval(u + static_cast<ftype>(0.001)), eval(u));
				p = AlgoUtil<TVert, ftype, dim>::normal(p);
			}
			else {
				p = AlgoUtil<TVert, ftype, dim>::sub(eval(u), eval(u - static_cast<ftype>(0.001)));
				p = AlgoUtil<TVert, ftype, dim>::normal(p);
			}
			return p;
		}
	};

	/*
	 * @struct PolynomialInterp 低于7个点的多项式插值
	 * @detail 因为点数太多会出现扭摆的现象（曲线不稳定），因此不建议使用多余7个点的插值
	 * @author lyc
	 * @date 2020.2.22
	 */
	template<class TVert, class ftype, int dim = 3>
	class PolynomialCurve
	{
		//using ftype = float;
		//static const int dim = 3;
		//using TVert = std::array<ftype, 3>;
	private:
		std::vector<TVert>coefficient;
	public: //插值算法
		/*
		 * @brief 多项式插值曲线
		 * @detail 适用于6个点一下的插值曲线(否则出现扭摆现象)
		 * @author lyc
		 * @date 2020.2.22
		 */
		bool interpPolynome(const std::vector<TVert>& pts, ParamType param = ParamType::idChordLength)
		{
			assert(pts.size() <= 10 && "pts size must Less than 11");
			std::vector<ftype> uparam;
			getCurveParam(pts, uparam, param);
			auto sparam = uparam;
			std::sort(sparam.begin(), sparam.end());
			sparam.erase(std::unique(sparam.begin(), sparam.end()), sparam.end());
			if (sparam.size() != uparam.size()) //具有重复点
				return false;
			int npt = static_cast<int>(pts.size());
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> mat;
			mat.resize(npt, npt);
			for (int i = 0; i < npt; ++i) {
				for (int j = 0; j < npt; ++j) {
					mat(i, j) = static_cast<ftype>(pow(uparam[i], j));
				}
			}
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> pst;
			pst.resize(npt, dim);
			for (int i = 0; i < npt; ++i) {
				for (int j = 0; j < dim; ++j) {
					pst(i, j) = pts[i][j];
				}
			}
			auto cmat = mat.transpose()*mat;
			auto cpst = mat.transpose()*pst;
			auto bp = cmat.ldlt();
			if (bp.info() != Eigen::Success) //矩阵求解不成功
				return false;
			auto cb = bp.solve(cpst);
			coefficient.clear();
			coefficient.resize(npt);
			for (int i = 0; i < npt; ++i) {
				for (int j = 0; j < dim; ++j) {
					coefficient[i][j] = cb(i, j);
				}
			}
			return true;
		}
		TVert eval(ftype u) {
			using namespace BasicTypeCore;
			assert((u >= 0 && u <= 1) && "u must be [0,1]");
			int ncoeff = static_cast<int>(coefficient.size());
			TVert pt;
			for (int i = 0; i < dim; ++i)
				pt[i] = 0;
			for (int i = 0; i < ncoeff; ++i) {
				TVert p = AlgoUtil<TVert, ftype, dim>::muilt(static_cast<ftype>(pow(u, i)), coefficient[i]);
				pt = AlgoUtil<TVert, ftype, dim>::add(pt, p);
			}
			return pt;
		}
	};
	template<class TVert, class ftype, int dim = 3>
	class SegmentCurveInterpolation
	{
		//using ftype = float;
		//const static int dim = 3;
		//using TVert = std::array<ftype, dim>;
		using nanokdtree = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<ftype, kdtreenode<TVert, ftype, dim> >,
			kdtreenode<TVert, ftype, dim>, dim>;
	private:
		nanokdtree* kdtree = nullptr;
		kdtreenode<TVert, ftype, dim> kdt;
	private:
		void creatSearch()
		{
			kdt.mpts.clear();
			kdt.param.clear();
			AlgoUtil<TVert, ftype, dim>::getAverage(500, kdt.param);
			for (auto fp : kdt.param) {
				kdt.mpts.push_back(eval(fp));
			}
			if (kdtree != nullptr)
				delete kdtree;
			kdtree = new nanokdtree(dim, kdt);
			kdtree->buildIndex();
		}
	private:
		std::vector<TVert> mpts;
		std::vector<ftype> mweight;
		std::vector<TVert> mnormal;
		std::vector<ftype> muparam;
		CurveType mtype;
	private:
		std::vector<weightHermiteCubic<TVert, ftype, dim>> hermite;
		std::vector<CubeicCurve<TVert, ftype, dim>> cubeic;
		void startEvalCurve(CurveType ct = CurveType::idCubeicCurve) {
			assert((ct == CurveType::idCubeicCurve || ct == CurveType::idHermiteCubic) && "segment only be cuebeic or hermite!");
			int npt = static_cast<int>(mpts.size() - 1);
			mtype = ct;
			if (ct == CurveType::idHermiteCubic) {
				hermite.clear();
				hermite.reserve(npt);
				for (int i = 0; i < npt; ++i) {
					weightHermiteCubic<TVert, ftype, dim> whc;
					whc.interpolation(mpts[i], mpts[i + 1], mnormal[i], mnormal[i + 1]);
					hermite.push_back(whc);
				}
			}
			else if (ct == CurveType::idCubeicCurve) {
				cubeic.clear();
				cubeic.reserve(npt);
				for (int i = 0; i < npt; ++i) {
					CubeicCurve<TVert, ftype, dim> cc;
					cc.interpolation(mpts[i], mpts[i + 1], mnormal[i], mnormal[i + 1]);
					cubeic.push_back(cc);
				}
			}
			creatSearch();
		}

	public:
		SegmentCurveInterpolation():kdtree(nullptr){}
		~SegmentCurveInterpolation() {}
	public:
		/*
		 * @brief 分段曲线C1加权差值插值
		 * @author lyc
		 * @date 2020.2.22
		 */
		void intelpolation(
			const std::vector<TVert>& pts,
			const std::vector<TVert>& normal,
			ParamType param = ParamType::idChordLength, CurveType ct = CurveType::idCubeicCurve) {
			assert(pts.size() == normal.size() && "pts size must be equal to the size of normal!");
			getCurveParam<TVert, ftype>(pts, muparam, param);
			mpts = pts;
			mnormal = normal;
			mweight.clear();
			mweight.reserve(pts.size());
			for (auto& ip : pts)
				mweight.push_back(1);
			startEvalCurve(ct);
		}

		/*
		 * @brief 分段曲线C1加权差值插值
		 * @author lyc
		 * @date 2020.2.22
		 */
		void intelpolation(
			const std::vector<TVert>& pts,
			const std::vector<ftype> weight,
			const std::vector<TVert>& normal, ParamType param = ParamType::idChordLength, CurveType ct = CurveType::idCubeicCurve) {
			assert((pts.size() == normal.size() && weight.size() == normal.size()) && "the size of pts,weight and normal must be same!");
			getCurveParam<TVert, ftype>(pts, muparam, param);
			mpts = pts;
			mnormal = normal;
			mweight = weight;
			startEvalCurve(ct);
		}
		/*
		 * @brief 分段曲线C1加权差值插值
		 * @author lyc
		 * @date 2020.2.22
		 */
		void intelpolation(const std::vector<TVert>& pts,
			DerivativeType dt = DerivativeType::idBessel,
			ParamType param = ParamType::idChordLength, CurveType ct = CurveType::idCubeicCurve) {
			mpts = pts;
			mweight.clear();
			int npt = static_cast<int>(pts.size());
			mweight.reserve(npt);
			mnormal.clear();
			mnormal.reserve(npt);
			getCurveParam<TVert, ftype, dim>(pts, muparam, param);
			for (int i = 0; i < npt; ++i) {
				mweight.push_back(1);
				mnormal.push_back(getCurveDeriative<TVert, ftype, dim>(pts, muparam, i, dt));
			}
			startEvalCurve(ct);
		}

	public:
		TVert eval(ftype u) {
			assert((u >= 0 && u <= 1) && "the value of u must be [0,1]");
			int ip = 0;
			int npt = static_cast<int>(mpts.size() - 1);
			for (int i = 0; i < npt; ++i) {
				if (u >= muparam[i] && u <= muparam[i + 1]) {
					ip = i;
					break;
				}
			}
			switch (mtype)
			{
			case FussenAlgo::CurveType::idHermiteCubic:
				return hermite[ip].evalArbit(u, muparam[ip], muparam[ip + 1], mweight[ip], mweight[ip + 1]);
			case FussenAlgo::CurveType::idCubeicCurve:
				return cubeic[ip].eval(u, muparam[ip], muparam[ip + 1]);
			}
		}
		//反算参数
		ftype reval(const TVert& pt) {
			if (kdtree == nullptr)
				return -1;
			size_t index;
			ftype dist;
			std::vector<ftype> data(dim, static_cast<ftype>(0.0f));
			for (int i = 0; i < dim; ++i) {
				data[i] = pt[i];
			}
			kdtree->knnSearch(data.data(), 1, &index, &dist);
			return kdt.param[index];
		}
		//寻找离pt最近的n个曲线点
		void reval(const TVert& pt, int n, std::vector<ftype>& param) {
			if (kdtree == nullptr)
				return;
			std::vector<size_t> indexs(n, 0);
			std::vector<ftype> dists(n, 0);
			std::vector<ftype> data(dim, static_cast<ftype>(0.0f));
			for (int i = 0; i < dim; ++i) {
				data[i] = pt[i];
			}
			kdtree->knnSearch(data.data(), n, indexs.data(), dists.data());
			param.clear();
			param.reserve(n);
			for (auto fp : indexs) {
				param.push_back(kdt.param[fp]);
			}
		}

		//因为标准借口不可靠
		TVert curveTanght(ftype u) {
			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			TVert p;
			if (u <= static_cast<ftype>(0.999)) {
				p = AlgoUtil<TVert, ftype, dim>::sub(eval(u + static_cast<ftype>(0.001)), eval(u));
				p = AlgoUtil<TVert, ftype, dim>::normal(p);
			}
			else {
				p = AlgoUtil<TVert, ftype, dim>::sub(eval(u), eval(u - static_cast<ftype>(0.001)));
				p = AlgoUtil<TVert, ftype, dim>::normal(p);
			}
			return p;
		}

	};


	template<class TVert, class ftype, int dim = 3>
	class BSpline {
		//using ftype = float;
		//const static int dim = 3;
		//using TVert = std::array<ftype, dim>;
		using SplineVert = Eigen::Spline<ftype, dim + 1>;
		using nanokdtree = nanoflann::KDTreeSingleIndexAdaptor<
			nanoflann::L2_Simple_Adaptor<ftype, kdtreenode<TVert, ftype, dim> >,
			kdtreenode<TVert, ftype, dim>, dim>;
	private:
		SplineVert spline;
	private:
		std::vector<TVert> control;
		std::vector<ftype> knote;
		size_t mdegree;
		bool misClosed; //曲线是否封闭
		nanokdtree* kdtree = nullptr;
		kdtreenode<TVert, ftype, dim> kdt;
	private:
		void creatSearch()
		{
			kdt.mpts.clear();
			kdt.param.clear();
			AlgoUtil<TVert, ftype, dim>::getAverage(200, kdt.param);
			for (auto fp : kdt.param) {
				kdt.mpts.push_back(eval(fp));
			}
			if (kdtree != nullptr)
				delete kdtree;
			kdtree = new nanokdtree(dim, kdt);
			kdtree->buildIndex();
		}
	public:
		BSpline() :kdtree(nullptr) {}
		~BSpline() {
			if (kdtree)
				delete kdtree;
			kdtree = nullptr;
		}
	public:
		std::vector<TVert> getControl()const { return control; }
		std::vector<ftype> getKnote()const { return knote; }
		size_t degree()const { return mdegree; }
		bool isClosed()const { return misClosed; }
	public:
		//曲线的值
		TVert eval(ftype u) {
			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			auto p = spline(u);
			TVert cp;
			cp[0] = p(1);
			cp[1] = p(2);
			cp[2] = p(3);
			return cp;
		}
		//曲线的微商
		TVert evaldiff(ftype u, int n) {

			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			assert(n >= 0 && n <= mdegree && "the value of n must be [0,degree]");
			auto p = spline.derivatives(u, n);
			TVert cp;
			cp[0] = p(1);
			cp[1] = p(2);
			cp[2] = p(3);
			cp = AlgoUtil<TVert, ftype, dim>::normal(cp);
			return cp;
		}
		//因为标准借口不可靠
		TVert curveTanght(ftype u) {
			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			TVert p;
			if (u<=static_cast<ftype>(0.999)) {
				p = AlgoUtil<TVert, ftype, dim>::sub(eval(u + static_cast<ftype>(0.001)),eval(u));
				p=AlgoUtil<TVert, ftype, dim>::normal(p);
			}
			else {
				p = AlgoUtil<TVert, ftype, dim>::sub(eval(u), eval(u- static_cast<ftype>(0.001)));
				p=AlgoUtil<TVert, ftype, dim>::normal(p);
			}
			return p;
		}

	public:
		//反算参数
		ftype reval(const TVert& pt) {
			if (kdtree == nullptr)
				return -1;
			size_t index;
			ftype dist;
			std::vector<ftype> data(dim, static_cast<ftype>(0.0f));
			for (int i = 0; i < dim; ++i) {
				data[i] = pt[i];
			}
			kdtree->knnSearch(data.data(), 1, &index, &dist);
			return kdt.param[index];
		}
		//寻找离pt最近的n个曲线点
		void reval(const TVert& pt, int n, std::vector<ftype>& param) {
			if (kdtree == nullptr)
				return;
			std::vector<size_t> indexs(n, 0);
			std::vector<ftype> dists(n, 0);
			std::vector<ftype> data(dim, static_cast<ftype>(0.0f));
			for (int i = 0; i < dim; ++i) {
				data[i] = pt[i];
			}
			kdtree->knnSearch(data.data(), n, indexs.data(), dists.data());
			param.clear();
			param.reserve(n);
			for (auto fp : indexs) {
				param.push_back(kdt.param[fp]);
			}
		}
	public:
		void  splitCurve(ftype umin, ftype umax, BSpline<TVert, ftype, dim>& split) {
			std::vector<ftype> param;
			distPoints<ftype>::getAverage(50, param);
			std::vector<TVert> pts;
			for (int i = 0; i < 50; ++i) {
				pts.push_back(eval(umin+(umax-umin)*param[i]));
			}
			split.interpolationSpline(pts);
		}

	public: //插值曲线
		/*
		 * @brief b样条插值曲线
		 * @param[in] pts 曲线上的采样点
		 * @param[in] degree degree次B样条
		 * @param[in] param 曲线参数化类型
		 */
		bool interpolationSpline(const std::vector<TVert>& ori_pts,
			size_t degree = 3, ParamType param = ParamType::idChordLength,
			bool isclosed = false, ftype alpha = static_cast<ftype>(0.8)) {
			std::vector<TVert> pts;
			AlgoUtil<TVert, ftype, dim>::repeat(ori_pts, pts);
			misClosed = isclosed;
			if (misClosed&&AlgoUtil<TVert, ftype, dim>::length(pts[0], pts[pts.size() - 1]) > BasicTypeCore::epslion)
				pts.push_back(pts[0]);
			if (pts.size() == 2)
				degree = 1;
			if (pts.size() < degree)
				return false;
			std::vector<ftype> uparam;
			int npt = static_cast<int>(pts.size());
			getCurveParam<TVert, ftype, dim>(pts, uparam, param, degree, alpha);
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> points;
			points.resize(dim + 1, npt);
			for (int i = 0; i < npt; ++i) {
				points(0, i) = uparam[i];
				for (int j = 1; j <= dim; ++j) {
					points(j, i) = pts[i][j - 1];
				}
			}
			spline = Eigen::SplineFitting<SplineVert>::Interpolate(points, degree);
			auto ckonte = spline.knots();
			int ncols = static_cast<int>(ckonte.cols());
			knote.clear();
			knote.reserve(ncols);
			for (int i = 0; i < ncols; ++i) {
				knote.push_back(ckonte(0, i));
			}
			auto ncontrol = spline.ctrls();
			control.clear();
			int nc = static_cast<int>(ncontrol.cols());
			control.reserve(nc);
			for (int i = 0; i < nc; ++i) {
				TVert vert;
				for (int j = 0; j < dim; ++j) {
					vert[j] = ncontrol(j + 1, i);
				}
				control.emplace_back(vert);
			}
			mdegree = static_cast<int>(spline.degree());
			creatSearch();
			return true;
		}
		/*
		 * @brief 带法（部分法向量也行）向量约束的样条插值曲线
		 * @param[in] pts 曲线上的采样点
		 * @param[in] pnormal 曲线的法向量
		 * @param[in] index  法向量对应的点的下标
		 * @param[in] degree degree次B样条
		 * @param[in] param 曲线参数化类型
		 */
		bool interpolationSpline(
			const std::vector<TVert>& ori_pts,
			const std::vector<TVert>& cpnormal, const std::vector<int>& cindex,
			ftype alpha = static_cast<ftype>(0.8),
			size_t degree = 3, bool isclose = false,
			ParamType param = ParamType::idChordLength)
		{
			std::vector<TVert> pts = ori_pts, pnormal = cpnormal;
			std::vector<int> index = cindex;
			misClosed = isclose;
			if (misClosed&&AlgoUtil<TVert, ftype, dim>::length(pts[0], pts[pts.size() - 1]) > BasicTypeCore::epslion) {
				pts.push_back(pts[0]);
			}
			if (pts.size() < degree || pts.size() < 2 || pnormal.size() != index.size())
				return false;
			std::vector<ftype> uparam;
			int npt = static_cast<int>(pts.size());
			getCurveParam<TVert, ftype, dim>(pts, uparam, param, degree, alpha);
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> points;
			points.resize(dim + 1, npt);
			for (int i = 0; i < npt; ++i) {
				points(0, i) = uparam[i];
				for (int j = 1; j <= dim; ++j) {
					points(j, i) = pts[i][j - 1];
				}
			}
			int nnormal = static_cast<int>(pnormal.size());
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> normal;
			normal.resize(dim + 1, nnormal);
			Eigen::VectorXi pindex(nnormal);
			for (int i = 0; i < nnormal; ++i) {
				normal(0, i) = uparam[index[i]];
				for (int j = 1; j <= dim; ++j) {
					normal(j, i) = pnormal[i][j - 1];
				}
				pindex(i) = index[i];
			}
			spline = Eigen::SplineFitting<SplineVert>::
				InterpolateWithDerivatives(points, normal, pindex, degree);
			auto ckonte = spline.knots();
			int ncols = static_cast<int>(ckonte.cols());
			knote.clear();
			knote.reserve(ncols);
			for (int i = 0; i < ncols; ++i) {
				knote.push_back(ckonte(0, i));
			}
			auto ncontrol = spline.ctrls();
			control.clear();
			int nc = static_cast<int>(ncontrol.cols());
			control.reserve(nc);
			for (int i = 0; i < nc; ++i)
			{
				TVert vert;
				for (int j = 0; j < dim; ++j) {
					vert[j] = ncontrol(j + 1, i);
				}
				control.emplace_back(vert);
			}
			mdegree = static_cast<int>(spline.degree());
			creatSearch();
			return true;
		}
	public:
		/*
		 * @brief 最小二乘逼近曲线（曲线不经过采样点，拟合一条最小二乘一样下的样条曲线）
		 * @param[in] pts 原始采样点，如果相邻点重合，会自动去重
		 * @param[in] mcontrol 拟合样条曲线控制点数，样条曲线控制点越多，越接近原始曲线，但是一般不超过n+degree+2个
		 * @param[in] degree 样条曲线阶数，一般取3，样条曲线次数太高会出现扭摆现象
		 * @param[in] alpha 参数因子
		 * @param[in] param 样条曲线的参数化方式，一般选用积累弦长参数化
		 */
		bool ApproximationSpline(const std::vector<TVert>& ori_pts, int mcontrol,
			size_t degree = 3, ftype alpha = static_cast<ftype>(0.8),
			ParamType param = ParamType::idChordLength)
		{

			using CAlgo = AlgoUtil<TVert, ftype, dim>;
			using Matrix = Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic>;
			std::vector<TVert> qts;
			CAlgo::repeat(ori_pts, qts);
			assert(qts.size() >= mcontrol + 1 && mcontrol + 1 >= degree && degree >= 1);
			mdegree = degree;
			int m = static_cast<int>(qts.size() - 1);
			std::vector<ftype> uparam;
			getCurveParam<TVert, ftype, dim>(qts, uparam, param, degree, alpha);
			Matrix vect;
			int n = mcontrol - 1;
			bool isok = true;
			do
			{
				if (!isok) {
					int nk = (m - mdegree - 1) / mdegree;
					if (n - nk > mdegree&&n < m)
						n -= nk;
					else if (n - 1 > mdegree)
						--n;
					else
						return false;
				}
				//计算节点向量
				knote.clear();
				knote.resize(n + mdegree + 2);
				for (int j = 0; j <= mdegree; ++j)
					knote[j] = static_cast<ftype>(0.0);
				ftype d = static_cast<ftype>(m + 1) / static_cast<ftype>(n - mdegree + 1);
				//这个理论存在争议，需要进一步探索
				for (int j = 1; j <= n - mdegree; ++j) {
					int i = static_cast<int>(j*d);
					ftype alpha = j * d - i;
					knote[mdegree + j] = (1 - alpha)*uparam[i - 1] + alpha * uparam[i];
				}
				for (int j = 1; j <= mdegree + 1; ++j)
					knote[n + j] = 1;
				//计算控制点
				control.clear();
				control.resize(n + 1);
				control[0] = qts[0];
				control[n] = qts[m];
				//计算control[i],i=1,2,---,mcontrl-1
				Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> matrix;
				matrix.resize(m - 1, n - 1);
				matrix.setZero();
				//这一段代码可以再进一步优化
				for (int i = 1; i <= m - 1; ++i) {
					for (int j = 1; j <= n - 1; ++j) {
						matrix(i - 1, j - 1) = CAlgo::bsplineOneBasis(j, mdegree, knote, uparam[i]);
					}
				}
				Matrix cmat = matrix.transpose()*matrix;
				//计算Rk
				std::vector<TVert> dR;
				dR.reserve(m - 1);
				for (int i = 1; i <= m - 1; ++i) {
					TVert r0 = qts[i];
					TVert r1 = CAlgo::muilt(
						CAlgo::bsplineOneBasis(0, mdegree, knote, uparam[i]), qts[0]);
					TVert r2 = CAlgo::muilt(
						CAlgo::bsplineOneBasis(n, mdegree, knote, uparam[i]), qts[m]);
					TVert r = CAlgo::sub(r0, r1);
					r = CAlgo::sub(r, r2);
					dR.emplace_back(r);
				}
				std::vector<TVert> R;
				R.reserve(n - 1);
				for (int i = 1; i <= n - 1; ++i) {
					TVert p = CAlgo::muilt(CAlgo::bsplineOneBasis(i, mdegree, knote, uparam[1]), dR[0]);
					for (int j = 2; j <= m - 1; ++j) {
						TVert cp = CAlgo::muilt(CAlgo::bsplineOneBasis(i, mdegree, knote, uparam[j]), dR[j - 1]);
						p = CAlgo::add(p, cp);
					}
					R.push_back(p);
				}
				Matrix vec;
				vec.resize(n - 1, dim + 1);
				for (int i = 0; i < n - 1; ++i) {
					vec(i, 0) = uparam[i];
					for (int j = 1; j <= dim; ++j)
						vec(i, j) = R[i][j - 1];
				}
				auto solver = cmat.householderQr();
				vect = solver.solve(vec);
				isok = vec.isApprox(cmat*vect);
			} while (!isok);
			for (int d = 1; d < dim + 1; ++d) {
				for (int i = 1; i < n; ++i)
					control[i][d - 1] = vect(i - 1, d);
			}
			typename Eigen::SplineTraits<SplineVert>::ControlPointVectorType ctt;
			typename Eigen::SplineTraits<SplineVert>::KnotVectorType kvt;
			ctt.resize(dim + 1, n + 1);
			for (int i = 0; i < n + 1; ++i) {
				ctt(0, i) = 1;
				for (int j = 1; j <= dim; ++j) {
					ctt(j, i) = control[i][j - 1];
				}
			}
			kvt.resize(1, n + mdegree + 2);
			for (int i = 0; i < n + mdegree + 2; ++i) {
				kvt(0, i) = knote[i];
			}
			spline = SplineVert(kvt, ctt);
			creatSearch();
			return true;
		}

		/*
		 * @brief 最小二乘逼近曲线（曲线不经过采样点，拟合一条最小二乘一样下的样条曲线）
		 * @param[in] pts 原始采样点，如果相邻点重合，会自动去重
		 * @param[in] mcontrol 拟合样条曲线控制点数，样条曲线控制点越多，越接近原始曲线，但是一般不超过n+degree+2个
		 * @param[in] degree 样条曲线阶数，一般取3，样条曲线次数太高会出现扭摆现象
		 * @param[in] alpha 参数因子
		 * @param[in] param 样条曲线的参数化方式，一般选用积累弦长参数化
		 */
		bool ApproximationSplineColse(const std::vector<TVert>& ori_pts, int mcontrol,
			size_t degree = 3, ftype alpha = static_cast<ftype>(0.8),
			ParamType param = ParamType::idChordLength) 
		{
			if (ApproximationSpline(ori_pts, mcontrol, degree, alpha, param)) {
				std::vector<TVert> pts;
				int nsize = (int)ori_pts.size();
				pts.reserve(nsize);
				std::vector<ftype> paramt;
				distPoints<ftype>::getAverage(nsize, paramt);
				for (int i = 0; i < nsize; ++i) {
					pts.emplace_back(eval(paramt[i]));
				}
				return interpolationSpline(pts, degree, param, true, alpha);
			}
			else
				return false;
		}


	};
	template<class TVert, class ftype, int dim = 3>
	class Nurbs {
		//using ftype = float;
		//const static int dim = 3;
		//using TVert = std::array<ftype, dim>;
		using algok = AlgoUtil<TVert, ftype, dim>;
	private:
		std::vector<ftype> knotes;
		std::vector<TVert> control;
		std::vector<ftype> weights;
		size_t degree;
	public:
		Nurbs() {};
		Nurbs(const BSpline<TVert, ftype, dim>& spline, std::vector<ftype>& cweights) {
			knotes = spline.getKnote();
			degree = spline.degree();
			control = spline.getControl();
			weights = cweights;
		}
	public:
		TVert eval(ftype u) {
			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			return algok::curvePoint(degree, knotes, control, weights, u);
		}
		TVert curveTanght(ftype u) {
			assert(u >= 0 && u <= 1 && "the value of u must be [0,1]!");
			return algok::curveTangent(degree, knotes, control, weights, u);
		}
	};

	//线统一输出asc格式
	namespace IO {

		template<class T>
		void debug(std::string filaname, const std::vector<T>& pts, int r = 255, int g = 0, int b = 0) {
			std::ofstream fcout(filaname);
			for (auto ip : pts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			fcout.close();
		}

		template<class TVert, class ftype, int dim = 3>
		void debug(std::string filaname, BSpline<TVert, ftype, dim>& spline,
			int r = 255, int g = 0, int b = 0)
		{
			int npts = 300;
			std::vector<ftype> param;
			distPoints<ftype>::getAverage(npts, param);
			std::ofstream fcout(filaname);
			std::vector<TVert> cpts;
			for (auto ip : param)
				cpts.emplace_back(spline.eval(ip));
			for (auto ip : cpts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			fcout.close();
		}


		template<class TVert, class ftype, int dim = 3>
		void debugNorm(std::string filaname, BSpline<TVert, ftype, dim>& spline,
			int r = 255, int g = 0, int b = 0,int np=20,int nr=120,int ng=0,int nb=0)
		{
			int npts = 300;
			std::vector<ftype> param,nparam;
			distPoints<ftype>::getAverage(npts, param);
			std::ofstream fcout(filaname);
			std::vector<TVert> cpts;
			for (auto ip : param)
				cpts.emplace_back(spline.eval(ip));
			for (auto ip : cpts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			distPoints<ftype>::getAverage(np, nparam);
			for (int i = 0; i < np; ++i) {
				std::vector<TVert> cpt;
				TVert tt = spline.eval(nparam[i]);
				cpt.push_back(tt);
				TVert norm = spline.curveTanght(nparam[i]);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				BSpline<TVert, ftype, dim> nspline;
				nspline.interpolationSpline(cpt);
				for (auto ip : param) {
					TVert tt = nspline.eval(ip);
					fcout << tt[0] << " " << tt[1] << " " << tt[2] << " " << nr << " " << ng << " " << nb << std::endl;
				}
			}
			fcout.close();
		}

		template<class TVert, class ftype, int dim = 3>
		void debugNorm(std::string filaname, CubeicCurve<TVert, ftype, dim>& cc,
			int r = 255, int g = 0, int b = 0, int np = 20, int nr = 120, int ng = 0, int nb = 0)
		{
			int npts = 300;
			std::vector<ftype> param, nparam;
			distPoints<ftype>::getAverage(npts, param);
			std::ofstream fcout(filaname);
			std::vector<TVert> cpts;
			for (auto ip : param)
				cpts.emplace_back(cc.eval(ip));
			for (auto ip : cpts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			distPoints<ftype>::getAverage(np, nparam);
			for (int i = 0; i < np; ++i) {
				std::vector<TVert> cpt;
				TVert tt = cc.eval(nparam[i]);
				cpt.push_back(tt);
				TVert norm = cc.curveTanght(nparam[i]);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				tt = AlgoUtil<TVert, ftype, dim>::add(tt, norm);
				cpt.push_back(tt);
				BSpline<TVert, ftype, dim> nspline;
				nspline.interpolationSpline(cpt);
				for (auto ip : param) {
					TVert tt = nspline.eval(ip);
					fcout << tt[0] << " " << tt[1] << " " << tt[2] << " " << nr << " " << ng << " " << nb << std::endl;
				}
			}
			fcout.close();
		}

		template<class TVert, class ftype, int dim = 3>
		void debug(std::string filaname, Nurbs<TVert, ftype, dim>& nurbs,
			int r = 255, int g = 0, int b = 0)
		{
			int npts = 300;
			std::vector<ftype> param;
			distPoints<ftype>::getAverage(npts, param);
			std::ofstream fcout(filaname);
			std::vector<TVert> cpts;
			for (auto ip : param)
				cpts.emplace_back(nurbs.eval(ip));
			for (auto ip : cpts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			fcout.close();
		}

		template<class TVert, class ftype, int dim = 3>
		void debug(std::string filaname, CubeicCurve<TVert, ftype, dim>& spline,
			int r = 255, int g = 0, int b = 0)
		{
			int npts = 200;
			std::vector<ftype> param;
			distPoints<ftype>::getAverage(npts, param);
			std::ofstream fcout(filaname);
			std::vector<TVert> cpts;
			for (auto ip : param)
				cpts.emplace_back(spline.eval(ip));
			for (auto ip : cpts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			fcout.close();
		}

		template<class TVert, class ftype, int dim = 3>
		void debug(std::string filaname, SegmentCurveInterpolation<TVert, ftype, dim>& spline,
			int r = 255, int g = 0, int b = 0)
		{
			int npts = 200;
			std::vector<ftype> param;
			distPoints<ftype>::getAverage(npts, param);
			std::ofstream fcout(filaname);
			std::vector<TVert> cpts;
			for (auto ip : param)
				cpts.emplace_back(spline.eval(ip));
			for (auto ip : cpts) {
				fcout << ip[0] << " " << ip[1] << " " << ip[2] << " " << r << " " << g << " " << b << std::endl;
			}
			fcout.close();
		}
	}
}
#pragma warning(pop)