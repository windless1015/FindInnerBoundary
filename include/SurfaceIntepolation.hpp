#pragma once
#ifndef SURFACE_INTEPOLTION_H
#define SURFACE_INTEPOLTION_H
#include <iostream>
#include <array>
#include <vector>
#include <assert.h>
#include <Eigen/Eigen>
using namespace std;
/*
 * @brief 一个仅头文件的曲面插值和逼近算法
 * @detail 
 *   ---Eigen--提供矩阵计算的算法（开源免费的）
 * @date 2020.3.16
 * @detail 
 * @modify
 *   ---2020.3.16---三角域超限插值算法
 */

namespace FussenAlgo
{
	namespace SurfaceAlgo
	{
		/*
		 * @brief 提供数值修正与检测算法
		 * @date 2020.3.16
		 */
		template <class Real>
		struct utils {
			static Real realThreshold() { return Real(1e-13); }
			static Real realPrecision() { return Real(1e-10); }
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
				else if (rl<=0 && rl>-realPrecision())
					return -realPrecision();
				else
					return rl;
			}
			//检查数值有效性
			static bool checkEval(Real rl) {
				return !std::isnan<Real>(rl)&&!std::isinf<Real>(rl);
			}
			//检查数值的零值
			static bool isZero(Real rl){
				return rl >= -realPrecision() && rl <= realPrecision();
			}
		};

		template<class TVert,class ftype>
		class CmpVec
		{
		public:
			explicit CmpVec(ftype _eps = FLT_MIN) : eps_(_eps) {}
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
		private:
			ftype eps_;
		};

		/*
		 * @brief 提供类型基础算法
		 * @date 2020.3.16
		 * @author lyc
		 */
		template<class TVert,class ftype=float,int dim=3>
		class AlgoUtil {
		public:
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
				ftype len = length(p);
				len=utils<ftype>::fixPos(len);
				return divd(p, len);
			}
		};
	}

	template<class TVert,class ftype>
	class TriangleTransfinite
	{
		//using ftype = float;
		//using TVert = std::array<ftype, 3>;
	protected:
		using AlgoUtil = SurfaceAlgo::AlgoUtil<TVert, ftype, 3>;
		/* @brief 计算三角形的面积 
		 *          the
		 *          /\
		 *         /  \
		 *        /    \
		 *     one------two
		 * @detail Area=(two-one)x(the-one)*0.5
		 */
		static ftype clacArea(const TVert& one, const TVert& two, const TVert& the) {
			TVert ot = AlgoUtil::sub(two, one);
			TVert oh = AlgoUtil::sub(the, one);
			TVert al = AlgoUtil::cross(ot, oh);
			return AlgoUtil::length(al)*static_cast<ftype>(0.5);
		}
	public:
		//点边插值算法

		//移动插值算法

		//面面插值算法

	};

	template<typename ftype>
	struct RbfPow3
	{
		 inline ftype f(const ftype& x) { return  x*x*x; }
		 inline ftype df(const ftype& x) { return ftype(3) * x * x; }
		 inline ftype ddf(const ftype& x) { return ftype(6) * x; }
	};

	template<class TVert, class ftype,class rbffun>
	class rbfImplicitSurface {
	private:
		std::vector<TVert> napts;
		std::vector<ftype> alpha;
	public:
		bool interpolation(std::vector<TVert> pts,
			std::vector<TVert>& normal) {
			assert(pts.size() == normal.size() && "the size of pts and normal must be eauql!");
			int npts = static_cast<int>(pts.size());
			rbffun fun;
			int cpts = 3 * npts+4;
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> cmat;
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> vec;
			napts.clear();
			napts.reserve(cpts);
			std::vector<ftype> vals;
			vals.reserve(cpts);
			for (int i = 0; i < npts; ++i) {
				napts.emplace_back(pts[i]);
				napts.emplace_back(AlgoUtil<TVert,ftype,3>::add(pts[i], normal[i]));
				napts.emplace_back(AlgoUtil<TVert,ftype,3>::sub(pts[i], normal[i]));
				vals.push_back(0);
				vals.push_back(1);
				vals.push_back(-1);
			}
			cmat.resize(cpts, cpts);
			vec.resize(cpts, 1);
			cmat.setZero();
			vec.setZero();
			int ccnp = static_cast<int>(napts.size());
#pragma omp parallel for
			for (int i = 0; i < ccnp; ++i) {
				TVert it = napts[i];
				for (int j = i; j < ccnp; ++j) {
					TVert jt = napts[j];
					ftype cr = AlgoUtil<TVert, ftype, 3>::length(it, jt);
					cmat(j,i)=cmat(i, j) = fun.f(cr);
				}
				cmat(ccnp, i) = cmat(i, ccnp) = 1;
				cmat(ccnp + 1, i) = cmat(i, ccnp + 1) = it[0];
				cmat(ccnp + 2, i) = cmat(i, ccnp + 2) = it[1];
				cmat(ccnp + 3, i) = cmat(i, ccnp + 3) = it[2];
				vec(i, 0) = vals[i];
			}
			alpha.clear();
			alpha.reserve(cpts);
			Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic> cval;
			auto solver = cmat.householderQr();
			cval = solver.solve(vec);
			for (int i = 0; i < cpts; ++i) {
				alpha.push_back(cval(i,0));
			}
			return true;
		}
	public:
		ftype eval(const TVert& pt) {
			ftype ft = static_cast<ftype>(0);
			int npt = static_cast<int>(napts.size());
			rbffun fun;
			for (int i = 0; i < npt;++i) {
				ftype dl = AlgoUtil<TVert,ftype,3>::length(pt, napts[i]);
				ftype cl=fun.f(dl);
				ft += alpha[i] * cl;
			}
			ft += alpha[npt];
			ft += alpha[npt + 1] * pt[0];
			ft += alpha[npt + 2] * pt[1];
			ft += alpha[npt + 3] * pt[2];
			return ft;
		}

		TVert grad(const TVert& pt) {
			TVert p;
			p[0] = static_cast<ftype>(0);
			p[1] = static_cast<ftype>(0);
			p[2] = static_cast<ftype>(0);
			int npt = static_cast<int>(napts.size());
			for (int i = 0; i < npt; ++i) {
				ftype dl = AlgoUtil<TVert, ftype, 3>::length(pt, napts[i]);
				p[0] = p[0] + 3 * alpha[i] * (pt[0] - napts[i][0])*dl;
				p[1] = p[1] + 3 * alpha[i] * (pt[1] - napts[i][1])*dl;
				p[2] = p[2] + 3 * alpha[i] * (pt[2] - napts[i][2])*dl;
			}
			p[0] += alpha[npt + 1];
			p[1] += alpha[npt + 2];
			p[2] += alpha[npt + 3];
			return p;
		}

		TVert newton(const TVert& src)
		{
			TVert st = src;
			TVert et = src;

		}


	};


}



#endif