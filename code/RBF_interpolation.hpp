#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <vector>
#include <iostream>

using namespace std;
namespace FussenAlgo
{
	template<class ftype>
	class RBFCore
	{
	public:
		using Matrix = Eigen::Matrix<ftype, Eigen::Dynamic, Eigen::Dynamic>;
	public:
		/*
		 * @brief RBF径向基函数插值
		 * @param[in] x 采样点坐标,每一行是一个数据
		 * @param[in] y 采样点值,函数值是一列多行
		 * @param[in] RBFFunction 插值函数
		 * @param[in] const_num 插值函数参数
		 * @param[in] RBFSmooth 平滑参数 0是不平滑
		 */
		static Matrix rbfcreate( Matrix& x,Matrix& y, 
			Matrix RBFFunction(Matrix r, ftype const_num),
			ftype RBFConstant, ftype RBFSmooth)
		{
			int nXDim = x.cols();
			int nYDim = y.cols();
			int nXCount = x.rows();
			int nYCount = y.rows();
			if (nXCount != nYCount)
				std::cerr << "x and y should have the same number of rows" << std::endl;
			if (nYDim != 1)
				std::cerr << "y should be n by 1 vector" << std::endl;
			Matrix A;
			A = rbfAssemble(x, RBFFunction, RBFConstant, RBFSmooth);
			Matrix b = Matrix::Zero(y.rows() + nXDim + 1, 1);
			b.topRows(y.rows()) = y;
			Matrix rbfcoeff;
			rbfcoeff = A.lu().solve(b);
			return rbfcoeff;
		}
	public:
		/*
		 * @brief rbf插值
		 * @param[in] nodes 采样点坐标
		 * @param[in] rbfcoeff rbf参数
		 * @param[in] x 待插值点坐标
		 * @param[in] RBFFunction 插值函数
		 * @param[in] RBFConstant 插值函数参数
		 */
		static Matrix rbfinterp(Matrix& nodes, Matrix& rbfcoeff, Matrix& x,
			Matrix RBFFunction(Matrix r, ftype const_num),
			ftype RBFConstant) {
			int dim =(int) nodes.cols();
			int n = (int)nodes.rows();
			int dimPoints = x.cols();
			int nPoints = x.rows();
			if (dim != dimPoints)
				std::cerr << "x should have the same number of rows as an array used to create RBF interpolation" << std::endl;
			Matrix r = Matrix::Zero(nPoints, n);
			Matrix temp_A;
			for (int i = 0; i < nPoints; i++)
			{
				for (int j = 0; j < n; j++)
				{
					r(i, j) = (x.row(i) - nodes.row(j)).norm();
				}
			}
			temp_A = RBFFunction(r, RBFConstant);
			Matrix P(x.rows(), x.cols() + 1);
			P.leftCols(1) = Matrix::Ones(x.rows(), 1);
			P.rightCols(x.cols()) = x;
			Matrix A(nPoints, n + x.cols() + 1);
			A.topLeftCorner(temp_A.rows(), temp_A.cols()) = temp_A;
			A.topRightCorner(P.rows(), P.cols()) = P;
			Matrix f;
			f = A * rbfcoeff;
			return f;
		}
		static ftype rbfcheck(Matrix& x, Matrix& y, Matrix& rbfcoeff,
			Matrix RBFFunction(Matrix r, ftype const_num),
			ftype RBFConstant)
		{
			Matrix S;
			S = rbfinterp(x, rbfcoeff, x, RBFFunction, RBFConstant);
			ftype maxdiff;
			maxdiff = (S - y).cwiseAbs().maxCoeff();
			return maxdiff;
		}
	public:
		static Matrix rbfphi_linear(Matrix r, ftype const_num) {
			Matrix u(r.rows(), r.cols());
			u = r;
			return u;
		}
		static Matrix rbfphi_cubic(Matrix r, ftype const_num) {
			Matrix u(r.rows(), r.cols());
			u = r.array().pow(3);
			return u;
		}
		static Matrix rbfphi_gaussian(Matrix r, ftype const_num) {
			Matrix u(r.rows(), r.cols());
			u = (-0.5*r.array().square() / (const_num*const_num)).array().exp();
			return u;
		}
		static Matrix rbfphi_multiquadrics(Matrix r, ftype const_num) {
			Matrix u(r.rows(), r.cols());
			u = (1.0 + r.array().square() / (const_num*const_num)).array().sqrt();
			return u;
		}
		static Matrix rbfphi_thinplate(Matrix r, ftype const_num) {
			Matrix u(r.rows(), r.cols());
			u = (r.array().square()).cwiseProduct((r.array() + 1).log());
			return u;
		}
	protected:
		static Matrix rbfAssemble(Matrix& x, Matrix RBFFunction(Matrix r, ftype const_num),
			ftype RBFConstant, ftype RBFSmooth) {
			int dim = x.cols();//维度数
			int n = x.rows();//采样点数
			Matrix r = Matrix::Zero(n, n);
			Matrix temp_A;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < i; j++)
				{
					r(i, j) = (x.row(i) - x.row(j)).norm();
					r(j, i) = r(i, j);
				}
			}
			temp_A = RBFFunction(r, RBFConstant);
			//平滑:使得采样点不再精确等于采样值
			for (int i = 0; i < n; i++)
			{
				temp_A(i, i) = temp_A(i, i) - RBFSmooth;
			}
			Matrix P(x.rows(), x.cols() + 1);
			P.leftCols(1) = Matrix::Ones(n, 1);
			P.rightCols(x.cols()) = x;
			Matrix A = Matrix::Zero(temp_A.rows() + P.cols(), temp_A.cols() + P.cols());
			A.topLeftCorner(temp_A.rows(), temp_A.cols()) = temp_A;
			A.topRightCorner(P.rows(), P.cols()) = P;
			A.bottomLeftCorner(P.cols(), P.rows()) = P.transpose();
			return A;
		}
	};






	/// @brief fitting surface on a cloud point and evaluating the implicit surface
	/// @tparam _Scalar : a base type float, double etc.
	/// @tparam _Dim    : integer of the dimension of the ambient space
	/// (for a implicit surface == 3)
	/// @tparam Rbf     : the class of a the radial basis
	/// must implement float Rbf::f(float) float Rbf::df(float) float Rbf::ddf(float)
	/// (see hrbf_phi_funcs.h for an example)
	/// @note Please go see http://eigen.tuxfamily.org to use matrix and vectors
	/// types of this lib. The documentation is pretty good.
	template<typename _Scalar, int _Dim, typename Rbf>
	class HRBF_fit
	{
	public:
		typedef _Scalar Scalar;
		enum { Dim = _Dim };

		typedef Eigen::Matrix<Scalar, Dim, Dim>                       MatrixDD;
		typedef Eigen::Matrix<Scalar, Dim, 1>                         Vector;
		typedef Eigen::Matrix<Scalar, Dim, Eigen::Dynamic>            MatrixDX;
		typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXX;
		typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1>              VectorX;

		HRBF_fit() {}

		// --------------------------------------------------------------------------

		/// Compute surface interpolation given a set of points and normals.
		/// This solve the lineat system of equation to find alpha scalars and beta
		/// beta vectors stored in '_alphas' and '_betas' attributes.
		void hermite_fit(const std::vector<Vector>& points,
			const std::vector<Vector>& normals)
		{
			assert(points.size() == normals.size());

			int nb_points = points.size();
			int nb_hrbf_constraints = (Dim + 1)*nb_points;
			int nb_constraints = nb_hrbf_constraints;
			int nb_coeffs = (Dim + 1)*nb_points;

			_node_centers.resize(Dim, nb_points);
			_betas.resize(Dim, nb_points);
			_alphas.resize(nb_points);

			// Assemble the "design" and "value" matrix and vector
			MatrixXX  D(nb_constraints, nb_coeffs);
			VectorX   f(nb_constraints);
			VectorX   x(nb_coeffs);

			// copy the node centers
			for (int i = 0; i < nb_points; ++i)
				_node_centers.col(i) = points[i];

			for (int i = 0; i < nb_points; ++i)
			{
				Vector p = points[i];
				Vector n = normals[i];

				int io = (Dim + 1) * i;
				f(io) = 0;
				f.template segment<Dim>(io + 1) = n;

				for (int j = 0; j < nb_points; ++j)
				{
					int jo = (Dim + 1) * j;
					Vector diff = p - _node_centers.col(j);
					Scalar l = diff.norm();
					if (l == 0) {
						D.template block<Dim + 1, Dim + 1>(io, jo).setZero();
					}
					else {
						Scalar w = Rbf::f(l);
						Scalar dw_l = Rbf::df(l) / l;
						Scalar ddw = Rbf::ddf(l);
						Vector g = diff * dw_l;
						D(io, jo) = w;
						D.row(io).template segment<Dim>(jo + 1) = g.transpose();
						D.col(jo).template segment<Dim>(io + 1) = g;
						D.template block<Dim, Dim>(io + 1, jo + 1) = (ddw - dw_l) / (l*l) * (diff * diff.transpose());
						D.template block<Dim, Dim>(io + 1, jo + 1).diagonal().array() += dw_l;
					}
				}
			}

			x = D.lu().solve(f);
			Eigen::Map< Eigen::Matrix<Scalar, Dim + 1, Eigen::Dynamic> > mx(x.data(), Dim + 1, nb_points);

			_alphas = mx.row(0);
			_betas = mx.template bottomRows<Dim>();
		}

		// -------------------------------------------------------------------------

		/// Evaluate potential f() at position 'x'
		Scalar eval(const Vector& x) const
		{
			Scalar ret = 0;
			int nb_nodes = _node_centers.cols();
			Vector diff;
			Scalar l;
			for (int i = 0; i < nb_nodes; ++i)
			{
				diff = x - _node_centers.col(i);
				l = diff.norm();

				if (l > 0.00001f) {
					ret += _alphas(i) * Rbf::f(l) + _betas.col(i).dot(diff) * Rbf::df(l) / l;
				}
			}
			return ret;
		}

		// -------------------------------------------------------------------------

		/// Evaluate gradient nabla f() at position 'x'
		Vector grad(const Vector& x) const
		{
			Vector grad = Vector::Zero();
			int nb_nodes = _node_centers.cols();
			for (int i = 0; i < nb_nodes; ++i)
			{
				Vector beta = _betas.col(i);
				float  alpha = _alphas(i);
				Vector diff = x - _node_centers.col(i);

				Vector diffNormalized = diff;
				float l = diff.norm();

				if (l > 0.00001f)
				{
					diffNormalized.normalize();
					float dphi = Rbf::df(l);
					float ddphi = Rbf::ddf(l);

					float alpha_dphi = alpha * dphi;

					float bDotd_l = beta.dot(diff) / l;
					float squared_l = diff.squaredNorm();

					grad += alpha_dphi * diffNormalized;
					grad += bDotd_l * (ddphi * diffNormalized - diff * dphi / squared_l) + beta * dphi / l;
				}
			}
			return grad;
		}


		// -------------------------------------------------------------------------

		Scalar eval_n_grad(const Vector& x, Vector& grad) const
		{
			Scalar ret = 0;
			grad = Vector::Zero();
			int nb_nodes = _node_centers.cols();

			Vector beta;
			float  alpha;
			Vector diff;
			Vector diff_normalized;
			float l;
			float dphi;
			float bDotd_l;
			for (int i = 0; i < nb_nodes; ++i)
			{
				beta = _betas.col(i);
				alpha = _alphas(i);
				diff = x - _node_centers.col(i);

				diff_normalized = diff;
				l = diff.norm();

				if (l > 0.00001f)
				{
					diff_normalized.normalize();
					dphi = Rbf::df(l);
					bDotd_l = beta.dot(diff) / l;

					grad += alpha * dphi * diff_normalized +
						bDotd_l * (Rbf::ddf(l) * diff_normalized - diff * dphi / diff.squaredNorm()) +
						beta * dphi / l;

					ret += alpha * Rbf::f(l) + bDotd_l * dphi;
				}
			}
			return ret;
		}


		// -------------------------------------------------------------------------

		inline int    nb_nodes()      const { return _node_centers.cols(); }
		inline Vector get_node(int i) const { return _node_centers.col(i); }

		// -------------------------------------------------------------------------

		/// Each column represents p_i:  VectorX pi = _node_centers.col(i);
		MatrixDX  _node_centers;
		/// Vector of scalar values alpha
		VectorX   _alphas;
		/// Each column represents beta_i:  VectorX bi = _betas.col(i);
		MatrixDX  _betas;

	}; // END HermiteRbfReconstruction Class =======================================

}