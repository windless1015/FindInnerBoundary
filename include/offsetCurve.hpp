#pragma once
#ifndef OFFSET_CURVE_HPP
#define OFFSET_CURVE_HPP
#include <cassert>
#include <cmath>
#include <iterator>
#include <cmath>
#include <iostream>
#include <array>
#include <cmath>
#include <limits>
#include <stack>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
using namespace std;
namespace FussenAlgo
{
	namespace cavc {
		//公共算法基础
		namespace AlgoIntel
		{
			template<class TVert>
			TVert Add(const TVert& p, const TVert& q) {
				TVert cp = { p[0] + q[0],p[1] + q[1],p[2] + q[2] };
				return cp;
			}
			template<class TVert>
			TVert SubT(const TVert& p, const TVert& q) {
				TVert tp;
				tp[0] = p[0] - q[0];
				tp[1] = p[1] - q[1];
				tp[2] = p[2] - q[2];
				return tp;
			}
			template<class TVert, class ftype>
			TVert DivD(const TVert& p, const ftype& q) {
				TVert tp;
				tp[0] = p[0] / q;
				tp[1] = p[1] / q;
				tp[2] = p[2] / q;
				return tp;
			}
			template<class TVert, class ftype>
			ftype dotproduct(const TVert& p, const TVert& q) {
				ftype cp = p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
				return cp;
			}
			template<class TVert, class ftype>
			TVert ndotproduct(const ftype& q, const TVert& p) {
				TVert tp;
				tp[0] = p[0] * q;
				tp[1] = p[1] * q;
				tp[2] = p[2] * q;
				return tp;
			}
		}
	}
	//from mathutils.h
	namespace cavc {
		namespace utils {
			template <typename Real>  Real realThreshold() { return Real(1e-13); }
			template <typename Real>  Real realPrecision() { return Real(1e-10); }
			template <typename Real>  Real pi() { return Real(3.14159265358979323846264338327950288); }
			template <typename Real>  Real tau() { return Real(2) * pi<Real>(); }
			template <typename Real> bool fuzzyEqual(Real x, Real y, Real epsilon = realThreshold<Real>()) {
				return std::abs(x - y) < epsilon;
			}

			template <typename Real>
			bool fuzzyInRange(Real minValue, Real value, Real maxValue, Real epsilon = realThreshold<Real>()) {
				return (value + epsilon > minValue) && (value < maxValue + epsilon);
			}

			/// Normalize radius to be between 0 and 2PI, e.g. -PI/4 becomes 7PI/8 and 5PI becomes PI.
			template <typename Real> Real normalizeRadians(Real angle) {
				if (angle >= Real(0) && angle <= tau<Real>()) {
					return angle;
				}

				return angle - std::floor(angle / tau<Real>()) * tau<Real>();
			}

			/// Returns the smaller difference between two angles, result is negative if angle2 < angle1.
			template <typename Real> Real deltaAngle(Real angle1, Real angle2) {
				Real diff = normalizeRadians(angle2 - angle1);
				if (diff > pi<Real>()) {
					diff -= tau<Real>();
				}

				return diff;
			}

			/// Tests if angle is between a start and end angle (counter clockwise start to end, inclusive).
			template <typename Real>
			bool angleIsBetween(Real startAngle, Real endAngle, Real testAngle,
				Real epsilon = realThreshold<Real>()) {
				Real endSweep = normalizeRadians(endAngle - startAngle);
				Real midSweep = normalizeRadians(testAngle - startAngle);

				return midSweep < endSweep + epsilon;
			}

			template <typename Real>
			bool angleIsWithinSweep(Real startAngle, Real sweepAngle, Real testAngle,
				Real epsilon = realThreshold<Real>()) {
				Real endAngle = startAngle + sweepAngle;
				if (sweepAngle < Real(0)) {
					return angleIsBetween(endAngle, startAngle, testAngle, epsilon);
				}

				return angleIsBetween(startAngle, endAngle, testAngle, epsilon);
			}

			/// Returns the solutions to for the quadratic equation -b +/- sqrt (b * b - 4 * a * c) / (2 * a).
			template <typename Real>
			std::pair<Real, Real> quadraticSolutions(Real a, Real b, Real c, Real discr) {
				// Function avoids loss in precision due to taking the difference of two floating point values
				// that are very near each other in value.
				// See:
				// https://math.stackexchange.com/questions/311382/solving-a-quadratic-equation-with-precision-when-using-floating-point-variables
				//assert(fuzzyEqual(b * b - Real(4) * a * c, discr) && "discriminate is not correct");
				Real sqrtDiscr = std::sqrt(discr);
				Real denom = Real(2) * a;
				Real sol1;
				if (b < Real(0)) {
					sol1 = (-b + sqrtDiscr) / denom;
				}
				else {
					sol1 = (-b - sqrtDiscr) / denom;
				}

				Real sol2 = (c / a) / sol1;

				return std::make_pair(sol1, sol2);
			}

			template <typename T> std::size_t nextWrappingIndex(std::size_t index, const T &container) {
				if (index == container.size() - 1) {
					return 0;
				}

				return index + 1;
			}

			template <typename T> std::size_t prevWrappingIndex(std::size_t index, const T &container) {
				if (index == 0) {
					return container.size() - 1;
				}

				return index - 1;
			}
		} // namespace utils
	} // namespace cavc
	//from vector.h
	namespace cavc {
		template <typename Real, std::size_t N> class Vector {
		public:
			Vector() = default;

			Vector(std::initializer_list<Real> values) {
				if (N == values.size()) {
					std::copy(values.begin(), values.end(), m_data.begin());
				}
				else if (N > values.size()) {
					std::copy(values.begin(), values.end(), m_data.begin());
					std::fill(m_data.begin() + values.size(), m_data.end(), Real(0));
				}
				else {
					std::copy(values.begin(), values.begin() + N, m_data.begin());
				}
			}

			Vector(Real x, Real y) {
				static_assert(N == 2, "constructor for Vector2 only");
				m_data[0] = x;
				m_data[1] = y;
			}

			Vector(Real x, Real y, Real z) {
				static_assert(N == 3, "constructor for Vector3 only");
				m_data[0] = x;
				m_data[1] = y;
				m_data[2] = z;
			}

			inline Real const &operator[](std::size_t i) const { return m_data[i]; }

			inline Real &operator[](std::size_t i) { return m_data[i]; }

			inline bool operator==(Vector const &vec) const { return m_data == vec.m_data; }

			inline bool operator!=(Vector const &vec) const { return m_data != vec.m_data; }

			inline bool operator<(Vector const &vec) const { return m_data < vec.m_data; }

			inline bool operator<=(Vector const &vec) const { return m_data <= vec.m_data; }

			inline bool operator>(Vector const &vec) const { return m_data > vec.m_data; }

			inline bool operator>=(Vector const &vec) const { return m_data >= vec.m_data; }

			void makeZero() { std::fill(m_data.begin(), m_data.end(), Real(0)); }

			void makeOnes() { std::fill(m_data.begin(), m_data.end(), Real(1)); }

			void makeUnit(std::size_t d) {
				std::fill(m_data.begin(), m_data.end(), Real(0));
				m_data[d] = Real(1);
			}

			static Vector zero() {
				Vector v;
				v.makeZero();
				return v;
			}

			static Vector ones() {
				Vector v;
				v.makeOnes();
				return v;
			}

			static Vector unit(std::size_t d) {
				Vector v;
				v.makeUnit(d);
				return v;
			}

			Real &x() {
				static_assert(N >= 1, "N >= 1 to accesss x");
				return m_data[0];
			}

			Real x() const {
				static_assert(N >= 1, "N >= 1 to accesss x");
				return m_data[0];
			}

			Real &y() {
				static_assert(N >= 2, "N >= 2 to accesss y");
				return m_data[1];
			}

			Real y() const {
				static_assert(N >= 2, "N >= 2 to accesss y");
				return m_data[1];
			}

			Real &z() {
				static_assert(N >= 3, "N >= 3 to accesss z");
				return m_data[2];
			}

			Real z() const {
				static_assert(N >= 3, "N >= 3 to accesss z");
				return m_data[2];
			}

		protected:
			std::array<Real, N> m_data;
		};

		template <typename Real, std::size_t N>
		bool fuzzyZero(Vector<Real, N> const &v, Real epsilon = utils::realThreshold<Real>()) {
			bool allCompAreZero = std::abs(v[0]) < epsilon;
			for (std::size_t i = 1; i < N; ++i) {
				allCompAreZero = allCompAreZero && std::abs(v[i]) < epsilon;
			}

			return allCompAreZero;
		}

		template <typename Real, std::size_t N>
		bool fuzzyEqual(Vector<Real, N> const &v1, Vector<Real, N> const &v2,
			Real epsilon = utils::realThreshold<Real>()) {
			for (std::size_t i = 0; i < N; ++i) {
				if (!utils::fuzzyEqual(v1[i], v2[i], epsilon)) {
					return false;
				}
			}

			return true;
		}

		template <typename Real, std::size_t N> Vector<Real, N> operator+(Vector<Real, N> const &v) {
			return v;
		}

		template <typename Real, std::size_t N> Vector<Real, N> operator-(Vector<Real, N> const &v) {
			Vector<Real, N> result;
			for (std::size_t i = 0; i < N; ++i) {
				result[i] = -v[i];
			}
			return result;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator+(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
			Vector<Real, N> result = v0;
			return result += v1;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator-(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
			Vector<Real, N> result;
			for (std::size_t i = 0; i < N; i++) {
				result[i] = v0[i] - v1[i];
			}
			return result;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator*(Vector<Real, N> const &v, Real scalar) {
			Vector<Real, N> result = v;
			return result *= scalar;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator*(Real scalar, Vector<Real, N> const &v) {
			Vector<Real, N> result = v;
			return result *= scalar;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator/(Vector<Real, N> const &v, Real scalar) {
			Vector<Real, N> result = v;
			return result /= scalar;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> &operator+=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
			for (std::size_t i = 0; i < N; ++i) {
				v0[i] += v1[i];
			}
			return v0;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> &operator-=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
			for (std::size_t i = 0; i < N; ++i) {
				v0[i] -= v1[i];
			}
			return v0;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> &operator*=(Vector<Real, N> &v, Real scalar) {
			for (std::size_t i = 0; i < N; ++i) {
				v[i] *= scalar;
			}
			return v;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> &operator/=(Vector<Real, N> &v, Real scalar) {
			if (scalar != Real(0)) {
				Real invScalar = Real(1) / scalar;
				for (std::size_t i = 0; i < N; ++i) {
					v[i] *= invScalar;
				}
			}
			else {
				for (std::size_t i = 0; i < N; ++i) {
					v[i] = Real(0);
				}
			}
			return v;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator*(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
			Vector<Real, N> result = v0;
			return result *= v1;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> operator/(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
			Vector<Real, N> result = v0;
			return result /= v1;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> &operator*=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
			for (std::size_t i = 0; i < N; ++i) {
				v0[i] *= v1[i];
			}
			return v0;
		}

		template <std::size_t N, typename Real>
		Vector<Real, N> &operator/=(Vector<Real, N> &v0, Vector<Real, N> const &v1) {
			for (std::size_t i = 0; i < N; ++i) {
				v0[i] /= v1[i];
			}
			return v0;
		}

		template <std::size_t N, typename Real>
		Real dot(Vector<Real, N> const &v0, Vector<Real, N> const &v1) {
			Real dot = v0[0] * v1[0];
			for (std::size_t i = 1; i < N; ++i) {
				dot += v0[i] * v1[i];
			}
			return dot;
		}

		template <std::size_t N, typename Real> Real length(Vector<Real, N> const &v) {
			assert(!fuzzyZero(v) && "length not defined for zero vector");
			return std::sqrt(dot(v, v));
		}

		template <std::size_t N, typename Real> Real normalize(Vector<Real, N> &v) {
			assert(!fuzzyZero(v) && "normalize not defined for zero vector");
			Real length = std::sqrt(dot(v, v));
			v /= length;
			return length;
		}
	} // namespace cavc
	//from vector2.h
	namespace cavc {
		template <typename Real> using Vector2 = Vector<Real, 2>;

		/// Perpendicular vector to v (rotating counter clockwise).
		template <typename Real> Vector2<Real> perp(Vector2<Real> const &v) {
			return Vector2<Real>{-v.y(), v.x()};
		}

		/// Normalized perpendicular vector to v (rotating counter clockwise).
		template <typename Real> Vector2<Real> unitPerp(Vector2<Real> const &v) {
			Vector2<Real> result{ -v.y(), v.x() };
			normalize(result);
			return result;
		}

		/// Perpendicular dot product. Equivalent to dot(v0, perp(v1)).
		template <typename Real> Real perpDot(Vector2<Real> const &v0, Vector2<Real> const &v1) {
			return v0.x() * v1.y() - v0.y() * v1.x();
		}

		/// Returns the distance squared between p0 and p1. Equivalent to dot(p1 - p0, p1 - p0).
		template <typename Real> Real distSquared(Vector2<Real> const &p0, Vector2<Real> const &p1) {
			Vector2<Real> d = p1 - p0;
			return dot(d, d);
		}

		/// Counter clockwise angle of the vector going from p0 to p1.
		template <typename Real> Real angle(Vector2<Real> const &p0, Vector2<Real> const &p1) {
			return std::atan2(p1.y() - p0.y(), p1.x() - p0.x());
		}

		/// Returns the midpoint between p0 and p1.
		template <typename Real> Vector2<Real> midpoint(Vector2<Real> const &p0, Vector2<Real> const &p1) {
			return Vector2<Real>{(p0.x() + p1.x()) / Real(2), (p0.y() + p1.y()) / Real(2)};
		}

		/// Computes the point on the circle with radius, center, and polar angle given.
		template <typename Real>
		Vector2<Real> pointOnCircle(Real radius, Vector2<Real> const &center, Real angle) {
			return Vector2<Real>{center.x() + radius * std::cos(angle),
				center.y() + radius * std::sin(angle)};
		}

		/// Returns the closest point that lies on the line segment from p0 to p1 to the point given.
		template <typename Real>
		Vector2<Real> closestPointOnLineSeg(Vector2<Real> const &p0, Vector2<Real> const &p1,
			Vector2<Real> const &point) {
			// Dot product used to find angles
			// See: http://geomalgorithms.com/a02-_lines.html
			Vector2<Real> v = p1 - p0;
			Vector2<Real> w = point - p0;
			Real c1 = dot(w, v);
			if (c1 <= Real(0)) {
				return p0;
			}

			Real c2 = dot(v, v);
			if (c2 <= c1) {
				return p1;
			}

			Real b = c1 / c2;
			return p0 + b * v;
		}

		/// Returns true if point is left of the line pointing in the direction of the vector (p1 - p0).
		template <typename Real>
		bool isLeft(Vector2<Real> const &p0, Vector2<Real> const &p1, Vector2<Real> const &point) {
			return (p1.x() - p0.x()) * (point.y() - p0.y()) - (p1.y() - p0.y()) * (point.x() - p0.x()) > 0.0;
		}

		/// Returns true if point is left or fuzzy coincident with the line pointing in the direction of the
		/// vector (p1 - p0).
		template <typename Real>
		bool isLeftOrCoincident(Vector2<Real> const &p0, Vector2<Real> const &p1,
			Vector2<Real> const &point, Real epsilon = utils::realThreshold<Real>()) {
			return (p1.x() - p0.x()) * (point.y() - p0.y()) - (p1.y() - p0.y()) * (point.x() - p0.x()) >
				-epsilon;
		}

		/// Returns true if point is right or fuzzy coincident with the line pointing in the direction of
		/// the vector (p1 - p0).
		template <typename Real>
		bool isRightOrCoincident(Vector2<Real> const &p0, Vector2<Real> const &p1,
			Vector2<Real> const &point, Real epsilon = utils::realThreshold<Real>()) {
			return (p1.x() - p0.x()) * (point.y() - p0.y()) - (p1.y() - p0.y()) * (point.x() - p0.x()) <
				epsilon;
		}

		/// Test if a point is within a arc sweep angle region defined by center, start, end, and bulge.
		template <typename Real>
		bool pointWithinArcSweepAngle(Vector2<Real> const &center, Vector2<Real> const &arcStart,
			Vector2<Real> const &arcEnd, Real bulge, Vector2<Real> const &point) {
			assert(std::abs(bulge) > utils::realThreshold<Real>() && "expected arc");
			assert(std::abs(bulge) <= Real(1) && "bulge should always be between -1 and 1");

			if (bulge > Real(0)) {
				return isLeftOrCoincident(center, arcStart, point) &&
					isRightOrCoincident(center, arcEnd, point);
			}

			return isRightOrCoincident(center, arcStart, point) && isLeftOrCoincident(center, arcEnd, point);
		}
	} // namespace cavc
	//from staticspatialindex.h
	namespace cavc {
		template <typename Real, std::size_t NodeSize = 16> class StaticSpatialIndex {
		public:
			StaticSpatialIndex(std::size_t numItems) {
				assert(numItems > 0 && "number of items must be greater than 0");
				static_assert(NodeSize >= 2 && NodeSize <= 65535, "node size must be between 2 and 65535");
				// calculate the total number of nodes in the R-tree to allocate space for
				// and the index of each tree level (used in search later)
				m_numItems = numItems;
				std::size_t n = numItems;
				std::size_t numNodes = numItems;
				m_levelBounds.push_back(n * 4);
				do {
					n = static_cast<std::size_t>(std::ceil(static_cast<float>(n) / NodeSize));
					numNodes += n;
					m_levelBounds.push_back(numNodes * 4);
				} while (n != 1);

				m_boxes.resize(numNodes * 4);
				m_indices.resize(numNodes);
				m_pos = 0;
				m_minX = std::numeric_limits<Real>::infinity();
				m_minY = std::numeric_limits<Real>::infinity();
				m_maxX = -1 * std::numeric_limits<Real>::infinity();
				m_maxY = -1 * std::numeric_limits<Real>::infinity();
			}

			Real MinX() const { return m_minX; }
			Real MinY() const { return m_minY; }
			Real MaxX() const { return m_maxX; }
			Real MaxY() const { return m_maxY; }

			void add(Real minX, Real minY, Real maxX, Real maxY) {
				std::size_t index = m_pos >> 2;
				m_indices[index] = index;
				m_boxes[m_pos++] = minX;
				m_boxes[m_pos++] = minY;
				m_boxes[m_pos++] = maxX;
				m_boxes[m_pos++] = maxY;

				if (minX < m_minX)
					m_minX = minX;
				if (minY < m_minY)
					m_minY = minY;
				if (maxX > m_maxX)
					m_maxX = maxX;
				if (maxY > m_maxY)
					m_maxY = maxY;
			}

			void finish() {
				assert(m_pos >> 2 == m_numItems && "added item count should equal static size given");
				Real width = m_maxX - m_minX;
				Real height = m_maxY - m_minY;
				std::vector<std::uint32_t> hilbertValues(m_numItems);

				std::size_t pos = 0;

				for (std::size_t i = 0; i < m_numItems; ++i) {
					pos = 4 * i;
					Real minX = m_boxes[pos++];
					Real minY = m_boxes[pos++];
					Real maxX = m_boxes[pos++];
					Real maxY = m_boxes[pos++];

					// hilbert max input value for x and y
					const Real hilbertMax = static_cast<Real>((1 << 16) - 1);
					// mapping the x and y coordinates of the center of the box to values in the range
					// [0 -> n - 1] such that the min of the entire set of bounding boxes maps to 0 and the max of
					// the entire set of bounding boxes maps to n - 1 our 2d space is x: [0 -> n-1] and
					// y: [0 -> n-1], our 1d hilbert curve value space is d: [0 -> n^2 - 1]
					Real x = std::floor(hilbertMax * ((minX + maxX) / 2 - m_minX) / width);
					std::uint32_t hx = static_cast<std::uint32_t>(x);
					Real y = std::floor(hilbertMax * ((minY + maxY) / 2 - m_minY) / height);
					std::uint32_t hy = static_cast<std::uint32_t>(y);
					hilbertValues[i] = hilbertXYToIndex(hx, hy);
				}

				// sort items by their Hilbert value (for packing later)
				sort(hilbertValues, m_boxes, m_indices, 0, m_numItems - 1);

				// generate nodes at each tree level, bottom-up
				pos = 0;
				for (std::size_t i = 0; i < m_levelBounds.size() - 1; i++) {
					auto end = m_levelBounds[i];

					// generate a parent node for each block of consecutive <nodeSize> nodes
					while (pos < end) {
						auto nodeMinX = std::numeric_limits<Real>::infinity();
						auto nodeMinY = std::numeric_limits<Real>::infinity();
						auto nodeMaxX = -1 * std::numeric_limits<Real>::infinity();
						auto nodeMaxY = -1 * std::numeric_limits<Real>::infinity();
						auto nodeIndex = pos;

						// calculate bbox for the new node
						for (std::size_t j = 0; j < NodeSize && pos < end; j++) {
							auto minX = m_boxes[pos++];
							auto minY = m_boxes[pos++];
							auto maxX = m_boxes[pos++];
							auto maxY = m_boxes[pos++];
							if (minX < nodeMinX)
								nodeMinX = minX;
							if (minY < nodeMinY)
								nodeMinY = minY;
							if (maxX > nodeMaxX)
								nodeMaxX = maxX;
							if (maxY > nodeMaxY)
								nodeMaxY = maxY;
						}

						// add the new node to the tree data
						m_indices[m_pos >> 2] = nodeIndex;
						m_boxes[m_pos++] = nodeMinX;
						m_boxes[m_pos++] = nodeMinY;
						m_boxes[m_pos++] = nodeMaxX;
						m_boxes[m_pos++] = nodeMaxY;
					}
				}
			}

			// Visit all the bounding boxes in the spatial index. Visitor function has the signature
			// void(Real xmin, Real ymin, Real xmax, Real ymax, std::size_t level).
			template <typename F> void visitBoundingBoxes(F &&visitor) const {
				std::size_t nodeIndex = m_boxes.size() - 4;
				std::size_t level = m_levelBounds.size() - 1;

				std::vector<std::size_t> stack;
				stack.reserve(16);

				bool done = false;
				while (!done) {
					auto end = std::min(nodeIndex + NodeSize * 4, m_levelBounds[level]);
					for (std::size_t pos = nodeIndex; pos < end; pos += 4) {
						auto index = m_indices[pos >> 2];
						visitor(m_boxes[pos], m_boxes[pos + 1], m_boxes[pos + 2], m_boxes[pos + 3], level);

						if (nodeIndex >= m_numItems * 4) {
							stack.push_back(index);
							stack.push_back(level - 1);
						}
					}

					if (stack.size() > 1) {
						level = stack.back();
						stack.pop_back();
						nodeIndex = stack.back();
						stack.pop_back();
					}
					else {
						done = true;
					}
				}
			}

			// Query the spatial index adding indexes to the results vector given.
			void query(Real minX, Real minY, Real maxX, Real maxY, std::vector<std::size_t> &results) const {
				auto visitor = [&](std::size_t index) {
					results.push_back(index);
					return true;
				};

				visitQuery(minX, minY, maxX, maxY, visitor);
			}

			// Query the spatial index, invoking a visitor function for each index that overlaps the bounding
			// box given. Visitor function has the signature bool(std::size_t index), if visitor returns false
			// the query stops early, otherwise the query continues.
			template <typename F>
			void visitQuery(Real minX, Real minY, Real maxX, Real maxY, F &&visitor) const {
				assert(m_pos == m_boxes.size() && "data not yet indexed - call Finish() before querying");

				auto nodeIndex = m_boxes.size() - 4;
				auto level = m_levelBounds.size() - 1;

				// stack for traversing nodes
				std::vector<std::size_t> stack;
				// reserve some space to avoid repeated small allocations
				stack.reserve(16);

				auto done = false;

				while (!done) {
					// find the end index of the node
					auto end = std::min(nodeIndex + NodeSize * 4, m_levelBounds[level]);

					// search through child nodes
					for (std::size_t pos = nodeIndex; pos < end; pos += 4) {
						auto index = m_indices[pos >> 2];
						// check if node bbox intersects with query bbox
						if (maxX < m_boxes[pos])
							continue; // maxX < nodeMinX
						if (maxY < m_boxes[pos + 1])
							continue; // maxY < nodeMinY
						if (minX > m_boxes[pos + 2])
							continue; // minX > nodeMaxX
						if (minY > m_boxes[pos + 3])
							continue; // minY > nodeMaxY

						if (nodeIndex < m_numItems * 4) {
							done = !visitor(index);
							if (done) {
								break;
							}
						}
						else {
							// push node index and level for further traversal
							stack.push_back(index);
							stack.push_back(level - 1);
						}
					}

					if (stack.size() > 1) {
						level = stack.back();
						stack.pop_back();
						nodeIndex = stack.back();
						stack.pop_back();
					}
					else {
						done = true;
					}
				}
			}

			static std::uint32_t hilbertXYToIndex(std::uint32_t x, std::uint32_t y) {
				std::uint32_t a = x ^ y;
				std::uint32_t b = 0xFFFF ^ a;
				std::uint32_t c = 0xFFFF ^ (x | y);
				std::uint32_t d = x & (y ^ 0xFFFF);

				std::uint32_t A = a | (b >> 1);
				std::uint32_t B = (a >> 1) ^ a;
				std::uint32_t C = ((c >> 1) ^ (b & (d >> 1))) ^ c;
				std::uint32_t D = ((a & (c >> 1)) ^ (d >> 1)) ^ d;

				a = A;
				b = B;
				c = C;
				d = D;
				A = ((a & (a >> 2)) ^ (b & (b >> 2)));
				B = ((a & (b >> 2)) ^ (b & ((a ^ b) >> 2)));
				C ^= ((a & (c >> 2)) ^ (b & (d >> 2)));
				D ^= ((b & (c >> 2)) ^ ((a ^ b) & (d >> 2)));

				a = A;
				b = B;
				c = C;
				d = D;
				A = ((a & (a >> 4)) ^ (b & (b >> 4)));
				B = ((a & (b >> 4)) ^ (b & ((a ^ b) >> 4)));
				C ^= ((a & (c >> 4)) ^ (b & (d >> 4)));
				D ^= ((b & (c >> 4)) ^ ((a ^ b) & (d >> 4)));

				a = A;
				b = B;
				c = C;
				d = D;
				C ^= ((a & (c >> 8)) ^ (b & (d >> 8)));
				D ^= ((b & (c >> 8)) ^ ((a ^ b) & (d >> 8)));

				a = C ^ (C >> 1);
				b = D ^ (D >> 1);

				std::uint32_t i0 = x ^ y;
				std::uint32_t i1 = b | (0xFFFF ^ (i0 | a));

				i0 = (i0 | (i0 << 8)) & 0x00FF00FF;
				i0 = (i0 | (i0 << 4)) & 0x0F0F0F0F;
				i0 = (i0 | (i0 << 2)) & 0x33333333;
				i0 = (i0 | (i0 << 1)) & 0x55555555;

				i1 = (i1 | (i1 << 8)) & 0x00FF00FF;
				i1 = (i1 | (i1 << 4)) & 0x0F0F0F0F;
				i1 = (i1 | (i1 << 2)) & 0x33333333;
				i1 = (i1 | (i1 << 1)) & 0x55555555;

				return (i1 << 1) | i0;
			}

		private:
			Real m_minX;
			Real m_minY;
			Real m_maxX;
			Real m_maxY;
			std::size_t m_numItems;
			std::vector<std::size_t> m_levelBounds;
			std::vector<Real> m_boxes;
			std::vector<std::size_t> m_indices;
			std::size_t m_pos;

			static void sort(std::vector<std::uint32_t> &values, std::vector<Real> &boxes,
				std::vector<std::size_t> &indices, std::size_t left, std::size_t right) {

				if (left >= right)
					return;

				auto pivot = values[(left + right) >> 1];
				auto i = left - 1;
				auto j = right + 1;

				while (true) {
					do
						i++;
					while (values[i] < pivot);
					do
						j--;
					while (values[j] > pivot);
					if (i >= j)
						break;
					swap(values, boxes, indices, i, j);
				}

				sort(values, boxes, indices, left, j);
				sort(values, boxes, indices, j + 1, right);
			}

			static void swap(std::vector<std::uint32_t> &values, std::vector<Real> &boxes,
				std::vector<std::size_t> &indices, std::size_t i, std::size_t j) {
				auto temp = values[i];
				values[i] = values[j];
				values[j] = temp;

				auto k = 4 * i;
				auto m = 4 * j;

				auto a = boxes[k];
				auto b = boxes[k + 1];
				auto c = boxes[k + 2];
				auto d = boxes[k + 3];
				boxes[k] = boxes[m];
				boxes[k + 1] = boxes[m + 1];
				boxes[k + 2] = boxes[m + 2];
				boxes[k + 3] = boxes[m + 3];
				boxes[m] = a;
				boxes[m + 1] = b;
				boxes[m + 2] = c;
				boxes[m + 3] = d;

				auto e = indices[i];
				indices[i] = indices[j];
				indices[j] = e;
			}
		};
	}
	//from intrcircle2circle2.h
	namespace cavc {
		enum class Circle2Circle2IntrType {
			// no intersect between circles
			NoIntersect,
			// one intersect between circles (tangent)
			OneIntersect,
			// two intersects between circles
			TwoIntersects,
			// circles are coincident
			Coincident
		};

		template <typename Real> struct IntrCircle2Circle2Result {
			// type of intersect
			Circle2Circle2IntrType intrType;
			// first intersect point if intrType is OneIntersect or TwoIntersects, undefined otherwise
			Vector2<Real> point1;
			// second intersect point if intrType is TwoIntersects, undefined otherwise
			Vector2<Real> point2;
		};

		// Find intersect between two circles in 2D.
		template <typename Real>
		IntrCircle2Circle2Result<Real> intrCircle2Circle2(Real radius1, Vector2<Real> const &center1,
			Real radius2, Vector2<Real> const &center2) {
			// Reference algorithm: http://paulbourke.net/geometry/circlesphere/

			IntrCircle2Circle2Result<Real> result;
			Vector2<Real> cv = center2 - center1;
			Real d2 = dot(cv, cv);
			Real d = std::sqrt(d2);
			if (d < utils::realThreshold<Real>()) {
				// same center position
				if (utils::fuzzyEqual(radius1, radius2)) {
					result.intrType = Circle2Circle2IntrType::Coincident;
				}
				else {
					result.intrType = Circle2Circle2IntrType::NoIntersect;
				}
			}
			else {
				// different center position
				if (d > radius1 + radius2 + utils::realThreshold<Real>() ||
					d + utils::realThreshold<Real>() < std::abs(radius1 - radius2)) {
					result.intrType = Circle2Circle2IntrType::NoIntersect;
				}
				else {
					Real rad1Sq = radius1 * radius1;
					Real a = (rad1Sq - radius2 * radius2 + d2) / (Real(2) * d);
					Vector2<Real> midPoint = center1 + a * cv / d;
					Real diff = rad1Sq - a * a;
					if (diff < Real(0)) {
						result.intrType = Circle2Circle2IntrType::OneIntersect;
						result.point1 = midPoint;
					}
					else {
						Real h = std::sqrt(diff);
						Real hOverD = h / d;
						Real xTerm = hOverD * cv.y();
						Real yTerm = hOverD * cv.x();
						Real x1 = midPoint.x() + xTerm;
						Real y1 = midPoint.y() - yTerm;
						Real x2 = midPoint.x() - xTerm;
						Real y2 = midPoint.y() + yTerm;
						result.point1 = Vector2<Real>(x1, y1);
						result.point2 = Vector2<Real>(x2, y2);
						if (fuzzyEqual(result.point1, result.point2)) {
							result.intrType = Circle2Circle2IntrType::OneIntersect;
						}
						else {
							result.intrType = Circle2Circle2IntrType::TwoIntersects;
						}
					}
				}
			}

			return result;
		}
	}
	//from intrlineseg2circle2.h
	namespace cavc {
		template <typename Real> struct IntrLineSeg2Circle2Result {
			// number of interescts found (0, 1, or 2)
			int numIntersects;
			// parametric value for first intersect (if numIntersects > 0) otherwise undefined
			Real t0;
			// parametric value for second intersect (if numintersects > 1) otherwise undefined
			Real t1;
		};

		// Gets the intersect between a segment and a circle, returning the parametric solution t to the
		// segment equation P(t) = v1 + t * (v2 - v1) for t = 0 to t = 1, if t < 0 or t > 1 then intersect
		// occurs only when extending the segment out past the points given (if t < 0 intersect nearest v1,
		// if t > 0 then intersect nearest v2), intersects are "sticky" and "snap" to tangent points, e.g. a
		// segment very close to being a tangent will be returned as a single intersect point
		template <typename Real>
		IntrLineSeg2Circle2Result<Real> intrLineSeg2Circle2(Vector2<Real> const &p0,
			Vector2<Real> const &p1, Real radius,
			Vector2<Real> const &circleCenter) {
			// This function solves for t by substituting the parametric equations for the segment x = v1.X +
			// t * (v2.X - v1.X) and y = v1.Y + t * (v2.Y - v1.Y) for t = 0 to t = 1 into the circle equation
			// (x-h)^2 + (y-k)^2 = r^2 and then solving the resulting equation in the form a*t^2 + b*t + c = 0
			// using the quadratic formula
			IntrLineSeg2Circle2Result<Real> result;
			Real dx = p1.x() - p0.x();
			Real dy = p1.y() - p0.y();
			Real h = circleCenter.x();
			Real k = circleCenter.y();

			Real a = dx * dx + dy * dy;
			if (std::abs(a) < utils::realThreshold<Real>()) {
				// v1 = v2, test if point is on the circle
				Real xh = p0.x() - h;
				Real yk = p0.y() - k;
				if (utils::fuzzyEqual(xh * xh + yk * yk, radius * radius)) {
					result.numIntersects = 1;
					result.t0 = Real(0);
				}
				else {
					result.numIntersects = 0;
				}
			}
			else {
				Real b = Real(2) * (dx * (p0.x() - h) + dy * (p0.y() - k));
				Real c = (p0.x() * p0.x() - 2.0 * h * p0.x() + h * h) +
					(p0.y() * p0.y() - 2.0 * k * p0.y() + k * k) - radius * radius;
				Real discr = b * b - 4.0 * a * c;

				if (std::abs(discr) < utils::realThreshold<Real>()) {
					// 1 solution (tangent line)
					result.numIntersects = 1;
					result.t0 = -b / (Real(2) * a);
				}
				else if (discr < Real(0)) {
					result.numIntersects = 0;
				}
				else {
					result.numIntersects = 2;
					std::pair<Real, Real> sols = utils::quadraticSolutions(a, b, c, discr);
					result.t0 = sols.first;
					result.t1 = sols.second;
				}
			}

			assert(result.numIntersects >= 0 && result.numIntersects <= 2);
			return result;
		}
	}
	//from intrlineseg3linese2.h
	namespace cavc {
		enum class LineSeg2LineSeg2IntrType {
			// no intersect (segments are parallel and not collinear)
			None,
			// true intersect between line segments
			True,
			// segments overlap each other by some amount
			Coincident,
			// false intersect between line segments (one or both of the segments must be extended)
			False
		};

		template <typename Real> struct IntrLineSeg2LineSeg2Result {
			// holds the type of intersect, if True or False then point holds the point that they intersect,
			// if True then t0 and t1 are undefined, if False then t0 is the parametric value of the first
			// segment and t1 is the parametric value of the second segment, if Coincident then point is
			// undefined and t0 holds the parametric value start of coincidence and t1 holds the parametric
			// value of the end of the coincidence for the first segments equation
			LineSeg2LineSeg2IntrType intrType;
			Real t0;
			Real t1;
			Vector2<Real> point;
		};

		template <typename Real>
		IntrLineSeg2LineSeg2Result<Real>
			intrLineSeg2LineSeg2(Vector2<Real> const &u1, Vector2<Real> const &u2, Vector2<Real> const &v1,
				Vector2<Real> const &v2) {
			// This implementation works by processing the segments in parametric equation form and using
			// perpendicular products
			// see: http://geomalgorithms.com/a05-_intersect-1.html and
			// http://mathworld.wolfram.com/PerpDotProduct.html

			IntrLineSeg2LineSeg2Result<Real> result;
			Vector2<Real> u = u2 - u1;
			Vector2<Real> v = v2 - v1;
			Real d = perpDot(u, v);

			Vector2<Real> w = u1 - v1;

			// Test if point is inside a segment, NOTE: assumes points are aligned
			auto isInSegment = [](Vector2<Real> const &pt, Vector2<Real> const &segStart,
				Vector2<Real> const &segEnd) {
				if (utils::fuzzyEqual(segStart.x(), segEnd.x())) {
					// vertical segment, test y coordinate
					auto minMax = std::minmax(segStart.y(), segEnd.y());
					return utils::fuzzyInRange(minMax.first, pt.y(), minMax.second);
				}

				// else just test x coordinate
				auto minMax = std::minmax(segStart.x(), segEnd.x());
				return utils::fuzzyInRange(minMax.first, pt.x(), minMax.second);
			};

			// threshold check here to avoid almost parallel lines resulting in very distant intersection
			if (std::abs(d) > utils::realThreshold<Real>()) {
				// segments not parallel or collinear
				result.t0 = perpDot(v, w) / d;
				result.t1 = perpDot(u, w) / d;
				result.point = v1 + result.t1 * v;
				if (result.t0 + utils::realThreshold<Real>() < Real(0) ||
					result.t0 > Real(1) + utils::realThreshold<Real>() ||
					result.t1 + utils::realThreshold<Real>() < Real(0) ||
					result.t1 > Real(1) + utils::realThreshold<Real>()) {
					result.intrType = LineSeg2LineSeg2IntrType::False;
				}
				else {
					result.intrType = LineSeg2LineSeg2IntrType::True;
				}
			}
			else {
				// segments are parallel or collinear
				Real a = perpDot(u, w);
				Real b = perpDot(v, w);
				// threshold check here, we consider almost parallel lines to be parallel
				if (std::abs(a) > utils::realThreshold<Real>() || std::abs(b) > utils::realThreshold<Real>()) {
					// parallel and not collinear so no intersect
					result.intrType = LineSeg2LineSeg2IntrType::None;
				}
				else {
					// either collinear or degenerate (segments are single points)
					bool uIsPoint = fuzzyEqual(u1, u2);
					bool vIsPoint = fuzzyEqual(v1, v2);
					if (uIsPoint && vIsPoint) {
						// both segments are just points
						if (fuzzyEqual(u1, v1)) {
							// same point
							result.point = u1;
							result.intrType = LineSeg2LineSeg2IntrType::True;
						}
						else {
							// distinct points
							result.intrType = LineSeg2LineSeg2IntrType::None;
						}

					}
					else if (uIsPoint) {
						if (isInSegment(u1, v1, v2)) {
							result.intrType = LineSeg2LineSeg2IntrType::True;
							result.point = u1;
						}
						else {
							result.intrType = LineSeg2LineSeg2IntrType::None;
						}

					}
					else if (vIsPoint) {
						if (isInSegment(v1, u1, u2)) {
							result.intrType = LineSeg2LineSeg2IntrType::True;
							result.point = v1;
						}
						else {
							result.intrType = LineSeg2LineSeg2IntrType::None;
						}
					}
					else {
						// neither segment is a point, check if they overlap
						Vector2<Real> w2 = u2 - v1;
						if (std::abs(v.x()) < utils::realThreshold<Real>()) {
							result.t0 = w.y() / v.y();
							result.t1 = w2.y() / v.y();
						}
						else {
							result.t0 = w.x() / v.x();
							result.t1 = w2.x() / v.x();
						}

						if (result.t0 > result.t1) {
							using std::swap;
							swap(result.t0, result.t1);
						}

						// using threshold check here to make intersect "sticky" to prefer considering it an
						// intersect
						if (result.t0 > Real(1) + utils::realThreshold<Real>() ||
							result.t1 + utils::realThreshold<Real>() < Real(0)) {
							// no overlap
							result.intrType = LineSeg2LineSeg2IntrType::None;
						}
						else {
							result.t0 = std::max(result.t0, Real(0));
							result.t1 = std::min(result.t1, Real(1));
							if (std::abs(result.t1 - result.t0) < utils::realThreshold<Real>()) {
								// intersect is a single point (segments line up end to end)
								result.intrType = LineSeg2LineSeg2IntrType::True;
								result.point = v1 + result.t0 * v;
							}
							else {
								result.intrType = LineSeg2LineSeg2IntrType::Coincident;
							}
						}
					}
				}
			}

			return result;
		}
	}
	//from polyline.h
	namespace cavc {

		template <typename Real> class PlineVertex {
		public:
			PlineVertex() = default;
			PlineVertex(Real x, Real y, Real bulge) : m_position(x, y), m_bulge(bulge) {}
			PlineVertex(Vector2<Real> position, Real bulge)
				: PlineVertex(position.x(), position.y(), bulge) {}

			Real x() const { return m_position.x(); }
			Real &x() { return m_position.x(); }

			Real y() const { return m_position.y(); }
			Real &y() { return m_position.y(); }

			Real bulge() const { return m_bulge; }
			Real &bulge() { return m_bulge; }

			bool bulgeIsZero(Real epsilon = utils::realPrecision<Real>()) const {
				return std::abs(m_bulge) < epsilon;
			}

			Vector2<Real> const &pos() const { return m_position; }
			Vector2<Real> &pos() { return m_position; }

		private:
			Vector2<Real> m_position;
			Real m_bulge;
		};

		template <typename Real> class Polyline {
		public:
			Polyline() : m_isClosed(false), m_vertexes() {}

			using PVertex = PlineVertex<Real>;

			inline PVertex const &operator[](std::size_t i) const { return m_vertexes[i]; }

			inline PVertex &operator[](std::size_t i) { return m_vertexes[i]; }

			bool isClosed() const { return m_isClosed; }
			bool &isClosed() { return m_isClosed; }

			void addVertex(Real x, Real y, Real bulge) { m_vertexes.emplace_back(x, y, bulge); }
			void addVertex(PVertex vertex) { addVertex(vertex.x(), vertex.y(), vertex.bulge()); }

			std::size_t size() const { return m_vertexes.size(); }

			PVertex const &lastVertex() const { return m_vertexes.back(); }
			PVertex &lastVertex() { return m_vertexes.back(); }

			std::vector<PVertex> &vertexes() { return m_vertexes; }
			std::vector<PVertex> const &vertexes() const { return m_vertexes; }

		private:
			bool m_isClosed;
			std::vector<PVertex> m_vertexes;
		};

		/// Iterate the segement indices of a polyline. f is invoked for each segment index pair, iteration
		/// stops when all indices have been iterated or f returns false. f signature is bool(std::size_t,
		/// std::size_t).
		template <typename Real, typename F> void iterateSegIndices(Polyline<Real> const &pline, F &&f) {
			std::size_t i;
			std::size_t j;
			if (pline.isClosed()) {
				i = 0;
				j = pline.vertexes().size() - 1;
			}
			else {
				i = 1;
				j = 0;
			}

			while (i < pline.vertexes().size() && f(j, i)) {
				j = i;
				i = i + 1;
			}
		}

		/// Result from computing the arc radius and arc center of a segment.
		template <typename Real> struct ArcRadiusAndCenter {
			Real radius;
			Vector2<Real> center;
		};

		/// Compute the arc radius and arc center of a arc segment defined by v1 to v2.
		template <typename Real>
		ArcRadiusAndCenter<Real> arcRadiusAndCenter(PlineVertex<Real> const &v1,
			PlineVertex<Real> const &v2) {
			assert(!v1.bulgeIsZero() && "v1 to v2 must be an arc");
			assert(!fuzzyEqual(v1.pos(), v2.pos()) && "v1 must not equal v2");

			// compute radius
			Real b = std::abs(v1.bulge());
			Vector2<Real> v = v2.pos() - v1.pos();
			Real d = length(v);
			Real r = d * (b * b + Real(1)) / (Real(4) * b);

			// compute center
			Real s = b * d / Real(2);
			Real m = r - s;
			Real offsX = -m * v.y() / d;
			Real offsY = m * v.x() / d;
			if (v1.bulge() < Real(0)) {
				offsX = -offsX;
				offsY = -offsY;
			}

			Vector2<Real> c{ v1.x() + v.x() / Real(2) + offsX, v1.y() + v.y() / Real(2) + offsY };
			return ArcRadiusAndCenter<Real>{r, c};
		}

		/// Result of splitting a segment v1 to v2.
		template <typename Real> struct SplitResult {
			/// Updated starting vertex.
			PlineVertex<Real> updatedStart;
			/// Vertex at the split point.
			PlineVertex<Real> splitVertex;
		};

		/// Split the segment defined by v1 to v2 at some point defined along it.
		template <typename Real>
		SplitResult<Real> splitAtPoint(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2,
			Vector2<Real> const &point) {
			SplitResult<Real> result;
			if (v1.bulgeIsZero()) {
				result.updatedStart = v1;
				result.splitVertex = PlineVertex<Real>(point, Real(0));
			}
			else if (fuzzyEqual(v1.pos(), v2.pos()) || fuzzyEqual(v1.pos(), point)) {
				result.updatedStart = PlineVertex<Real>(point, Real(0));
				result.splitVertex = PlineVertex<Real>(point, v1.bulge());
			}
			else if (fuzzyEqual(v2.pos(), point)) {
				result.updatedStart = v1;
				result.splitVertex = PlineVertex<Real>(v2.pos(), Real(0));
			}
			else {
				auto radiusAndCenter = arcRadiusAndCenter(v1, v2);
				Vector2<Real> arcCenter = radiusAndCenter.center;
				Real a = angle(arcCenter, point);
				Real arcStartAngle = angle(arcCenter, v1.pos());
				Real theta1 = utils::deltaAngle(arcStartAngle, a);
				Real bulge1 = std::tan(theta1 / Real(4));
				Real arcEndAngle = angle(arcCenter, v2.pos());
				Real theta2 = utils::deltaAngle(a, arcEndAngle);
				Real bulge2 = std::tan(theta2 / Real(4));

				result.updatedStart = PlineVertex<Real>(v1.pos(), bulge1);
				result.splitVertex = PlineVertex<Real>(point, bulge2);
			}

			return result;
		}

		/// Axis aligned bounding box (AABB).
		template <typename Real> struct AABB {
			Real xMin;
			Real yMin;
			Real xMax;
			Real yMax;

			void expand(Real val) {
				xMin -= val;
				yMin -= val;
				xMax += val;
				yMax += val;
			}
		};

		/// Compute the extents of a polyline.
		template <typename Real> AABB<Real> extents(Polyline<Real> const &pline) {
			if (pline.vertexes().size() < 2) {
				return AABB<Real>{pline[0].x(), pline[0].y(), pline[0].x(), pline[0].y()};
			}

			AABB<Real> result{ std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(),
							  std::numeric_limits<Real>::min(), std::numeric_limits<Real>::min() };

			auto visitor = [&](std::size_t i, std::size_t j) {
				PlineVertex<Real> const &v1 = pline[i];
				if (v1.bulgeIsZero()) {
					if (v1.x() < result.xMin)
						result.xMin = v1.x();
					if (v1.y() < result.yMin)
						result.yMin = v1.y();
					if (v1.x() > result.yMax)
						result.xMax = v1.x();
					if (v1.y() > result.yMax)
						result.yMax = v1.y();
				}
				else {
					PlineVertex<Real> const &v2 = pline[j];
					auto arc = arcRadiusAndCenter(v1, v2);

					Real startAngle = angle(arc.center, v1.pos());
					Real endAngle = angle(arc.center, v2.pos());
					Real sweepAngle = utils::deltaAngle(startAngle, endAngle);

					Real arcXMin, arcYMin, arcXMax, arcYMax;

					// crosses PI/2
					if (utils::angleIsWithinSweep(startAngle, sweepAngle, Real(0.5) * utils::pi<Real>())) {
						arcYMax = arc.center.y() + arc.radius;
					}
					else {
						arcYMax = std::max(v1.y(), v2.y());
					}

					// crosses PI
					if (utils::angleIsWithinSweep(startAngle, sweepAngle, utils::pi<Real>())) {
						arcXMin = arc.center.x() - arc.radius;
					}
					else {
						arcXMin = std::min(v1.x(), v2.x());
					}

					// crosses 3PI/2
					if (utils::angleIsWithinSweep(startAngle, sweepAngle, Real(1.5) * utils::pi<Real>())) {
						arcYMin = arc.center.y() - arc.radius;
					}
					else {
						arcYMin = std::min(v1.y(), v2.y());
					}

					// crosses 2PI
					if (utils::angleIsWithinSweep(startAngle, sweepAngle, Real(2) * utils::pi<Real>())) {
						arcXMax = arc.center.x() + arc.radius;
					}
					else {
						arcXMax = std::max(v1.x(), v2.x());
					}

					if (arcXMin < result.xMin)
						result.xMin = arcXMin;
					if (arcYMin < result.yMin)
						result.yMin = arcYMin;
					if (arcXMax > result.yMax)
						result.xMax = arcXMax;
					if (arcYMax > result.yMax)
						result.yMax = arcYMax;
				}

				// return true to iterate all segments
				return true;
			};

			iterateSegIndices(pline, visitor);

			return result;
		}

		/// Compute the closest point on a segment defined by v1 to v2 to the point given.
		template <typename Real>
		Vector2<Real> closestPointOnSeg(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2,
			Vector2<Real> const &point) {
			if (v1.bulgeIsZero()) {
				return closestPointOnLineSeg(v1.pos(), v2.pos(), point);
			}

			auto arc = arcRadiusAndCenter(v1, v2);

			if (fuzzyEqual(point, arc.center)) {
				// avoid normalizing zero length vector (point is at center, just return start point)
				return v1.pos();
			}

			if (pointWithinArcSweepAngle(arc.center, v1.pos(), v2.pos(), v1.bulge(), point)) {
				// closest point is on the arc
				Vector2<Real> vToPoint = point - arc.center;
				normalize(vToPoint);
				return arc.radius * vToPoint + arc.center;
			}

			// else closest point is one of the ends
			Real dist1 = distSquared(v1.pos(), point);
			Real dist2 = distSquared(v2.pos(), point);
			if (dist1 < dist2) {
				return v1.pos();
			}

			return v2.pos();
		}

		/// Result from computing the closest point and index.
		template <typename Real> struct ClosestPointAndIndex {
			std::size_t index;
			Vector2<Real> point;
			Real distance;
		};

		/// Computes the closest point and starting vertex index of the segment that point lies on to the
		/// point given.
		template <typename Real>
		ClosestPointAndIndex<Real> closestPointAndIndex(Polyline<Real> const &pline,
			Vector2<Real> const &point) {
			assert(pline.vertexes().size > 0 && "empty polyline has no closest point");
			ClosestPointAndIndex<Real> result;
			if (pline.vertexes().size() == 1) {
				result.index = 0;
				result.point = pline[0];
				result.distance = length(point - pline[0].pos());
				return result;
			}

			result.distance = std::numeric_limits<Real>::infinity();

			auto visitor = [&](std::size_t i, std::size_t j) {
				Vector2<Real> cp = closestPointOnSeg(pline[i], pline[j], point);
				auto diffVec = point - cp;
				Real dist2 = dot(diffVec, diffVec);
				if (dist2 < result.distance) {
					result.index = i;
					result.point = cp;
					result.distance = dist2;
				}

				// iterate all segments
				return true;
			};

			iterateSegIndices(pline, visitor);
			// we used the squared distance while iterating and comparing, take sqrt for actual distance
			result.distance = std::sqrt(result.distance);
			return result;
		}

		/// Returns a new polyline with all arc segments converted to line segments, error is the maximum
		/// distance from any line segment to the arc it is approximating. Line segments are circumscribed
		/// by the arc (all end points lie on the arc path).
		template <typename Real>
		Polyline<Real> convertArcsToLines(Polyline<Real> const &pline, Real error) {
			cavc::Polyline<Real> result;
			result.isClosed() = pline.isClosed();
			auto visitor = [&](std::size_t i, std::size_t j) {
				const auto &v1 = pline[i];
				const auto &v2 = pline[j];
				if (v1.bulgeIsZero()) {
					result.addVertex(v1);
				}
				else {

					auto arc = arcRadiusAndCenter(v1, v2);
					auto startAngle = angle(arc.center, v1.pos());
					auto endAngle = angle(arc.center, v2.pos());
					Real deltaAngle = std::abs(cavc::utils::deltaAngle(startAngle, endAngle));

					error = std::abs(error);
					Real segmentSubAngle = std::abs(Real(2) * std::acos(Real(1) - error / arc.radius));
					std::size_t segmentCount = static_cast<std::size_t>(std::ceil(deltaAngle / segmentSubAngle));
					// update segment subangle for equal length segments
					segmentSubAngle = deltaAngle / segmentCount;

					if (v1.bulge() < Real(0)) {
						segmentSubAngle = -segmentSubAngle;
					}
					// add the start point
					result.addVertex(v1.x(), v1.y(), 0.0);
					// add the remaining points
					for (std::size_t i = 1; i < segmentCount; ++i) {
						Real angle = i * segmentSubAngle + startAngle;
						result.addVertex(arc.radius * std::cos(angle) + arc.center.x(),
							arc.radius * std::sin(angle) + arc.center.y(), 0);
					}
				}

				return true;
			};

			iterateSegIndices(pline, visitor);
			if (!pline.isClosed()) {
				result.addVertex(pline.lastVertex());
			}

			return result;
		}

		/// Computes a fast approximate AABB of a segment described by v1 to v2, bounding box may be larger
		/// than the true bounding box for the segment
		template <typename Real>
		AABB<Real> createFastApproxBoundingBox(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
			AABB<Real> result;
			if (v1.bulgeIsZero()) {
				if (v1.x() < v2.x()) {
					result.xMin = v1.x();
					result.xMax = v2.x();
				}
				else {
					result.xMin = v2.x();
					result.xMax = v1.x();
				}

				if (v1.y() < v2.y()) {
					result.yMin = v1.y();
					result.yMax = v2.y();
				}
				else {
					result.yMin = v2.y();
					result.yMax = v1.y();
				}

				return result;
			}

			// For arcs we don't compute the actual extents which requires slow trig functions, instead we
			// create an approximate bounding box from the rectangle formed by extending the chord by the
			// sagitta, NOTE: this approximate bounding box is always equal to or bigger than the true
			// bounding box
			Real b = v1.bulge();
			Real offsX = b * (v2.y() - v1.y()) / Real(2);
			Real offsY = -b * (v2.x() - v1.x()) / Real(2);

			Real pt1X = v1.x() + offsX;
			Real pt2X = v2.x() + offsX;
			Real pt1Y = v1.y() + offsY;
			Real pt2Y = v2.y() + offsY;

			Real endPointXMin, endPointXMax;
			if (v1.x() < v2.x()) {
				endPointXMin = v1.x();
				endPointXMax = v2.x();
			}
			else {
				endPointXMin = v2.x();
				endPointXMax = v1.x();
			}

			Real ptXMin, ptXMax;
			if (pt1X < pt2X) {
				ptXMin = pt1X;
				ptXMax = pt2X;
			}
			else {
				ptXMin = pt2X;
				ptXMax = pt1X;
			}

			Real endPointYMin, endPointYMax;
			if (v1.y() < v2.y()) {
				endPointYMin = v1.y();
				endPointYMax = v2.y();
			}
			else {
				endPointYMin = v2.y();
				endPointYMax = v1.y();
			}

			Real ptYMin, ptYMax;
			if (pt1Y < pt2Y) {
				ptYMin = pt1Y;
				ptYMax = pt2Y;
			}
			else {
				ptYMin = pt2Y;
				ptYMax = pt1Y;
			}

			result.xMin = std::min(endPointXMin, ptXMin);
			result.yMin = std::min(endPointYMin, ptYMin);
			result.xMax = std::max(endPointXMax, ptXMax);
			result.yMax = std::max(endPointYMax, ptYMax);
			return result;
		}

		/// Creates an approximate spatial index for all the segments in the polyline given using
		/// createFastApproxBoundingBox.
		template <typename Real>
		StaticSpatialIndex<Real> createApproxSpatialIndex(Polyline<Real> const &pline) {
			assert(pline.size() > 1 && "need at least 2 vertexes to form segments for spatial index");

			std::size_t segmentCount = pline.isClosed() ? pline.size() : pline.size() - 1;
			StaticSpatialIndex<Real> result(segmentCount);

			for (std::size_t i = 0; i < pline.size() - 1; ++i) {
				AABB<Real> approxBB = createFastApproxBoundingBox(pline[i], pline[i + 1]);
				result.add(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax);
			}

			if (pline.isClosed()) {
				// add final segment from last to first
				AABB<Real> approxBB = createFastApproxBoundingBox(pline.lastVertex(), pline[0]);
				result.add(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax);
			}

			result.finish();

			return result;
		}

		/// Inverts the direction of the polyline given. If polyline is closed then this just changes the
		/// direction from clockwise to counter clockwise, if polyline is open then the starting vertex
		/// becomes the end vertex and the end vertex becomes the starting vertex.
		template <typename Real> void invertDirection(Polyline<Real> &pline) {
			if (pline.size() < 2) {
				return;
			}
			std::reverse(std::begin(pline.vertexes()), std::end(pline.vertexes()));

			// shift and negate bulge (to maintain same geometric path)
			Real firstBulge = pline[0].bulge();

			for (std::size_t i = 1; i < pline.size(); ++i) {
				pline[i - 1].bulge() = -pline[i].bulge();
			}

			pline.lastVertex().bulge() = -firstBulge;
		}

		/// Returns a new polyline with all singularities (repeating vertex positions) from the polyline
		/// given removed.
		template <typename Real>
		Polyline<Real> pruneSingularities(Polyline<Real> const &pline,
			Real epsilon = utils::realPrecision<Real>()) {
			Polyline<Real> result;
			result.isClosed() = pline.isClosed();

			if (pline.size() == 0) {
				return result;
			}

			// allocate up front (most of the time the number of repeated positions are much less than the
			// total number of vertexes so we're not using very much more memory than required)
			result.vertexes().reserve(pline.size());

			result.addVertex(pline[0]);

			for (std::size_t i = 1; i < pline.size(); ++i) {
				if (fuzzyEqual(result.lastVertex().pos(), pline[i].pos(), epsilon)) {
					result.lastVertex().bulge() = pline[i].bulge();
				}
				else {
					result.addVertex(pline[i]);
				}
			}

			if (result.isClosed() && result.size() > 1) {
				if (fuzzyEqual(result.lastVertex().pos(), result[0].pos(), epsilon)) {
					result.vertexes().pop_back();
				}
			}

			return result;
		}

		/// Compute the area of a closed polyline, assumes no self intersects, returns positive number if
		/// polyline direction is counter clockwise, negative if clockwise, zero if not closed
		template <typename Real> Real area(Polyline<Real> const &pline) {
			// Implementation notes:
			// Using the shoelace formula (https://en.wikipedia.org/wiki/Shoelace_formula) modified to support
			// arcs defined by a bulge value. The shoelace formula returns a negative value for clockwise
			// oriented polygons and positive value for counter clockwise oriented polygons. The area of each
			// circular segment defined by arcs is then added if it is a counter clockwise arc or subtracted
			// if it is a clockwise arc. The area of the circular segments are computed by finding the area of
			// the arc sector minus the area of the triangle defined by the chord and center of circle.
			// See https://en.wikipedia.org/wiki/Circular_segment
			if (!pline.isClosed() || pline.size() < 2) {
				return Real(0);
			}

			Real doubleEdgeAreaTotal = Real(0);
			Real doubleArcSegAreaTotal = Real(0);

			auto visitor = [&](std::size_t i, std::size_t j) {
				doubleEdgeAreaTotal += pline[i].x() * pline[j].y() - pline[i].y() * pline[j].x();
				if (!pline[i].bulgeIsZero()) {
					// add segment area
					Real b = std::abs(pline[i].bulge());
					Real sweepAngle = Real(4) * std::atan(b);
					Real triangleBase = length(pline[j].pos() - pline[i].pos());
					Real radius = triangleBase * (b * b + Real(1)) / (Real(4) * b);
					Real sagitta = b * triangleBase / Real(2);
					Real triangleHeight = radius - sagitta;
					Real doubleSectorArea = sweepAngle * radius * radius;
					Real doubleTriangleArea = triangleBase * triangleHeight;
					Real doubleArcSegArea = doubleSectorArea - doubleTriangleArea;
					if (pline[i].bulge() < Real(0)) {
						doubleArcSegArea = -doubleArcSegArea;
					}

					doubleArcSegAreaTotal += doubleArcSegArea;
				}

				// iterate all segments
				return true;
			};

			iterateSegIndices(pline, visitor);

			return (doubleEdgeAreaTotal + doubleArcSegAreaTotal) / Real(2);
		}

		/// Compute the winding number for the point in relation to the polyline. If polyline is open and
		/// the first vertex does not overlap the last vertex then 0 is always returned. This algorithm is
		/// adapted from http://geomalgorithms.com/a03-_inclusion.html to support arc segments. NOTE: The
		/// result is not defined if the point lies ontop of the polyline.
		template <typename Real>
		int windingNumber(Polyline<Real> const &pline, Vector2<Real> const &point) {
			if (!pline.isClosed() && !(fuzzyEqual(pline[0].pos(), pline.lastVertex().pos()))) {
				return 0;
			}

			int windingNumber = 0;

			auto lineVisitor = [&](const auto &v1, const auto &v2) {
				if (v1.y() <= point.y()) {
					if (v2.y() > point.y() && isLeft(v1.pos(), v2.pos(), point))
						// left and upward crossing
						windingNumber += 1;
				}
				else if (v2.y() <= point.y() && !(isLeft(v1.pos(), v2.pos(), point))) {
					// right and downward crossing
					windingNumber -= 1;
				}
			};

			// Helper function to determine if point is inside an arc sector area
			auto distToArcCenterLessThanRadius = [](const auto &v1, const auto &v2, const auto &pt) {
				auto arc = arcRadiusAndCenter(v1, v2);
				Real dist2 = distSquared(arc.center, pt);
				return dist2 < arc.radius * arc.radius;
			};

			auto arcVisitor = [&](const auto &v1, const auto &v2) {
				bool isCCW = v1.bulge() > Real(0);
				bool pointIsLeft = isLeft(v1.pos(), v2.pos(), point);

				if (v1.y() <= point.y()) {
					if (v2.y() > point.y()) {
						// upward crossing of arc chord
						if (isCCW) {
							if (pointIsLeft) {
								// counter clockwise arc left of chord
								windingNumber += 1;
							}
							else {
								// counter clockwise arc right of chord
								if (distToArcCenterLessThanRadius(v1, v2, point)) {
									windingNumber += 1;
								}
							}
						}
						else {
							if (pointIsLeft) {
								// clockwise arc left of chord
								if (!distToArcCenterLessThanRadius(v1, v2, point)) {
									windingNumber += 1;
								}
							}
							// else clockwise arc right of chord, no crossing
						}
					}
					else {
						// not crossing arc chord and chord is below, check if point is inside arc sector
						if (isCCW && !pointIsLeft) {
							if (v2.x() < point.x() && point.x() < v1.x() &&
								distToArcCenterLessThanRadius(v1, v2, point)) {
								windingNumber += 1;
							}
						}
						else if (!isCCW && pointIsLeft) {
							if (v1.x() < point.x() && point.x() < v2.x() &&
								distToArcCenterLessThanRadius(v1, v2, point)) {
								windingNumber -= 1;
							}
						}
					}
				}
				else {
					if (v2.y() <= point.y()) {
						// downward crossing of arc chord
						if (isCCW) {
							if (!pointIsLeft) {
								// counter clockwise arc right of chord
								if (!distToArcCenterLessThanRadius(v1, v2, point)) {
									windingNumber -= 1;
								}
							}
							// else counter clockwise arc left of chord, no crossing
						}
						else {
							if (pointIsLeft) {
								// clockwise arc left of chord
								if (distToArcCenterLessThanRadius(v1, v2, point)) {
									windingNumber -= 1;
								}
							}
							else {
								// clockwise arc right of chord
								windingNumber -= 1;
							}
						}
					}
					else {
						// not crossing arc chord and chord is above, check if point is inside arc sector
						if (isCCW && !pointIsLeft) {
							if (v1.x() < point.x() && point.x() < v2.x() &&
								distToArcCenterLessThanRadius(v1, v2, point)) {
								windingNumber += 1;
							}
						}
						else if (!isCCW && pointIsLeft) {
							if (v2.x() < point.x() && point.x() < v1.x() &&
								distToArcCenterLessThanRadius(v1, v2, point)) {
								windingNumber -= 1;
							}
						}
					}
				}
			};

			auto visitor = [&](std::size_t i, std::size_t j) {
				const auto &v1 = pline[i];
				const auto &v2 = pline[j];
				if (v1.bulgeIsZero()) {
					lineVisitor(v1, v2);
				}
				else {
					arcVisitor(v1, v2);
				}

				return true;
			};

			iterateSegIndices(pline, visitor);

			return windingNumber;
		}

		/// Represents a raw polyline offset segment.
		template <typename Real> struct PlineOffsetSegment {
			PlineVertex<Real> v1;
			PlineVertex<Real> v2;
			Vector2<Real> origV2Pos;
			bool collapsedArc;
		};

		/// Creates all the raw polyline offset segments.
		template <typename Real>
		std::vector<PlineOffsetSegment<Real>> createUntrimmedOffsetSegments(Polyline<Real> const &pline,
			Real offset) {
			std::size_t segmentCount = pline.isClosed() ? pline.size() : pline.size() - 1;

			std::vector<PlineOffsetSegment<Real>> result;
			result.reserve(segmentCount);

			auto lineVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
				result.emplace_back();
				PlineOffsetSegment<Real> &seg = result.back();
				seg.collapsedArc = false;
				seg.origV2Pos = v2.pos();
				Vector2<Real> edge = v2.pos() - v1.pos();
				Vector2<Real> offsetV = offset * unitPerp(edge);
				seg.v1.pos() = v1.pos() + offsetV;
				seg.v1.bulge() = v1.bulge();
				seg.v2.pos() = v2.pos() + offsetV;
				seg.v2.bulge() = v2.bulge();
			};

			auto arcVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
				auto arc = arcRadiusAndCenter(v1, v2);
				Real offs = v1.bulge() < Real(0) ? offset : -offset;
				Real radiusAfterOffset = arc.radius + offs;
				Vector2<Real> v1ToCenter = v1.pos() - arc.center;
				normalize(v1ToCenter);
				Vector2<Real> v2ToCenter = v2.pos() - arc.center;
				normalize(v2ToCenter);

				result.emplace_back();
				PlineOffsetSegment<Real> &seg = result.back();
				seg.origV2Pos = v2.pos();
				seg.v1.pos() = offs * v1ToCenter + v1.pos();
				seg.v2.pos() = offs * v2ToCenter + v2.pos();
				seg.v2.bulge() = v2.bulge();

				if (radiusAfterOffset < utils::realThreshold<Real>()) {
					// collapsed arc, offset arc start and end points towards arc center and turn into line
					// handles case where offset vertexes are equal and simplifies path for clipping algorithm
					seg.collapsedArc = true;
					seg.v1.bulge() = Real(0);
				}
				else {
					seg.collapsedArc = false;
					seg.v1.bulge() = v1.bulge();
				}
			};

			auto offsetVisitor = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
				if (v1.bulgeIsZero()) {
					lineVisitor(v1, v2);
				}
				else {
					arcVisitor(v1, v2);
				}
			};

			for (std::size_t i = 1; i < pline.size(); ++i) {
				offsetVisitor(pline[i - 1], pline[i]);
			}

			if (pline.isClosed()) {
				offsetVisitor(pline.lastVertex(), pline[0]);
			}

			return result;
		}

		namespace detail {
			struct IndexPairHash {
				std::size_t operator()(std::pair<std::size_t, std::size_t> const &pair) const {
					return std::hash<std::size_t>()(pair.first) ^ std::hash<std::size_t>()(pair.second);
				}
			};

			template <typename Real>
			Vector2<Real> arcTangentVector(Vector2<Real> const &arcCenter, bool isCCW,
				Vector2<Real> const &pointOnArc) {
				if (isCCW) {
					// counter clockwise, rotate vector from center to pointOnArc 90 degrees
					return Vector2<Real>(arcCenter.y() - pointOnArc.y(), pointOnArc.x() - arcCenter.x());
				}

				// clockwise, rotate vector from center to pointOnArc -90 degrees
				return Vector2<Real>(pointOnArc.y() - arcCenter.y(), arcCenter.x() - pointOnArc.x());
			}

			template <typename Real> bool falseIntersect(Real t) { return t < 0.0 || t > 1.0; }

			template <typename Real> Real segX(Vector2<Real> const &p0, Vector2<Real> const &p1, Real t) {
				return p0.x() + t * (p1.x() - p0.x());
			}
			template <typename Real> Real segY(Vector2<Real> const &p0, Vector2<Real> const &p1, Real t) {
				return p0.y() + t * (p1.y() - p0.y());
			}

			template <typename Real>
			Vector2<Real> pointFromParametric(Vector2<Real> const &p0, Vector2<Real> const &p1, Real t) {
				return Vector2<Real>(segX(p0, p1, t), segY(p0, p1, t));
			}

			template <typename Real>
			Vector2<Real> segMidpoint(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
				if (v1.bulgeIsZero()) {
					return midpoint(v1.pos(), v2.pos());
				}

				auto arc = arcRadiusAndCenter(v1, v2);
				Real a1 = angle(arc.center, v1.pos());
				Real a2 = angle(arc.center, v2.pos());
				Real angleOffset = std::abs(utils::deltaAngle(a1, a2) / Real(2));
				// use arc direction to determine offset sign to robustly handle half circles
				Real midAngle = v1.bulge() > Real(0) ? a1 + angleOffset : a1 - angleOffset;
				return pointOnCircle(arc.radius, arc.center, midAngle);
			}

			template <typename Real>
			void addOrReplaceIfSamePos(Polyline<Real> &pline, PlineVertex<Real> const &vertex,
				Real epsilon = utils::realPrecision<Real>()) {
				if (pline.size() == 0) {
					pline.addVertex(vertex);
					return;
				}

				if (fuzzyEqual(pline.lastVertex().pos(), vertex.pos(), epsilon)) {
					pline.lastVertex().bulge() = vertex.bulge();
					return;
				}

				pline.addVertex(vertex);
			}
			// Gets the bulge to describe the arc going from start point to end point with the given arc center
			// and curve orientation, if orientation is negative then bulge is negative otherwise it is positive
			template <typename Real>
			Real bulgeForConnection(Vector2<Real> const &arcCenter, Vector2<Real> const &sp,
				Vector2<Real> const &ep, bool isCCW) {
				Real a1 = angle(arcCenter, sp);
				Real a2 = angle(arcCenter, ep);
				Real absSweepAngle = std::abs(utils::deltaAngle(a1, a2));
				Real absBulge = std::tan(absSweepAngle / Real(4));
				if (isCCW) {
					return absBulge;
				}

				return -absBulge;
			}

			template <typename Real>
			void lineToLineJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
				bool connectionArcsAreCCW, Polyline<Real> &result) {
				const auto &v1 = s1.v1;
				const auto &v2 = s1.v2;
				const auto &u1 = s2.v1;
				const auto &u2 = s2.v2;
				assert(v1.bulgeIsZero() && u1.bulgeIsZero() && "both pairs should be lines");

				auto connectUsingArc = [&] {
					auto const &arcCenter = s1.origV2Pos;
					auto const &sp = v2.pos();
					auto const &ep = u1.pos();
					Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
					addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
					addOrReplaceIfSamePos(result, PlineVertex<Real>(ep, Real(0)));
				};

				if (s1.collapsedArc || s2.collapsedArc) {
					// connecting to/from collapsed arc, always connect using arc
					connectUsingArc();
				}
				else {
					auto intrResult = intrLineSeg2LineSeg2(v1.pos(), v2.pos(), u1.pos(), u2.pos());

					switch (intrResult.intrType) {
					case LineSeg2LineSeg2IntrType::None:
						// just join with straight line
						addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
						addOrReplaceIfSamePos(result, u1);
						break;
					case LineSeg2LineSeg2IntrType::True:
						addOrReplaceIfSamePos(result, PlineVertex<Real>(intrResult.point, Real(0)));
						break;
					case LineSeg2LineSeg2IntrType::Coincident:
						addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
						break;
					case LineSeg2LineSeg2IntrType::False:
						if (intrResult.t0 > Real(1) && falseIntersect(intrResult.t1)) {
							// extend and join the lines together using an arc
							connectUsingArc();
						}
						else {
							addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
							addOrReplaceIfSamePos(result, u1);
						}
						break;
					}
				}
			}

			template <typename Real>
			void lineToArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
				bool connectionArcsAreCCW, Polyline<Real> &result) {

				const auto &v1 = s1.v1;
				const auto &v2 = s1.v2;
				const auto &u1 = s2.v1;
				const auto &u2 = s2.v2;
				assert(v1.bulgeIsZero() && !u1.bulgeIsZero() &&
					"first pair should be arc, second pair should be line");

				auto connectUsingArc = [&] {
					auto const &arcCenter = s1.origV2Pos;
					auto const &sp = v2.pos();
					auto const &ep = u1.pos();
					Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
					addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
					addOrReplaceIfSamePos(result, u1);
				};

				const auto arc = arcRadiusAndCenter(u1, u2);

				auto processIntersect = [&](Real t, Vector2<Real> const &intersect) {
					const bool trueSegIntersect = !falseIntersect(t);
					const bool trueArcIntersect =
						pointWithinArcSweepAngle(arc.center, u1.pos(), u2.pos(), u1.bulge(), intersect);
					if (trueSegIntersect && trueArcIntersect) {
						// trim at intersect
						Real a = angle(arc.center, intersect);
						Real arcEndAngle = angle(arc.center, u2.pos());
						Real theta = utils::deltaAngle(a, arcEndAngle);
						// ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
						// which case we do not want to update the bulge)
						if ((theta > Real(0)) == (u1.bulge() > Real(0))) {
							addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, std::tan(theta / Real(4))));
						}
						else {
							addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, u1.bulge()));
						}
					}
					else if (t > Real(1) && !trueArcIntersect) {
						connectUsingArc();
					}
					else if (s1.collapsedArc) {
						// collapsed arc connecting to arc, connect using arc
						connectUsingArc();
					}
					else {
						// connect using line
						addOrReplaceIfSamePos(result, PlineVertex<Real>(v2.pos(), Real(0)));
						addOrReplaceIfSamePos(result, u1);
					}
				};

				auto intrResult = intrLineSeg2Circle2(v1.pos(), v2.pos(), arc.radius, arc.center);
				if (intrResult.numIntersects == 0) {
					connectUsingArc();
				}
				else if (intrResult.numIntersects == 1) {
					processIntersect(intrResult.t0, pointFromParametric(v1.pos(), v2.pos(), intrResult.t0));
				}
				else {
					assert(intrResult.numIntersects == 2);
					// always use intersect closest to original point
					Vector2<Real> i1 = pointFromParametric(v1.pos(), v2.pos(), intrResult.t0);
					Real dist1 = distSquared(i1, s1.origV2Pos);
					Vector2<Real> i2 = pointFromParametric(v1.pos(), v2.pos(), intrResult.t1);
					Real dist2 = distSquared(i2, s1.origV2Pos);

					if (dist1 < dist2) {
						processIntersect(intrResult.t0, i1);
					}
					else {
						processIntersect(intrResult.t1, i2);
					}
				}
			}

			template <typename Real>
			void arcToLineJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
				bool connectionArcsAreCCW, Polyline<Real> &result) {

				const auto &v1 = s1.v1;
				const auto &v2 = s1.v2;
				const auto &u1 = s2.v1;
				const auto &u2 = s2.v2;
				assert(!v1.bulgeIsZero() && u1.bulgeIsZero() &&
					"first pair should be line, second pair should be arc");

				auto connectUsingArc = [&] {
					auto const &arcCenter = s1.origV2Pos;
					auto const &sp = v2.pos();
					auto const &ep = u1.pos();
					Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
					addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
					addOrReplaceIfSamePos(result, u1);
				};

				const auto arc = arcRadiusAndCenter(v1, v2);

				auto processIntersect = [&](Real t, Vector2<Real> const &intersect) {
					const bool trueSegIntersect = !falseIntersect(t);
					const bool trueArcIntersect =
						pointWithinArcSweepAngle(arc.center, v1.pos(), v2.pos(), v1.bulge(), intersect);
					if (trueSegIntersect && trueArcIntersect) {
						// modify previous bulge and trim at intersect
						PlineVertex<Real> &prevVertex = result.lastVertex();

						if (!prevVertex.bulgeIsZero()) {
							Real a = angle(arc.center, intersect);
							auto prevArc = arcRadiusAndCenter(prevVertex, v2);
							Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
							Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, a);

							// ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
							// which case we do not want to update the bulge)
							if ((updatedPrevTheta > Real(0)) == (prevVertex.bulge() > Real(0))) {
								result.lastVertex().bulge() = std::tan(updatedPrevTheta / Real(4));
							}
						}

						addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, Real(0)));

					}
					else {
						connectUsingArc();
					}
				};

				auto intrResult = intrLineSeg2Circle2(u1.pos(), u2.pos(), arc.radius, arc.center);
				if (intrResult.numIntersects == 0) {
					connectUsingArc();
				}
				else if (intrResult.numIntersects == 1) {
					processIntersect(intrResult.t0, pointFromParametric(u1.pos(), u2.pos(), intrResult.t0));
				}
				else {
					assert(intrResult.numIntersects == 2);
					const auto &origPoint = s2.collapsedArc ? u1.pos() : s1.origV2Pos;
					Vector2<Real> i1 = pointFromParametric(u1.pos(), u2.pos(), intrResult.t0);
					Real dist1 = distSquared(i1, origPoint);
					Vector2<Real> i2 = pointFromParametric(u1.pos(), u2.pos(), intrResult.t1);
					Real dist2 = distSquared(i2, origPoint);

					if (dist1 < dist2) {
						processIntersect(intrResult.t0, i1);
					}
					else {
						processIntersect(intrResult.t1, i2);
					}
				}
			}

			template <typename Real>
			void arcToArcJoin(PlineOffsetSegment<Real> const &s1, PlineOffsetSegment<Real> const &s2,
				bool connectionArcsAreCCW, Polyline<Real> &result) {

				const auto &v1 = s1.v1;
				const auto &v2 = s1.v2;
				const auto &u1 = s2.v1;
				const auto &u2 = s2.v2;
				assert(!v1.bulgeIsZero() && !u1.bulgeIsZero() && "both pairs should be arcs");

				const auto arc1 = arcRadiusAndCenter(v1, v2);
				const auto arc2 = arcRadiusAndCenter(u1, u2);

				auto connectUsingArc = [&] {
					auto const &arcCenter = s1.origV2Pos;
					auto const &sp = v2.pos();
					auto const &ep = u1.pos();
					Real bulge = bulgeForConnection(arcCenter, sp, ep, connectionArcsAreCCW);
					addOrReplaceIfSamePos(result, PlineVertex<Real>(sp, bulge));
					addOrReplaceIfSamePos(result, u1);
				};

				auto processIntersect = [&](Vector2<Real> const &intersect) {
					const bool trueArcIntersect1 =
						pointWithinArcSweepAngle(arc1.center, v1.pos(), v2.pos(), v1.bulge(), intersect);
					const bool trueArcIntersect2 =
						pointWithinArcSweepAngle(arc2.center, u1.pos(), u2.pos(), u1.bulge(), intersect);

					if (trueArcIntersect1 && trueArcIntersect2) {
						// modify previous bulge and trim at intersect
						PlineVertex<Real> &prevVertex = result.lastVertex();
						if (!prevVertex.bulgeIsZero()) {
							Real a1 = angle(arc1.center, intersect);
							auto prevArc = arcRadiusAndCenter(prevVertex, v2);
							Real prevArcStartAngle = angle(prevArc.center, prevVertex.pos());
							Real updatedPrevTheta = utils::deltaAngle(prevArcStartAngle, a1);

							// ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
							// which case we do not want to update the bulge)
							if ((updatedPrevTheta > Real(0)) == (prevVertex.bulge() > Real(0))) {
								result.lastVertex().bulge() = std::tan(updatedPrevTheta / Real(4));
							}
						}

						// add the vertex at our current trim/join point
						Real a2 = angle(arc2.center, intersect);
						Real endAngle = angle(arc2.center, u2.pos());
						Real theta = utils::deltaAngle(a2, endAngle);

						// ensure the sign matches (may get flipped if intersect is at the very end of the arc, in
						// which case we do not want to update the bulge)
						if ((theta > Real(0)) == (u1.bulge() > Real(0))) {
							addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, std::tan(theta / Real(4))));
						}
						else {
							addOrReplaceIfSamePos(result, PlineVertex<Real>(intersect, u1.bulge()));
						}

					}
					else {
						connectUsingArc();
					}
				};

				const auto intrResult = intrCircle2Circle2(arc1.radius, arc1.center, arc2.radius, arc2.center);
				switch (intrResult.intrType) {
				case Circle2Circle2IntrType::NoIntersect:
					connectUsingArc();
					break;
				case Circle2Circle2IntrType::OneIntersect:
					processIntersect(intrResult.point1);
					break;
				case Circle2Circle2IntrType::TwoIntersects: {
					Real dist1 = distSquared(intrResult.point1, s1.origV2Pos);
					Real dist2 = distSquared(intrResult.point2, s1.origV2Pos);
					if (dist1 < dist2) {
						processIntersect(intrResult.point1);
					}
					else {
						processIntersect(intrResult.point2);
					}
				} break;
				case Circle2Circle2IntrType::Coincident:
					// same constant arc radius and center, just add the vertex (nothing to trim/extend)
					addOrReplaceIfSamePos(result, u1);
					break;
				}
			}

			enum class PlineSegIntrType {
				NoIntersect,
				TangentIntersect,
				OneIntersect,
				TwoIntersects,
				SegmentOverlap,
				ArcOverlap
			};

			template <typename Real> struct IntrPlineSegsResult {
				PlineSegIntrType intrType;
				Vector2<Real> point1;
				Vector2<Real> point2;
			};

			template <typename Real>
			IntrPlineSegsResult<Real> intrPlineSegs(PlineVertex<Real> const &v1, PlineVertex<Real> const &v2,
				PlineVertex<Real> const &u1, PlineVertex<Real> const &u2) {
				IntrPlineSegsResult<Real> result;
				const bool vIsLine = v1.bulgeIsZero();
				const bool uIsLine = u1.bulgeIsZero();

				// helper function to process line arc intersect
				auto processLineArcIntr = [&result](Vector2<Real> const &p0, Vector2<Real> const &p1,
					PlineVertex<Real> const &a1, PlineVertex<Real> const &a2) {
					auto arc = arcRadiusAndCenter(a1, a2);
					auto intrResult = intrLineSeg2Circle2(p0, p1, arc.radius, arc.center);

					// helper function to test and get point within arc sweep
					auto pointInSweep = [&](Real t) {
						if (t + utils::realThreshold<Real>() < Real(0) ||
							t > Real(1) + utils::realThreshold<Real>()) {
							return std::make_pair(false, Vector2<Real>());
						}

						Vector2<Real> p = pointFromParametric(p0, p1, t);
						bool withinSweep = pointWithinArcSweepAngle(arc.center, a1.pos(), a2.pos(), a1.bulge(), p);
						return std::make_pair(withinSweep, p);
					};

					if (intrResult.numIntersects == 0) {
						result.intrType = PlineSegIntrType::NoIntersect;
					}
					else if (intrResult.numIntersects == 1) {
						auto p = pointInSweep(intrResult.t0);
						if (p.first) {
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = p.second;
						}
						else {
							result.intrType = PlineSegIntrType::NoIntersect;
						}
					}
					else {
						assert(intrResult.numIntersects == 2);
						auto p1 = pointInSweep(intrResult.t0);
						auto p2 = pointInSweep(intrResult.t1);

						if (p1.first && p2.first) {
							result.intrType = PlineSegIntrType::TwoIntersects;
							result.point1 = p1.second;
							result.point2 = p2.second;
						}
						else if (p1.first) {
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = p1.second;
						}
						else if (p2.first) {
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = p2.second;
						}
						else {
							result.intrType = PlineSegIntrType::NoIntersect;
						}
					}
				};

				if (vIsLine && uIsLine) {
					auto intrResult = intrLineSeg2LineSeg2(v1.pos(), v2.pos(), u1.pos(), u2.pos());
					switch (intrResult.intrType) {
					case LineSeg2LineSeg2IntrType::None:
						result.intrType = PlineSegIntrType::NoIntersect;
						break;
					case LineSeg2LineSeg2IntrType::True:
						result.intrType = PlineSegIntrType::OneIntersect;
						result.point1 = intrResult.point;
						break;
					case LineSeg2LineSeg2IntrType::Coincident:
						result.intrType = PlineSegIntrType::SegmentOverlap;
						// build points from parametric parameters for v1->v2
						result.point1 = pointFromParametric(v1.pos(), v2.pos(), intrResult.t0);
						result.point2 = pointFromParametric(v1.pos(), v2.pos(), intrResult.t1);
						break;
					case LineSeg2LineSeg2IntrType::False:
						result.intrType = PlineSegIntrType::NoIntersect;
						break;
					}

				}
				else if (vIsLine) {
					processLineArcIntr(v1.pos(), v2.pos(), u1, u2);
				}
				else if (uIsLine) {
					processLineArcIntr(u1.pos(), u2.pos(), v1, v2);
				}
				else {
					auto arc1 = arcRadiusAndCenter(v1, v2);
					auto arc2 = arcRadiusAndCenter(u1, u2);

					auto startAndSweepAngle = [](Vector2<Real> const &sp, Vector2<Real> const &center, Real bulge) {
						Real startAngle = utils::normalizeRadians(angle(center, sp));
						Real sweepAngle = Real(4) * std::atan(bulge);
						return std::make_pair(startAngle, sweepAngle);
					};

					auto bothArcsSweepPoint = [&](Vector2<Real> const &pt) {
						return pointWithinArcSweepAngle(arc1.center, v1.pos(), v2.pos(), v1.bulge(), pt) &&
							pointWithinArcSweepAngle(arc2.center, u1.pos(), u2.pos(), u1.bulge(), pt);
					};

					auto intrResult = intrCircle2Circle2(arc1.radius, arc1.center, arc2.radius, arc2.center);

					switch (intrResult.intrType) {
					case Circle2Circle2IntrType::NoIntersect:
						result.intrType = PlineSegIntrType::NoIntersect;
						break;
					case Circle2Circle2IntrType::OneIntersect:
						if (bothArcsSweepPoint(intrResult.point1)) {
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = intrResult.point1;
						}
						else {
							result.intrType = PlineSegIntrType::NoIntersect;
						}
						break;
					case Circle2Circle2IntrType::TwoIntersects: {
						const bool pt1InSweep = bothArcsSweepPoint(intrResult.point1);
						const bool pt2InSweep = bothArcsSweepPoint(intrResult.point2);
						if (pt1InSweep && pt2InSweep) {
							result.intrType = PlineSegIntrType::TwoIntersects;
							result.point1 = intrResult.point1;
							result.point2 = intrResult.point2;
						}
						else if (pt1InSweep) {
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = intrResult.point1;
						}
						else if (pt2InSweep) {
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = intrResult.point2;
						}
						else {
							result.intrType = PlineSegIntrType::NoIntersect;
						}
					} break;
					case Circle2Circle2IntrType::Coincident:
						// determine if arcs overlap along their sweep
						// start and sweep angles
						auto arc1StartAndSweep = startAndSweepAngle(v1.pos(), arc1.center, v1.bulge());
						auto arc2StartAndSweep = startAndSweepAngle(u1.pos(), arc2.center, u1.bulge());
						// end angles (start + sweep)
						auto arc1End = arc1StartAndSweep.first + arc1StartAndSweep.second;
						auto arc2End = arc2StartAndSweep.first + arc2StartAndSweep.second;

						if (utils::fuzzyEqual(arc1StartAndSweep.first, arc2End)) {
							// only end points touch at start of arc1
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = v1.pos();
						}
						else if (utils::fuzzyEqual(arc2StartAndSweep.first, arc1End)) {
							// only end points touch at start of arc2
							result.intrType = PlineSegIntrType::OneIntersect;
							result.point1 = u1.pos();
						}
						else {
							const bool arc2StartsInArc1Sweep = utils::angleIsWithinSweep(
								arc1StartAndSweep.first, arc1StartAndSweep.second, arc2StartAndSweep.first);
							const bool arc2EndsInArc1Sweep =
								utils::angleIsWithinSweep(arc1StartAndSweep.first, arc1StartAndSweep.second, arc2End);
							if (arc2StartsInArc1Sweep && arc2EndsInArc1Sweep) {
								// arc2 is fully overlapped by arc1
								result.intrType = PlineSegIntrType::ArcOverlap;
								result.point1 = u1.pos();
								result.point2 = u2.pos();
							}
							else if (arc2StartsInArc1Sweep) {
								// overlap from arc2 start to arc1 end
								result.intrType = PlineSegIntrType::ArcOverlap;
								result.point1 = u1.pos();
								result.point2 = v2.pos();
							}
							else if (arc2EndsInArc1Sweep) {
								// overlap from arc1 start to arc2 end
								result.intrType = PlineSegIntrType::ArcOverlap;
								result.point1 = v1.pos();
								result.point2 = u2.pos();
							}
							else {
								const bool arc1StartsInArc2Sweep = utils::angleIsWithinSweep(
									arc2StartAndSweep.first, arc2StartAndSweep.second, arc1StartAndSweep.first);
								if (arc1StartsInArc2Sweep) {
									result.intrType = PlineSegIntrType::ArcOverlap;
									result.point1 = v1.pos();
									result.point2 = v2.pos();
								}
								else {
									result.intrType = PlineSegIntrType::NoIntersect;
								}
							}
						}

						break;
					}
				}

				return result;
			}

			template <typename Real>
			void offsetCircleIntersectsWithPline(Polyline<Real> const &pline, Real offset,
				Vector2<Real> const &circleCenter,
				StaticSpatialIndex<Real> const &spatialIndex,
				std::vector<std::pair<std::size_t, Vector2<Real>>> &output) {

				const Real circleRadius = std::abs(offset);

				std::vector<std::size_t> queryResults;

				spatialIndex.query(circleCenter.x() - circleRadius, circleCenter.y() - circleRadius,
					circleCenter.x() + circleRadius, circleCenter.y() + circleRadius,
					queryResults);

				auto validLineSegIntersect = [](Real t) {
					return !falseIntersect(t) && std::abs(t) > utils::realPrecision<Real>();
				};

				auto validArcSegIntersect = [](Vector2<Real> const &arcCenter, Vector2<Real> const &arcStart,
					Vector2<Real> const &arcEnd, Real bulge,
					Vector2<Real> const &intrPoint) {
					return !fuzzyEqual(arcStart, intrPoint, utils::realPrecision<Real>()) &&
						pointWithinArcSweepAngle(arcCenter, arcStart, arcEnd, bulge, intrPoint);
				};

				for (std::size_t sIndex : queryResults) {
					PlineVertex<Real> const &v1 = pline[sIndex];
					PlineVertex<Real> const &v2 = pline[sIndex + 1];
					if (v1.bulgeIsZero()) {
						IntrLineSeg2Circle2Result<Real> intrResult =
							intrLineSeg2Circle2(v1.pos(), v2.pos(), circleRadius, circleCenter);
						if (intrResult.numIntersects == 0) {
							continue;
						}
						else if (intrResult.numIntersects == 1) {
							if (validLineSegIntersect(intrResult.t0)) {
								output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), intrResult.t0));
							}
						}
						else {
							assert(intrResult.numIntersects == 2);
							if (validLineSegIntersect(intrResult.t0)) {
								output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), intrResult.t0));
							}
							if (validLineSegIntersect(intrResult.t1)) {
								output.emplace_back(sIndex, pointFromParametric(v1.pos(), v2.pos(), intrResult.t1));
							}
						}
					}
					else {
						auto arc = arcRadiusAndCenter(v1, v2);
						IntrCircle2Circle2Result<Real> intrResult =
							intrCircle2Circle2(arc.radius, arc.center, circleRadius, circleCenter);
						switch (intrResult.intrType) {
						case Circle2Circle2IntrType::NoIntersect:
							break;
						case Circle2Circle2IntrType::OneIntersect:
							if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intrResult.point1)) {
								output.emplace_back(sIndex, intrResult.point1);
							}
							break;
						case Circle2Circle2IntrType::TwoIntersects:
							if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intrResult.point1)) {
								output.emplace_back(sIndex, intrResult.point1);
							}
							if (validArcSegIntersect(arc.center, v1.pos(), v2.pos(), v1.bulge(), intrResult.point2)) {
								output.emplace_back(sIndex, intrResult.point2);
							}
							break;
						case Circle2Circle2IntrType::Coincident:
							break;
						}
					}
				}
			}

			/// Function to test if a point is a valid distance from the original polyline.
			template <typename Real, std::size_t N>
			bool pointValidForOffset(Polyline<Real> const &pline, Real offset,
				StaticSpatialIndex<Real, N> const &spatialIndex,
				Vector2<Real> const &point, Real offsetTol = Real(1e-3)) {
				const Real absOffset = std::abs(offset) - offsetTol;
				const Real minDist = absOffset * absOffset;

				bool pointValid = true;

				auto visitor = [&](std::size_t i) {
					std::size_t j = utils::nextWrappingIndex(i, pline.vertexes());
					auto closestPoint = closestPointOnSeg(pline[i], pline[j], point);
					Real dist = distSquared(closestPoint, point);
					pointValid = dist > minDist;
					return pointValid;
				};

				spatialIndex.visitQuery(point.x() - absOffset, point.y() - absOffset, point.x() + absOffset,
					point.y() + absOffset, visitor);
				return pointValid;
			}

		} // namespace detail

		/// Creates the raw offset polyline.
		template <typename Real>
		Polyline<Real> createRawOffsetPline(Polyline<Real> const &pline, Real offset) {

			Polyline<Real> result;
			if (pline.size() < 2) {
				return result;
			}

			result.vertexes().reserve(pline.size());
			result.isClosed() = pline.isClosed();

			std::vector<PlineOffsetSegment<Real>> rawOffsets = createUntrimmedOffsetSegments(pline, offset);
			if (rawOffsets.size() == 0) {
				return result;
			}

			const bool connectionArcsAreCCW = offset < Real(0);

			auto joinResultVisitor = [connectionArcsAreCCW](PlineOffsetSegment<Real> const &s1,
				PlineOffsetSegment<Real> const &s2,
				Polyline<Real> &result) {
				const bool s1IsLine = s1.v1.bulgeIsZero();
				const bool s2IsLine = s2.v1.bulgeIsZero();
				if (s1IsLine && s2IsLine) {
					detail::lineToLineJoin(s1, s2, connectionArcsAreCCW, result);
				}
				else if (s1IsLine) {
					detail::lineToArcJoin(s1, s2, connectionArcsAreCCW, result);
				}
				else if (s2IsLine) {
					detail::arcToLineJoin(s1, s2, connectionArcsAreCCW, result);
				}
				else {
					detail::arcToArcJoin(s1, s2, connectionArcsAreCCW, result);
				}
			};

			result.addVertex(rawOffsets[0].v1);

			for (std::size_t i = 1; i < rawOffsets.size(); ++i) {
				const auto &seg1 = rawOffsets[i - 1];
				const auto &seg2 = rawOffsets[i];
				joinResultVisitor(seg1, seg2, result);
			}

			if (pline.isClosed() && result.size() > 1) {
				// joining segments at vertex indexes (n, 0) and (0, 1)
				const auto &s1 = rawOffsets.back();
				const auto &s2 = rawOffsets[0];

				// temp polyline to capture results of joining (to avoid mutating result)
				Polyline<Real> closingPartResult;
				closingPartResult.addVertex(result.lastVertex());
				joinResultVisitor(s1, s2, closingPartResult);

				// update last vertexes
				result.lastVertex() = closingPartResult[0];
				for (std::size_t i = 1; i < closingPartResult.size(); ++i) {
					result.addVertex(closingPartResult[i]);
				}
				result.vertexes().pop_back();

				// update first vertex
				const Vector2<Real> &updatedFirstPos = closingPartResult.lastVertex().pos();
				if (result[0].bulgeIsZero()) {
					// just update position
					result[0].pos() = updatedFirstPos;
				}
				else {
					// update position and bulge
					const auto arc = arcRadiusAndCenter(result[0], result[1]);
					const Real a1 = angle(arc.center, updatedFirstPos);
					const Real a2 = angle(arc.center, result[1].pos());
					const Real updatedTheta = utils::deltaAngle(a1, a2);
					if ((updatedTheta < Real(0) && result[0].bulge() > Real(0)) ||
						(updatedTheta > Real(0) && result[0].bulge() < Real(0))) {
						// first vertex not valid, just update its position to be removed later
						result[0].pos() = updatedFirstPos;
					}
					else {
						// update position and bulge
						result[0].pos() = updatedFirstPos;
						result[0].bulge() = std::tan(updatedTheta / Real(4));
					}
				}

				// must do final singularity prune between first and second vertex after joining curves (n, 0)
				// and (0, 1)
				if (result.size() > 1) {
					if (fuzzyEqual(result[0].pos(), result[1].pos(), utils::realPrecision<Real>())) {
						result.vertexes().erase(result.vertexes().begin());
					}
				}
			}
			else {
				detail::addOrReplaceIfSamePos(result, rawOffsets.back().v2);
			}

			return result;
		}

		/// Enum to represent the type of polyline intersect. Simple is a basic intersect, tangent is a
		/// tangent intersect (between an arc and line or two arcs) and coincident represents when two
		/// segments overlap.
		enum class PlineIntersectType { Simple, Tangent, Coincident };

		/// Represents a polyline intersect.
		template <typename Real> struct PlineIntersect {
			// index of start vertex of first segment
			std::size_t sIndex1;
			// index of start vertex of second segment
			std::size_t sIndex2;
			// intersect position
			Vector2<Real> pos;
			// type of intersect
			PlineIntersectType intrType;
			PlineIntersect(std::size_t si1, std::size_t si2, Vector2<Real> p, PlineIntersectType iType)
				: sIndex1(si1), sIndex2(si2), pos(p), intrType(iType) {}
		};

		/// Finds all local self intersects of the polyline, local self intersects are defined as between
		/// two polyline segments that share a vertex. NOTES:
		/// - Singularities (repeating vertexes) are returned as coincident intersects
		template <typename Real>
		void localSelfIntersects(Polyline<Real> const &pline, std::vector<PlineIntersect<Real>> &output) {
			if (pline.size() < 2) {
				return;
			}

			if (pline.size() == 2) {
				if (pline.isClosed()) {
					// check if overlaps on itself from vertex 1 to vertex 2
					if (utils::fuzzyEqual(pline[0].bulge(), -pline[1].bulge())) {
						output.emplace_back(0, 1, pline[1].pos(), PlineIntersectType::Coincident);
						output.emplace_back(1, 0, pline[0].pos(), PlineIntersectType::Coincident);
					}
				}
				return;
			}

			auto testAndAddIntersect = [&](std::size_t i, std::size_t j, std::size_t k) {
				const PlineVertex<Real> &v1 = pline[i];
				const PlineVertex<Real> &v2 = pline[j];
				const PlineVertex<Real> &v3 = pline[k];
				// testing intersection between v1->v2 and v2->v3 segments

				if (fuzzyEqual(v1.pos(), v2.pos(), utils::realPrecision<Real>())) {
					// singularity
					output.emplace_back(i, j, v1.pos(), PlineIntersectType::Coincident);
				}
				else {
					using namespace detail;
					IntrPlineSegsResult<Real> intrResult = intrPlineSegs(v1, v2, v2, v3);
					switch (intrResult.intrType) {
					case PlineSegIntrType::NoIntersect:
						break;
					case PlineSegIntrType::TangentIntersect:
						if (!fuzzyEqual(intrResult.point1, v2.pos(), utils::realPrecision<Real>())) {
							output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Tangent);
						}
						break;
					case PlineSegIntrType::OneIntersect:
						if (!fuzzyEqual(intrResult.point1, v2.pos(), utils::realPrecision<Real>())) {
							output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Simple);
						}
						break;
					case PlineSegIntrType::TwoIntersects:
						if (!fuzzyEqual(intrResult.point1, v2.pos(), utils::realPrecision<Real>())) {
							output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Simple);
						}
						if (!fuzzyEqual(intrResult.point2, v2.pos(), utils::realPrecision<Real>())) {
							output.emplace_back(i, j, intrResult.point2, PlineIntersectType::Simple);
						}
						break;
					case PlineSegIntrType::SegmentOverlap:
					case PlineSegIntrType::ArcOverlap:
						output.emplace_back(i, j, intrResult.point1, PlineIntersectType::Coincident);
						break;
					}
				}
			};

			for (std::size_t i = 2; i < pline.size(); ++i) {
				testAndAddIntersect(i - 2, i - 1, i);
			}

			if (pline.isClosed()) {
				// we tested for intersect between segments at indexes 0->1, 1->2 and everything up to and
				// including (count-3)->(count-2), (count-2)->(count-1), polyline is closed so now test
				// [(count-2)->(count-1), (count-1)->0] and [(count-1)->0, 0->1]
				testAndAddIntersect(pline.size() - 2, pline.size() - 1, 0);
				testAndAddIntersect(pline.size() - 1, 0, 1);
			}
		}

		/// Finds all global self intersects of the polyline, global self intersects are defined as all
		/// intersects between polyline segments that DO NOT share a vertex (use the localSelfIntersects
		/// function to find those). A spatial index is used to minimize the intersect comparisons required,
		/// the spatial index should hold bounding boxes for all of the polyline's segments.
		/// NOTES:
		/// - We never include intersects at a segment's start point, the matching intersect from the
		/// previous segment's end point is included (no sense in including both)
		template <typename Real, std::size_t N>
		void globalSelfIntersects(Polyline<Real> const &pline, std::vector<PlineIntersect<Real>> &output,
			StaticSpatialIndex<Real, N> const &spatialIndex) {
			if (pline.size() < 3) {
				return;
			}
			using namespace detail;

			std::unordered_set<std::pair<std::size_t, std::size_t>, IndexPairHash> visitedSegmentPairs;

			auto visitor = [&](std::size_t i, std::size_t j) {
				const PlineVertex<Real> &v1 = pline[i];
				const PlineVertex<Real> &v2 = pline[j];
				AABB<Real> envelope = createFastApproxBoundingBox(v1, v2);
				envelope.expand(utils::realThreshold<Real>());
				auto indexVisitor = [&](std::size_t hitIndexStart) {
					std::size_t hitIndexEnd = utils::nextWrappingIndex(hitIndexStart, pline);
					// skip/filter already visited intersects
					// skip local segments
					if (i == hitIndexStart || i == hitIndexEnd || j == hitIndexStart || j == hitIndexEnd) {
						return true;
					}
					// skip reversed segment order (would end up comparing the same segments)
					if (visitedSegmentPairs.find(std::make_pair(hitIndexStart, i)) != visitedSegmentPairs.end()) {
						return true;
					}

					// add the segment pair we're visiting now
					visitedSegmentPairs.emplace(i, hitIndexStart);

					const PlineVertex<Real> &u1 = pline[hitIndexStart];
					const PlineVertex<Real> &u2 = pline[hitIndexEnd];

					auto intrAtStartPt = [&](Vector2<Real> const &intr) {
						return fuzzyEqual(v1.pos(), intr) || fuzzyEqual(u1.pos(), intr);
					};

					IntrPlineSegsResult<Real> intrResult = intrPlineSegs(v1, v2, u1, u2);
					switch (intrResult.intrType) {
					case PlineSegIntrType::NoIntersect:
						break;
					case PlineSegIntrType::TangentIntersect:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Tangent);
						}
						break;
					case PlineSegIntrType::OneIntersect:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Simple);
						}
						break;
					case PlineSegIntrType::TwoIntersects:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Simple);
						}
						if (!intrAtStartPt(intrResult.point2)) {
							output.emplace_back(i, hitIndexStart, intrResult.point2, PlineIntersectType::Simple);
						}
						break;
					case PlineSegIntrType::SegmentOverlap:
					case PlineSegIntrType::ArcOverlap:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i, hitIndexStart, intrResult.point1, PlineIntersectType::Coincident);
						}
						if (!intrAtStartPt(intrResult.point2)) {
							output.emplace_back(i, hitIndexStart, intrResult.point2, PlineIntersectType::Coincident);
						}
						break;
					}

					// visit the entire query
					return true;
				};

				spatialIndex.visitQuery(envelope.xMin, envelope.yMin, envelope.xMax, envelope.yMax,
					indexVisitor);

				// visit all pline indexes
				return true;
			};

			iterateSegIndices(pline, visitor);
		}

		/// Finds all self intersects of the polyline (equivalent to calling localSelfIntersects and
		/// globalSelfIntersects).
		template <typename Real, std::size_t N>
		void allSelfIntersects(Polyline<Real> const &pline, std::vector<PlineIntersect<Real>> &output,
			StaticSpatialIndex<Real, N> const &spatialIndex) {
			localSelfIntersects(pline, output);
			globalSelfIntersects(pline, output, spatialIndex);
		}

		/// Represents an open polyline slice of the raw offset polyline.
		template <typename Real> struct OpenPolylineSlice {
			std::size_t intrStartIndex;
			Polyline<Real> pline;
			OpenPolylineSlice() = default;
			OpenPolylineSlice(std::size_t sIndex, Polyline<Real> slice)
				: intrStartIndex(sIndex), pline(std::move(slice)) {}
		};

		/// Finds all intersects between pline1 and pline2.
		template <typename Real, std::size_t N>
		void intersects(Polyline<Real> const &pline1, Polyline<Real> const &pline2,
			StaticSpatialIndex<Real, N> const &pline1SpatialIndex,
			std::vector<PlineIntersect<Real>> &output) {
			std::vector<std::size_t> queryResults;
			auto pline2SegVisitor = [&](std::size_t i2, std::size_t j2) {
				PlineVertex<Real> const &p2v1 = pline2[i2];
				PlineVertex<Real> const &p2v2 = pline2[j2];

				queryResults.clear();

				AABB<Real> bb = createFastApproxBoundingBox(p2v1, p2v2);
				pline1SpatialIndex.query(bb.xMin, bb.yMin, bb.xMax, bb.yMax, queryResults);

				using namespace detail;

				for (std::size_t i1 : queryResults) {
					std::size_t j1 = utils::nextWrappingIndex(i1, pline1);
					PlineVertex<Real> const &p1v1 = pline1[i1];
					PlineVertex<Real> const &p1v2 = pline1[j1];

					auto intrAtStartPt = [&](Vector2<Real> const &intr) {
						return fuzzyEqual(p1v1.pos(), intr) || fuzzyEqual(p2v1.pos(), intr);
					};

					auto intrResult = intrPlineSegs(p1v1, p1v2, p2v1, p2v2);
					switch (intrResult.intrType) {
					case PlineSegIntrType::NoIntersect:
						break;
					case PlineSegIntrType::TangentIntersect:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i1, i2, intrResult.point1, PlineIntersectType::Tangent);
						}
						break;
					case PlineSegIntrType::OneIntersect:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i1, i2, intrResult.point1, PlineIntersectType::Simple);
						}
						break;
					case PlineSegIntrType::TwoIntersects:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i1, i2, intrResult.point1, PlineIntersectType::Simple);
						}
						if (!intrAtStartPt(intrResult.point2)) {
							output.emplace_back(i1, i2, intrResult.point2, PlineIntersectType::Simple);
						}
						break;
					case PlineSegIntrType::SegmentOverlap:
					case PlineSegIntrType::ArcOverlap:
						if (!intrAtStartPt(intrResult.point1)) {
							output.emplace_back(i1, i2, intrResult.point1, PlineIntersectType::Coincident);
						}
						if (!intrAtStartPt(intrResult.point2)) {
							output.emplace_back(i1, i2, intrResult.point2, PlineIntersectType::Coincident);
						}
						break;
					}
				}

				// visit all indexes
				return true;
			};

			iterateSegIndices(pline2, pline2SegVisitor);
		}

		/// Slices a raw offset polyline at all of its self intersects.
		template <typename Real>
		std::vector<OpenPolylineSlice<Real>> sliceAtIntersects(Polyline<Real> const &originalPline,
			Polyline<Real> const &rawOffsetPline,
			Real offset) {
			assert(originalPline.isClosed() && "use dual slice at intersects for open polylines");

			std::vector<OpenPolylineSlice<Real>> result;
			if (rawOffsetPline.size() < 2) {
				return result;
			}

			StaticSpatialIndex<Real> origPlineSpatialIndex = createApproxSpatialIndex(originalPline);
			StaticSpatialIndex<Real> rawOffsetPlineSpatialIndex = createApproxSpatialIndex(rawOffsetPline);

			std::vector<PlineIntersect<Real>> selfIntersects;
			allSelfIntersects(rawOffsetPline, selfIntersects, rawOffsetPlineSpatialIndex);

			if (selfIntersects.size() == 0) {
				// no intersects, test that all vertexes are valid distance from original polyline
				bool isValid = std::all_of(rawOffsetPline.vertexes().begin(), rawOffsetPline.vertexes().end(),
					[&](PlineVertex<Real> const &v) {
					return detail::pointValidForOffset(originalPline, offset,
						origPlineSpatialIndex, v.pos());
				});

				if (!isValid) {
					return result;
				}

				// copy and convert raw offset into open polyline
				result.emplace_back(std::numeric_limits<std::size_t>::max(), rawOffsetPline);
				result.back().pline.isClosed() = false;
				result.back().pline.addVertex(rawOffsetPline[0]);
				result.back().pline.lastVertex().bulge() = Real(0);
				return result;
			}

			std::unordered_map<std::size_t, std::vector<Vector2<Real>>> intersectsLookup;
			intersectsLookup.reserve(2 * selfIntersects.size());

			for (PlineIntersect<Real> const &si : selfIntersects) {
				intersectsLookup[si.sIndex1].push_back(si.pos);
				intersectsLookup[si.sIndex2].push_back(si.pos);
			}

			// sort intersects by distance from start vertex
			for (auto &kvp : intersectsLookup) {
				Vector2<Real> startPos = rawOffsetPline[kvp.first].pos();
				auto cmp = [&](Vector2<Real> const &si1, Vector2<Real> const &si2) {
					return distSquared(si1, startPos) < distSquared(si2, startPos);
				};
				std::sort(kvp.second.begin(), kvp.second.end(), cmp);
			}

			auto intersectsOrigPline = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
				AABB<Real> approxBB = createFastApproxBoundingBox(v1, v2);
				bool hasIntersect = false;
				auto visitor = [&](std::size_t i) {
					using namespace detail;
					std::size_t j = utils::nextWrappingIndex(i, originalPline);
					IntrPlineSegsResult<Real> intrResult =
						intrPlineSegs(v1, v2, originalPline[i], originalPline[j]);
					hasIntersect = intrResult.intrType != PlineSegIntrType::NoIntersect;
					return !hasIntersect;
				};

				origPlineSpatialIndex.visitQuery(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax,
					visitor);

				return hasIntersect;
			};

			for (auto const &kvp : intersectsLookup) {
				// start index for the slice we're about to build
				std::size_t sIndex = kvp.first;
				// self intersect list for this start index
				std::vector<Vector2<Real>> const &siList = kvp.second;

				const auto &startVertex = rawOffsetPline[sIndex];
				std::size_t nextIndex = utils::nextWrappingIndex(sIndex, rawOffsetPline);
				const auto &endVertex = rawOffsetPline[nextIndex];

				if (siList.size() != 1) {
					// build all the segments between the N intersects in siList (N > 1)
					SplitResult<Real> firstSplit = splitAtPoint(startVertex, endVertex, siList[0]);
					auto prevVertex = firstSplit.splitVertex;
					for (std::size_t i = 1; i < siList.size(); ++i) {
						SplitResult<Real> split = splitAtPoint(prevVertex, endVertex, siList[i]);
						// update prevVertex for next loop iteration
						prevVertex = split.splitVertex;
						// skip if they're ontop of each other
						if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
							utils::realPrecision<Real>())) {
							continue;
						}

						// test start point
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							split.updatedStart.pos())) {
							continue;
						}

						// test end point
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							split.splitVertex.pos())) {
							continue;
						}

						// test mid point
						auto midpoint = detail::segMidpoint(split.updatedStart, split.splitVertex);
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, midpoint)) {
							continue;
						}

						// test intersection with original polyline
						if (intersectsOrigPline(split.updatedStart, split.splitVertex)) {
							continue;
						}

						result.emplace_back();
						result.back().intrStartIndex = sIndex;
						result.back().pline.addVertex(split.updatedStart);
						result.back().pline.addVertex(split.splitVertex);
					}
				}

				// build the segment between the last intersect in siList and the next intersect found

				// check that the first point is valid
				if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, siList.back())) {
					continue;
				}

				SplitResult<Real> split = splitAtPoint(startVertex, endVertex, siList.back());
				Polyline<Real> currPline;
				currPline.addVertex(split.splitVertex);

				std::size_t index = nextIndex;
				bool isValidPline = true;
				while (true) {
					// check that vertex point is valid
					if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
						rawOffsetPline[index].pos())) {
						isValidPline = false;
						break;
					}

					// check that the segment does not intersect original polyline
					if (intersectsOrigPline(currPline.lastVertex(), rawOffsetPline[index])) {
						isValidPline = false;
						break;
					}

					// add vertex
					detail::addOrReplaceIfSamePos(currPline, rawOffsetPline[index]);

					// check if segment that starts at vertex we just added has an intersect
					auto nextIntr = intersectsLookup.find(index);
					if (nextIntr != intersectsLookup.end()) {
						// there is an intersect, slice is done, check if final segment is valid

						// check intersect pos is valid (which will also be end vertex position)
						Vector2<Real> const &intersectPos = nextIntr->second[0];
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							intersectPos)) {
							isValidPline = false;
							break;
						}

						std::size_t nextIndex = utils::nextWrappingIndex(index, rawOffsetPline);
						SplitResult<Real> split =
							splitAtPoint(currPline.lastVertex(), rawOffsetPline[nextIndex], intersectPos);

						PlineVertex<Real> endVertex = PlineVertex<Real>(intersectPos, Real(0));
						// check mid point is valid
						Vector2<Real> mp = detail::segMidpoint(split.updatedStart, endVertex);
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, mp)) {
							isValidPline = false;
							break;
						}

						// trim last added vertex and add final intersect position
						currPline.lastVertex() = split.updatedStart;
						detail::addOrReplaceIfSamePos(currPline, endVertex);

						break;
					}
					// else there is not an intersect, increment index and continue
					index = utils::nextWrappingIndex(index, rawOffsetPline);
				}

				if (isValidPline && currPline.size() > 1) {
					result.emplace_back(sIndex, std::move(currPline));
				}
			}

			return result;
		}

		/// Slices a raw offset polyline at all of its self intersects and intersects with its dual.
		template <typename Real>
		std::vector<OpenPolylineSlice<Real>>
			dualSliceAtIntersects(Polyline<Real> const &originalPline, Polyline<Real> const &rawOffsetPline,
				Polyline<Real> const &dualRawOffsetPline, Real offset) {
			std::vector<OpenPolylineSlice<Real>> result;
			if (rawOffsetPline.size() < 2) {
				return result;
			}

			StaticSpatialIndex<Real> origPlineSpatialIndex = createApproxSpatialIndex(originalPline);
			StaticSpatialIndex<Real> rawOffsetPlineSpatialIndex = createApproxSpatialIndex(rawOffsetPline);

			std::vector<PlineIntersect<Real>> selfIntersects;
			allSelfIntersects(rawOffsetPline, selfIntersects, rawOffsetPlineSpatialIndex);

			std::vector<PlineIntersect<Real>> dualIntersects;
			intersects(rawOffsetPline, dualRawOffsetPline, rawOffsetPlineSpatialIndex, dualIntersects);

			std::unordered_map<std::size_t, std::vector<Vector2<Real>>> intersectsLookup;

			if (!originalPline.isClosed()) {
				// find intersects between circles generated at original open polyline end points and raw offset
				// polyline
				std::vector<std::pair<std::size_t, Vector2<Real>>> intersects;
				detail::offsetCircleIntersectsWithPline(rawOffsetPline, offset, originalPline[0].pos(),
					rawOffsetPlineSpatialIndex, intersects);
				detail::offsetCircleIntersectsWithPline(rawOffsetPline, offset,
					originalPline.lastVertex().pos(),
					rawOffsetPlineSpatialIndex, intersects);
				intersectsLookup.reserve(2 * selfIntersects.size() + intersects.size());
				for (auto const &pair : intersects) {
					intersectsLookup[pair.first].push_back(pair.second);
				}
			}
			else {
				intersectsLookup.reserve(2 * selfIntersects.size());
			}

			for (PlineIntersect<Real> const &si : selfIntersects) {
				intersectsLookup[si.sIndex1].push_back(si.pos);
				intersectsLookup[si.sIndex2].push_back(si.pos);
			}

			for (PlineIntersect<Real> const &intr : dualIntersects) {
				intersectsLookup[intr.sIndex1].push_back(intr.pos);
			}

			if (intersectsLookup.size() == 0) {
				// no intersects, test that all vertexes are valid distance from original polyline
				bool isValid = std::all_of(rawOffsetPline.vertexes().begin(), rawOffsetPline.vertexes().end(),
					[&](PlineVertex<Real> const &v) {
					return detail::pointValidForOffset(originalPline, offset,
						origPlineSpatialIndex, v.pos());
				});

				if (!isValid) {
					return result;
				}

				// copy and convert raw offset into open polyline
				result.emplace_back(std::numeric_limits<std::size_t>::max(), rawOffsetPline);
				result.back().pline.isClosed() = false;
				if (originalPline.isClosed()) {
					result.back().pline.addVertex(rawOffsetPline[0]);
					result.back().pline.lastVertex().bulge() = Real(0);
				}
				return result;
			}

			// sort intersects by distance from start vertex
			for (auto &kvp : intersectsLookup) {
				Vector2<Real> startPos = rawOffsetPline[kvp.first].pos();
				auto cmp = [&](Vector2<Real> const &si1, Vector2<Real> const &si2) {
					return distSquared(si1, startPos) < distSquared(si2, startPos);
				};
				std::sort(kvp.second.begin(), kvp.second.end(), cmp);
			}

			auto intersectsOrigPline = [&](PlineVertex<Real> const &v1, PlineVertex<Real> const &v2) {
				AABB<Real> approxBB = createFastApproxBoundingBox(v1, v2);
				bool intersects = false;
				auto visitor = [&](std::size_t i) {
					using namespace detail;
					std::size_t j = utils::nextWrappingIndex(i, originalPline);
					IntrPlineSegsResult<Real> intrResult =
						intrPlineSegs(v1, v2, originalPline[i], originalPline[j]);
					intersects = intrResult.intrType != PlineSegIntrType::NoIntersect;
					return !intersects;
				};

				origPlineSpatialIndex.visitQuery(approxBB.xMin, approxBB.yMin, approxBB.xMax, approxBB.yMax,
					visitor);

				return intersects;
			};

			if (!originalPline.isClosed()) {
				// build first open polyline that ends at the first intersect since we will not wrap back to
				// capture it as in the case of a closed polyline
				Polyline<Real> firstPline;
				std::size_t index = 0;
				while (true) {
					auto iter = intersectsLookup.find(index);
					if (iter == intersectsLookup.end()) {
						// no intersect found, test segment will be valid before adding the vertex
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							rawOffsetPline[index].pos())) {
							break;
						}

						// index check (only test segment if we're not adding the first vertex)
						if (index != 0 && intersectsOrigPline(firstPline.lastVertex(), rawOffsetPline[index])) {
							break;
						}

						detail::addOrReplaceIfSamePos(firstPline, rawOffsetPline[index]);
					}
					else {
						// intersect found, test segment will be valid before finishing first open polyline
						Vector2<Real> const &intersectPos = iter->second[0];
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							intersectPos)) {
							break;
						}

						SplitResult<Real> split =
							splitAtPoint(rawOffsetPline[index], rawOffsetPline[index + 1], intersectPos);

						PlineVertex<Real> endVertex = PlineVertex<Real>(intersectPos, Real(0));
						auto midpoint = detail::segMidpoint(split.updatedStart, endVertex);
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, midpoint)) {
							break;
						}

						if (intersectsOrigPline(split.updatedStart, endVertex)) {
							break;
						}

						detail::addOrReplaceIfSamePos(firstPline, split.updatedStart);
						detail::addOrReplaceIfSamePos(firstPline, endVertex);
						result.emplace_back(0, std::move(firstPline));
						break;
					}

					index += 1;
				}
			}

			for (auto const &kvp : intersectsLookup) {
				// start index for the slice we're about to build
				std::size_t sIndex = kvp.first;
				// self intersect list for this start index
				std::vector<Vector2<Real>> const &siList = kvp.second;

				const auto &startVertex = rawOffsetPline[sIndex];
				std::size_t nextIndex = utils::nextWrappingIndex(sIndex, rawOffsetPline);
				const auto &endVertex = rawOffsetPline[nextIndex];

				if (siList.size() != 1) {
					// build all the segments between the N intersects in siList (N > 1)
					SplitResult<Real> firstSplit = splitAtPoint(startVertex, endVertex, siList[0]);
					auto prevVertex = firstSplit.splitVertex;
					for (std::size_t i = 1; i < siList.size(); ++i) {
						SplitResult<Real> split = splitAtPoint(prevVertex, endVertex, siList[i]);
						// update prevVertex for next loop iteration
						prevVertex = split.splitVertex;
						// skip if they're ontop of each other
						if (fuzzyEqual(split.updatedStart.pos(), split.splitVertex.pos(),
							utils::realPrecision<Real>())) {
							continue;
						}

						// test start point
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							split.updatedStart.pos())) {
							continue;
						}

						// test end point
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							split.splitVertex.pos())) {
							continue;
						}

						// test mid point
						auto midpoint = detail::segMidpoint(split.updatedStart, split.splitVertex);
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, midpoint)) {
							continue;
						}

						// test intersection with original polyline
						if (intersectsOrigPline(split.updatedStart, split.splitVertex)) {
							continue;
						}

						result.emplace_back();
						result.back().intrStartIndex = sIndex;
						result.back().pline.addVertex(split.updatedStart);
						result.back().pline.addVertex(split.splitVertex);
					}
				}

				// build the segment between the last intersect in siList and the next intersect found

				// check that the first point is valid
				if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, siList.back())) {
					continue;
				}

				SplitResult<Real> split = splitAtPoint(startVertex, endVertex, siList.back());
				Polyline<Real> currPline;
				currPline.addVertex(split.splitVertex);

				std::size_t index = nextIndex;
				bool isValidPline = true;
				while (true) {
					// check that vertex point is valid
					if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
						rawOffsetPline[index].pos())) {
						isValidPline = false;
						break;
					}

					// check that the segment does not intersect original polyline
					if (intersectsOrigPline(currPline.lastVertex(), rawOffsetPline[index])) {
						isValidPline = false;
						break;
					}

					// add vertex
					detail::addOrReplaceIfSamePos(currPline, rawOffsetPline[index]);

					// check if segment that starts at vertex we just added has an intersect
					auto nextIntr = intersectsLookup.find(index);
					if (nextIntr != intersectsLookup.end()) {
						// there is an intersect, slice is done, check if final segment is valid

						// check intersect pos is valid (which will also be end vertex position)
						Vector2<Real> const &intersectPos = nextIntr->second[0];
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex,
							intersectPos)) {
							isValidPline = false;
							break;
						}

						std::size_t nextIndex = utils::nextWrappingIndex(index, rawOffsetPline);
						SplitResult<Real> split =
							splitAtPoint(currPline.lastVertex(), rawOffsetPline[nextIndex], intersectPos);

						PlineVertex<Real> endVertex = PlineVertex<Real>(intersectPos, Real(0));
						// check mid point is valid
						Vector2<Real> mp = detail::segMidpoint(split.updatedStart, endVertex);
						if (!detail::pointValidForOffset(originalPline, offset, origPlineSpatialIndex, mp)) {
							isValidPline = false;
							break;
						}

						// trim last added vertex and add final intersect position
						currPline.lastVertex() = split.updatedStart;
						detail::addOrReplaceIfSamePos(currPline, endVertex);

						break;
					}
					// else there is not an intersect, increment index and continue
					if (index == rawOffsetPline.size() - 1) {
						if (originalPline.isClosed()) {
							// wrap index
							index = 0;
						}
						else {
							// open polyline, we're done
							break;
						}
					}
					else {
						index += 1;
					}
				}

				if (isValidPline && currPline.size() > 1) {
					result.emplace_back(sIndex, std::move(currPline));
				}
			}

			return result;
		}

		/// Stitches raw offset polyline slices together, discarding any that are not valid.
		template <typename Real>
		std::vector<Polyline<Real>>
			stitchSlicesTogether(std::vector<OpenPolylineSlice<Real>> const &slices, bool closedPolyline,
				std::size_t origMaxIndex, Real joinThreshold = utils::realPrecision<Real>()) {
			std::vector<Polyline<Real>> result;
			if (slices.size() == 0) {
				return result;
			}

			if (slices.size() == 1) {
				result.emplace_back(slices[0].pline);
				if (closedPolyline) {
					result.back().isClosed() = true;
					result.back().vertexes().pop_back();
				}

				return result;
			}

			// load spatial index with all start points
			StaticSpatialIndex<Real> spatialIndex(slices.size());

			for (const auto &slice : slices) {
				auto const &point = slice.pline[0].pos();
				spatialIndex.add(point.x() - joinThreshold, point.y() - joinThreshold,
					point.x() + joinThreshold, point.y() + joinThreshold);
			}

			spatialIndex.finish();

			std::vector<bool> visitedIndexes(slices.size(), false);
			std::vector<std::size_t> queryResults;
			for (std::size_t i = 0; i < slices.size(); ++i) {
				if (visitedIndexes[i]) {
					continue;
				}

				visitedIndexes[i] = true;

				Polyline<Real> currPline;
				currPline.isClosed() = closedPolyline;
				std::size_t currIndex = i;
				auto const &initialStartPoint = slices[i].pline[0].pos();

				while (true) {
					const std::size_t currLoopStartIndex = slices[currIndex].intrStartIndex;
					auto const &currSlice = slices[currIndex].pline;
					auto const &currEndPoint = slices[currIndex].pline.lastVertex().pos();
					currPline.vertexes().insert(currPline.vertexes().end(), currSlice.vertexes().begin(),
						currSlice.vertexes().end());
					queryResults.clear();
					spatialIndex.query(currEndPoint.x() - joinThreshold, currEndPoint.y() - joinThreshold,
						currEndPoint.x() + joinThreshold, currEndPoint.y() + joinThreshold,
						queryResults);

					queryResults.erase(std::remove_if(queryResults.begin(), queryResults.end(),
						[&](std::size_t index) { return visitedIndexes[index]; }),
						queryResults.end());

					auto indexDistAndEqualInitial = [&](std::size_t index) {
						auto const &slice = slices[index];
						std::size_t indexDist;
						if (currLoopStartIndex <= slice.intrStartIndex) {
							indexDist = slice.intrStartIndex - currLoopStartIndex;
						}
						else {
							// forward wrapping distance (distance to end + distance to index)
							indexDist = origMaxIndex - currLoopStartIndex + slice.intrStartIndex;
						}

						bool equalToInitial = fuzzyEqual(slice.pline.lastVertex().pos(), initialStartPoint,
							utils::realPrecision<Real>());

						return std::make_pair(indexDist, equalToInitial);
					};

					std::sort(queryResults.begin(), queryResults.end(),
						[&](std::size_t index1, std::size_t index2) {
						auto distAndEqualInitial1 = indexDistAndEqualInitial(index1);
						auto distAndEqualInitial2 = indexDistAndEqualInitial(index2);
						if (distAndEqualInitial1.first == distAndEqualInitial2.first) {
							// index distances are equal, compare on position being equal to initial start
							// (testing index1 < index2, we want the longest closed loop possible)
							return distAndEqualInitial1.second < distAndEqualInitial2.second;
						}

						return distAndEqualInitial1.first < distAndEqualInitial2.first;
					});

					if (queryResults.size() == 0) {
						// we're done
						if (currPline.size() > 1) {
							if (closedPolyline && fuzzyEqual(currPline[0].pos(), currPline.lastVertex().pos(),
								utils::realPrecision<Real>())) {
								currPline.vertexes().pop_back();
							}
							result.emplace_back(std::move(currPline));
						}
						break;
					}

					// else continue stitching
					visitedIndexes[queryResults[0]] = true;
					currPline.vertexes().pop_back();
					currIndex = queryResults[0];
				}
			}

			return result;
		}

		/// Creates the paralell offset polylines to the polyline given.
		template <typename Real>
		std::vector<Polyline<Real>> parallelOffset(Polyline<Real> const &pline, Real offset,
			bool hasSelfIntersects = false) {
			if (pline.size() < 2) {
				return std::vector<Polyline<Real>>();
			}
			auto rawOffset = createRawOffsetPline(pline, offset);
			if (pline.isClosed() && !hasSelfIntersects) {
				auto slices = sliceAtIntersects(pline, rawOffset, offset);
				return stitchSlicesTogether(slices, pline.isClosed(), rawOffset.size() - 1);
			}

			// not closed polyline or has self intersects, must apply dual clipping
			auto dualRawOffset = createRawOffsetPline(pline, -offset);
			auto slices = dualSliceAtIntersects(pline, rawOffset, dualRawOffset, offset);
			return stitchSlicesTogether(slices, pline.isClosed(), rawOffset.size() - 1);
		}

	}


	namespace offset {
		using namespace cavc;
		using namespace cavc::AlgoIntel;
		/*
		 * @brief 2d曲线的偏置,TVert重载了[]运算符
		 * @detail 默认沿着法向进行偏置,注意bulge是指曲线的凸度，当为多段线时候，均为0
		 * @param[in] pts 曲线原始采样点
		 * @param[in] offset 偏置距离
		 * @param[in] hasSelfIntersects 是否做自交处理
		 * @param[in] isClosed 原始曲线是否是闭合曲线
		 * @param[out] offsetpts 偏置后的曲线
		 * @return 偏置是否成功
		 * @date 2020.2.18
		 * @author huting,shangTao,lyc
		 */
		template<class TVert, class ftype, int dim = 2>
		bool offsetCurve2D(
			const std::vector<TVert>& pts, const std::vector<ftype>& bulge,
			ftype offset, bool hasSelfIntersects, bool isClosed,
			std::vector<std::vector<TVert>>& offsetpts)
		{
			using namespace cavc;
			static_assert(dim == 2, "only for 2d curve!");
			if (pts.size() != bulge.size())
				return false;
			int npt = static_cast<int>(pts.size());
			Polyline<ftype> input;
			for (int i = 0; i < npt; ++i) {
				input.addVertex(pts[i][0], pts[i][1], bulge[i]);
			}
			input.isClosed() = isClosed;
			std::vector<Polyline<ftype>> results = parallelOffset(input, offset, hasSelfIntersects);
			offsetpts.clear();
			offsetpts.resize(results.size());
			int nresult = static_cast<int>(results.size());
			for (int i = 0; i < nresult; ++i) {
				auto& p = results[i];
				offsetpts[i].reserve(p.size());
				for (int j = 0; j < p.size(); ++j) {
					offsetpts[i].push_back({ p[j].x(),p[j].y() });
				}
			}
			return true;
		}

		/*
		 * @brief 2d曲线的偏置,TVert重载了[]运算符
		 * @detail 默认沿着法向进行偏置,注意凸度均为0
		 * @param[in] pts 曲线原始采样点
		 * @param[in] offset 偏置距离
		 * @param[in] hasSelfIntersects 是否做自交处理
		 * @param[in] isClosed 原始曲线是否是闭合曲线
		 * @param[out] offsetpts 偏置后的曲线
		 * @return 偏置是否成功
		 * @date 2020.2.18
		 * @author huting,shangTao,lyc
		 */
		template<class TVert, class ftype, int dim = 2>
		bool offsetCurve2D(
			const std::vector<TVert>& pts,
			ftype offset, bool hasSelfIntersects, bool isClosed,
			std::vector<std::vector<TVert>>& offsetpts)
		{
			using namespace cavc;
			static_assert(dim == 2, "only for 2d curve!");
			int npt = static_cast<int>(pts.size());
			Polyline<ftype> input;
			for (int i = 0; i < npt; ++i) {
				input.addVertex(pts[i][0], pts[i][1], 0);
			}
			input.isClosed() = isClosed;
			std::vector<Polyline<ftype>> results = parallelOffset(input, offset, hasSelfIntersects);
			offsetpts.clear();
			offsetpts.resize(results.size());
			int nresult = static_cast<int>(results.size());
			for (int i = 0; i < nresult; ++i) {
				auto& p = results[i];
				offsetpts[i].reserve(p.size());
				for (int j = 0; j < p.size(); ++j) {
					offsetpts[i].push_back({ p[j].x(),p[j].y() });
				}
			}
			return true;
		}
		template<class TVert>
		TVert crossproduct(const TVert& p, const TVert& q) {
			TVert cp = { p[1] * q[2] - p[2] * q[1],
				p[2] * q[0] - p[0] * q[2],
				p[0] * q[1] - p[1] * q[0] };
			return cp;
		}

		/*
		 * @brief 偏置2_5D曲线(以及可以转换为2.5D的曲线),TVert重载了[]运算符
		 * @detail 在空间某些简单曲线可以投影到空间平面上，形成3d平面曲线的曲线称为2.5d曲线
		 * @param[in] pts 曲线原始采样点
		 * @param[in] offset 偏置距离
		 * @param[in] hasSelfIntersects 是否做自交处理
		 * @param[in] isClosed 原始曲线是否是闭合曲线
		 * @param[out] offsetpts 偏置后的曲线
		 * @return 偏置是否成功
		 * @date 2020.2.18
		 * @author huting,shangTao,lyc
		 */
		template<class TVert, class ftype, int dim = 3>
		bool offsetCurve2_5D(
			const std::vector<TVert>& pts, const TVert& norm,
			ftype offset, bool hasSelfIntersects, bool isClosed,
			std::vector<std::vector<TVert>>& offsetpts)
		{
			using namespace cavc;
			using namespace cavc::AlgoIntel;
			static_assert(dim == 3, "only for 3d points");
			if (pts.size() < 2)
				return false;
			int npt = static_cast<int>(pts.size());
			//参数化曲线
			TVert udir = SubT(pts[1], pts[0]);
			ftype ulen = static_cast<ftype>(1 /
				sqrt(udir[0] * udir[0] + udir[1] * udir[1] + udir[2] * udir[2]));
			udir = ndotproduct<TVert, ftype>(ulen, udir);
			TVert vdir = crossproduct(udir, norm);
			ftype vlen = static_cast<ftype>(1 /
				sqrt(vdir[0] * vdir[0] + vdir[1] * vdir[1] + vdir[2] * vdir[2]));
			auto wdir = norm;
			ftype clen = static_cast<ftype>(1 /
				sqrt(wdir[0] * wdir[0] + wdir[1] * wdir[1] + wdir[2] * wdir[2]));
			wdir = ndotproduct<TVert, ftype>(clen, wdir);
			vdir = ndotproduct<TVert, ftype>(vlen, vdir);
			Polyline<ftype> input;
			std::vector<ftype> wpts;
			for (int i = 0; i < npt; ++i) {
				TVert p = SubT(pts[i], pts[0]);
				ftype u = dotproduct<TVert, ftype>(p, udir);
				ftype v = dotproduct<TVert, ftype>(p, vdir);
				ftype w = dotproduct<TVert, ftype>(p, wdir);
				wpts.push_back(w);
				input.addVertex(u, v, 0);
			}
			input.isClosed() = isClosed;
			std::vector<Polyline<ftype>> results = parallelOffset(input, offset, hasSelfIntersects);
			offsetpts.clear();
			offsetpts.resize(results.size());
			int nres = static_cast<int>(results.size());
			for (int i = 0; i < nres; ++i) {
				offsetpts[i].reserve(results[i].size());
				for (int j = 0; j < results[i].size(); ++j) {
					ftype u = results[i][j].x();
					ftype v = results[i][j].y();
					TVert p = Add(ndotproduct(u, udir), ndotproduct(v, vdir));
					TVert q = Add(pts[0], p);
					offsetpts[i].push_back(q);
				}
			}
			return true;
		}
	}



}















#endif