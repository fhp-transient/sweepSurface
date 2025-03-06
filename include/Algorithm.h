#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <random>

#include "glm/glm.hpp"
#include "myutil.h"
#include "tinynurbs/core/curve.h"
#include "tinynurbs/core/surface.h"
#include "tinynurbs/util/array2.h"
#include "tinynurbs/util/util.h"
#include <eigen/Eigen/Dense>
#include <glm/ext/scalar_constants.hpp>
#include <map>
#include <optional>
#include <tinynurbs/tinynurbs.h>
#include <vector>

template<typename T>
unsigned int findSpan(unsigned int degree, const std::vector<T> &knots, T u)
{
    // index of last control point
    int n = static_cast<int>(knots.size()) - degree - 2;
    assert(n >= 0);

    // For u that is equal to last knot value
    if (equalTo(u, knots[n + 1]))
    {
        return n;
    }

    // For values of u that lies outside the domain
    if (equalTo(u, knots[n + 1]) ||
        u > (knots[n + 1] - std::numeric_limits<T>::epsilon()))
    {
        return n;
    }
    if (u < (knots[degree] + std::numeric_limits<T>::epsilon()))
    {
        return degree;
    }

    //// Binary search
    unsigned int low  = degree;
    unsigned int high = n + 1;
    return std::distance(knots.begin(), std::upper_bound(knots.begin() + low, knots.begin() + high, u)) - 1;
}

///basis
template<typename T>
std::vector<T> basisFuns(const unsigned int deg, const int i, const std::vector<T> &knots, T u)
{
    std::vector<T> N;
    N.resize(deg + 1, T(0));
    std::vector<T> left, right;
    left.resize(deg + 1, static_cast<T>(0.0));
    right.resize(deg + 1, static_cast<T>(0.0));
    T saved = 0.0, temp = 0.0;

    N[0] = 1.0;

    for (int j = 1; j <= static_cast<int>(deg); j++)
    {
        left[j]  = u - knots[i + 1 - j];
        right[j] = knots[i + j] - u;
        saved    = 0.0;
        for (int r = 0; r < j; r++)
        {
            temp  = N[r] / (right[r + 1] + left[j - r]);
            N[r]  = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
    return N;
}

/// get curve point
template<int dim, typename T>
glm::vec<dim, T> getCurvePoint(unsigned int degree, const std::vector<T> &knots,
                               const std::vector<glm::vec<dim, T>> &control_points, T u)
{
    // Initialize result to 0s
    glm::vec<dim, T> point(T(0));

    // Find span and corresponding non-zero basis functions
    int span         = findSpan(degree, knots, u);
    std::vector<T> N = basisFuns(degree, span, knots, u);

    // Compute point
    for (unsigned int j = 0; j <= degree; j++)
    {
        point += static_cast<T>(N[j]) * control_points[span - degree + j];
    }
    return point;
}

template<typename T>
glm::vec<3, T> rationalCurvePoint(const tinynurbs::RationalCurve<T> &crv, T u)
{
    typedef glm::vec<4, T> tvecnp1;

    // Compute homogenous coordinates of control points
    std::vector<tvecnp1> Cw;
    Cw.reserve(crv.control_points.size());
    for (size_t i = 0; i < crv.control_points.size(); i++)
    {
        Cw.push_back(tvecnp1(tinynurbs::util::cartesianToHomogenous(crv.control_points[i], crv.weights[i])));
    }

    // Compute point using homogenous coordinates
    tvecnp1 pointw = getCurvePoint(crv.degree, crv.knots, Cw, u);

    // Convert back to cartesian coordinates
    return tinynurbs::util::homogenousToCartesian(pointw);
}

/// get surface point
template<int dim, typename T>
glm::vec<dim, T> getSurfacePoint(unsigned int degree_u, unsigned int degree_v,
                                 const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                                 const tinynurbs::array2<glm::vec<dim, T>> &control_points, T u, T v)
{
    // Initialize result to 0s
    glm::vec<dim, T> point(T(0.0));

    // Find span and non-zero basis functions
    int span_u        = findSpan(degree_u, knots_u, u);
    int span_v        = findSpan(degree_v, knots_v, v);
    std::vector<T> Nu = basisFuns(degree_u, span_u, knots_u, u);
    std::vector<T> Nv = basisFuns(degree_v, span_v, knots_v, v);
    unsigned int uind = span_u - degree_u;

    for (int l = 0; l <= degree_v; l++)
    {
        glm::vec<dim, T> temp(0.0);
        unsigned int vind = span_v - degree_v + l;
        for (int k = 0; k <= degree_u; k++)
        {
            temp += static_cast<T>(Nu[k]) * control_points(uind + k, vind);
        }

        point += static_cast<T>(Nv[l]) * temp;
    }
    return point;
}

template<typename T>
glm::vec<3, T> rationalSurfacePoint(const tinynurbs::RationalSurface<T> &srf, T u, T v)
{
    typedef glm::vec<4, T> tvecnp1;

    // Compute homogenous coordinates of control points
    tinynurbs::array2<tvecnp1> Cw;
    Cw.resize(srf.control_points.rows(), srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); i++)
    {
        for (int j = 0; j < srf.control_points.cols(); j++)
        {
            Cw(i, j) =
                    tvecnp1(tinynurbs::util::cartesianToHomogenous(srf.control_points(i, j), srf.weights(i, j)));
        }
    }

    // Compute point using homogenous coordinates
    tvecnp1 pointw = getSurfacePoint(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, Cw, u, v);

    // Convert back to cartesian coordinates
    return tinynurbs::util::homogenousToCartesian(pointw);
}

///derBasis
template<typename T>
tinynurbs::array2<T> dersBasisFuns(unsigned int deg, int span, const std::vector<T> &knots, T u,
                                   int num_ders)
{
    std::vector<T> left, right;
    left.resize(deg + 1, 0.0);
    right.resize(deg + 1, 0.0);
    T saved = 0.0, temp = 0.0;

    tinynurbs::array2<T> ndu(deg + 1, deg + 1);
    ndu(0, 0) = 1.0;

    for (int j = 1; j <= static_cast<int>(deg); j++)
    {
        left[j]  = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        saved    = 0.0;

        for (int r = 0; r < j; r++)
        {
            // Lower triangle
            ndu(j, r) = right[r + 1] + left[j - r];
            temp      = ndu(r, j - 1) / ndu(j, r);
            // Upper triangle
            ndu(r, j) = saved + right[r + 1] * temp;
            saved     = left[j - r] * temp;
        }

        ndu(j, j) = saved;
    }

    tinynurbs::array2<T> ders(num_ders + 1, deg + 1, T(0));

    for (int j = 0; j <= static_cast<int>(deg); j++)
    {
        ders(0, j) = ndu(j, deg);
    }

    tinynurbs::array2<T> a(2, deg + 1);

    for (int r = 0; r <= static_cast<int>(deg); r++)
    {
        int s1  = 0;
        int s2  = 1;
        a(0, 0) = 1.0;

        for (int k = 1; k <= num_ders; k++)
        {
            T d    = 0.0;
            int rk = r - k;
            int pk = deg - k;
            int j1 = 0;
            int j2 = 0;

            if (r >= k)
            {
                a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
                d        = a(s2, 0) * ndu(rk, pk);
            }

            if (rk >= -1)
            {
                j1 = 1;
            }
            else
            {
                j1 = -rk;
            }

            if (r - 1 <= pk)
            {
                j2 = k - 1;
            }
            else
            {
                j2 = deg - r;
            }

            for (int j = j1; j <= j2; j++)
            {
                a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
                d += a(s2, j) * ndu(rk + j, pk);
            }

            if (r <= pk)
            {
                a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
                d += a(s2, k) * ndu(r, pk);
            }

            ders(k, r) = d;

            int temp = s1;
            s1       = s2;
            s2       = temp;
        }
    }

    T fac = static_cast<T>(deg);
    for (int k = 1; k <= num_ders; k++)
    {
        for (int j = 0; j <= static_cast<int>(deg); j++)
        {
            ders(k, j) *= fac;
        }
        fac *= static_cast<T>(deg - k);
    }

    return ders;
}

template<int dim, typename T>
std::vector<glm::vec<dim, T>> getCurveDerivatives(unsigned int degree, const std::vector<T> &knots,
                                                  const std::vector<glm::vec<dim, T>> &control_points,
                                                  int num_ders, T u)
{
    typedef glm::vec<dim, T> tvecn;
    using std::vector;

    std::vector<glm::vec<dim, T>> curve_ders;
    curve_ders.resize(num_ders + 1);

    // Assign higher order derivatives to zero
    for (unsigned int k = degree + 1; k <= num_ders; k++)
    {
        curve_ders[k] = tvecn(0.0);
    }

    // Find the span and corresponding non-zero basis functions & derivatives
    int span                  = findSpan(degree, knots, u);
    tinynurbs::array2<T> ders = dersBasisFuns<T>(degree, span, knots, u, num_ders);

    // Compute first num_ders derivatives
    int du = num_ders < static_cast<int>(degree) ? num_ders : static_cast<int>(degree);
    for (int k = 0; k <= du; k++)
    {
        curve_ders[k] = tvecn(0.0);
        for (int j = 0; j <= static_cast<int>(degree); j++)
        {
            curve_ders[k] += static_cast<T>(ders(k, j)) * control_points[span - degree + j];
        }
    }
    return curve_ders;
}

template<typename T>
std::vector<glm::vec<3, T>> rationalCurveDerivatives(const tinynurbs::RationalCurve<T> &crv, int num_ders, T u)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    std::vector<tvecn> curve_ders;
    curve_ders.reserve(num_ders + 1);

    // Compute homogenous coordinates of control points
    std::vector<tvecnp1> Cw;
    Cw.reserve(crv.control_points.size());
    for (size_t i = 0; i < crv.control_points.size(); i++)
    {
        Cw.push_back(tinynurbs::util::cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
    }

    // Derivatives of Cw
    std::vector<tvecnp1> Cwders = getCurveDerivatives(crv.degree, crv.knots, Cw, num_ders, u);

    // Split Cwders into coordinates and weights
    std::vector<tvecn> Aders;
    std::vector<T> wders;
    for (const auto &val: Cwders)
    {
        Aders.push_back(tinynurbs::util::truncateHomogenous(val));
        wders.push_back(val.w);
    }

    // Compute rational derivatives
    for (int k = 0; k <= num_ders; k++)
    {
        tvecn v = Aders[k];
        for (int i = 1; i <= k; i++)
        {
            v -= static_cast<T>(tinynurbs::util::binomial(k, i)) * wders[i] * curve_ders[k - i];
        }
        curve_ders.push_back(v / wders[0]);
    }
    return curve_ders;
}

template<int dim, typename T>
tinynurbs::array2<glm::vec<dim, T>> getSurfaceDerivatives(unsigned int degree_u, unsigned int degree_v,
                                                          const std::vector<T> &knots_u,
                                                          const std::vector<T> &knots_v,
                                                          const tinynurbs::array2<glm::vec<dim, T>> &control_points,
                                                          unsigned int num_ders, T u, T v)
{
    tinynurbs::array2<glm::vec<dim, T>> surf_ders(num_ders + 1, num_ders + 1, glm::vec<dim, T>(0.0));

    // Set higher order derivatives to 0
    for (unsigned int k = degree_u + 1; k <= num_ders; k++)
    {
        for (unsigned int l = degree_v + 1; l <= num_ders; l++)
        {
            surf_ders(k, l) = glm::vec<dim, T>(0.0);
        }
    }

    // Find span and basis function derivatives
    int span_u                  = findSpan(degree_u, knots_u, u);
    int span_v                  = findSpan(degree_v, knots_v, v);
    tinynurbs::array2<T> ders_u = dersBasisFuns(degree_u, span_u, knots_u, u, num_ders);
    tinynurbs::array2<T> ders_v = dersBasisFuns(degree_v, span_v, knots_v, v, num_ders);

    // Number of non-zero derivatives is <= degree
    unsigned int du = std::min(num_ders, degree_u);
    unsigned int dv = std::min(num_ders, degree_v);

    std::vector<glm::vec<dim, T>> temp;
    temp.resize(degree_v + 1);
    // Compute derivatives
    for (int k = 0; k <= du; k++)
    {
        for (int s = 0; s <= degree_v; s++)
        {
            temp[s] = glm::vec<dim, T>(0.0);
            for (int r = 0; r <= degree_u; r++)
            {
                temp[s] += static_cast<T>(ders_u(k, r)) *
                           control_points(span_u - degree_u + r, span_v - degree_v + s);
            }
        }

        const unsigned int dd = std::min(num_ders - k, dv);

        for (int l = 0; l <= dd; l++)
        {
            for (int s = 0; s <= degree_v; s++)
            {
                surf_ders(k, l) += ders_v(l, s) * temp[s];
            }
        }
    }
    return surf_ders;
}

template<typename T>
tinynurbs::array2<glm::vec<3, T>> rationalSurfaceDerivatives(const tinynurbs::RationalSurface<T> &srf, int num_ders,
                                                             T u, T v)
{
    using namespace std;
    using namespace glm;

    typedef vec<3, T> tvecn;
    typedef vec<4, T> tvecnp1;

    tinynurbs::array2<tvecnp1> homo_cp;
    homo_cp.resize(srf.control_points.rows(), srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); ++i)
    {
        for (int j = 0; j < srf.control_points.cols(); ++j)
        {
            homo_cp(i, j) =
                    tinynurbs::util::cartesianToHomogenous(srf.control_points(i, j), srf.weights(i, j));
        }
    }

    tinynurbs::array2<tvecnp1> homoDers = getSurfaceDerivatives(
            srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, homo_cp, num_ders, u, v);

    tinynurbs::array2<tvecn> Aders;
    Aders.resize(num_ders + 1, num_ders + 1);
    for (int i = 0; i < homoDers.rows(); ++i)
    {
        for (int j = 0; j < homoDers.cols(); ++j)
        {
            Aders(i, j) = tinynurbs::util::truncateHomogenous(homoDers(i, j));
        }
    }

    tinynurbs::array2<tvecn> surfDers(num_ders + 1, num_ders + 1);
    for (int k = 0; k < num_ders + 1; ++k)
    {
        for (int l = 0; l < num_ders - k + 1; ++l)
        {
            auto der = Aders(k, l);

            for (int j = 1; j < l + 1; ++j)
            {
                der -= static_cast<T>(tinynurbs::util::binomial(l, j)) * homoDers(0, j).w * surfDers(k, l - j);
            }

            for (int i = 1; i < k + 1; ++i)
            {
                der -= static_cast<T>(tinynurbs::util::binomial(k, i)) * homoDers(i, 0).w * surfDers(k - i, l);

                tvecn tmp(static_cast<T>(0.0));
                for (int j = 1; j < l + 1; ++j)
                {
                    tmp -= static_cast<T>(tinynurbs::util::binomial(l, j)) * homoDers(i, j).w *
                           surfDers(k - 1, l - j);
                }

                der -= static_cast<T>(tinynurbs::util::binomial(k, i)) * tmp;
            }

            der *= 1 / homoDers(0, 0).w;
            surfDers(k, l) = der;
        }
    }
    return surfDers;
}

template<typename T>
glm::vec<3, T> getSurfaceNormal(const tinynurbs::RationalSurface<T> &srf, T u, T v)
{
    tinynurbs::array2<glm::vec<3, T>> surfaceDers = rationalSurfaceDerivatives(srf, 1, u, v);
    glm::vec<3, T> n                              = glm::cross(surfaceDers(0, 1), surfaceDers(1, 0));

    T n_len = glm::length(n);
    if (!equalTo(n_len, T(0)))
    {
        n /= n_len;
    }
    return n;
}

template<typename T>
void refineKnotVector(const tinynurbs::RationalCurve<T> &curve, std::vector<T> &insertKnotElements,
                      tinynurbs::RationalCurve<T> &result)
{
    typedef glm::vec<4, T> tvecnp1;
    typedef glm::vec<3, T> tvecn;

    int degree                = curve.degree;
    std::vector<T> knotVector = curve.knots;
    // std::vector<XYZW> controlPoints = curve.control_points;
    std::vector<tvecnp1> Cw;
    Cw.reserve(curve.control_points.size());
    for (size_t i = 0; i < curve.control_points.size(); i++)
    {
        Cw.push_back(tvecnp1(tinynurbs::util::cartesianToHomogenous(curve.control_points[i], curve.weights[i])));
    }

    const int n = Cw.size() - 1;
    const int m = n + degree + 1;
    int r       = insertKnotElements.size() - 1;

    const int a = findSpan(degree, knotVector, insertKnotElements[0]);
    const int b = findSpan(degree, knotVector, insertKnotElements[r]) + 1;

    std::vector<T> insertedKnotVector(m + r + 2);
    for (int j = 0; j <= a; j++)
    {
        insertedKnotVector[j] = knotVector[j];
    }
    for (int j = b + degree; j <= m; j++)
    {
        insertedKnotVector[j + r + 1] = knotVector[j];
    }

    std::vector<tvecnp1> updatedControlPoints(n + r + 2);
    for (int j = 0; j <= a - degree; j++)
    {
        updatedControlPoints[j] = Cw[j];
    }
    for (int j = b - 1; j <= n; j++)
    {
        updatedControlPoints[j + r + 1] = Cw[j];
    }

    int i = b + degree - 1;
    int k = b + degree + r;
    for (int j = r; j >= 0; j--)
    {
        while (insertKnotElements[j] <= knotVector[i] && i > a)
        {
            updatedControlPoints[k - degree - 1] = Cw[i - degree - 1];
            insertedKnotVector[k]                = knotVector[i];
            k                                    = k - 1;
            i                                    = i - 1;
        }

        updatedControlPoints[k - degree - 1] = updatedControlPoints[k - degree];
        for (int l = 1; l <= degree; l++)
        {
            int ind = k - degree + l;
            T alpha = insertedKnotVector[k + l] - insertKnotElements[j];
            if (equalTo(abs(alpha), (T) (0.0)))
            {
                updatedControlPoints[ind - 1] = updatedControlPoints[ind];
            }
            else
            {
                alpha                         = alpha / (insertedKnotVector[k + l] - knotVector[i - degree + l]);
                updatedControlPoints[ind - 1] = alpha * updatedControlPoints[ind - 1] +
                                                updatedControlPoints[ind] * (T) (1.0 - alpha);
            }
        }
        insertedKnotVector[k] = insertKnotElements[j];
        k                     = k - 1;
    }

    std::vector<T> updatedWeights;
    std::vector<tvecn> updatedPoints;
    tinynurbs::util::homogenousToCartesian(updatedControlPoints, updatedPoints, updatedWeights);

    result.degree         = degree;
    result.knots          = insertedKnotVector;
    result.control_points = updatedPoints;
    result.weights        = updatedWeights;
}

template<typename T>
void elevateDegree(const tinynurbs::RationalCurve<T> &curve, int times, tinynurbs::RationalCurve<T> &result)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;
    const T MaxDistance = 1E9;

    int degree                = curve.degree;
    std::vector<T> knotVector = curve.knots;
    // std::vector<tvecnp1> controlPoints = curve.ControlPoints;
    std::vector<tvecnp1> controlPoints;
    controlPoints.reserve(curve.control_points.size());
    for (size_t i = 0; i < curve.control_points.size(); i++)
    {
        controlPoints.push_back(
                tvecnp1(tinynurbs::util::cartesianToHomogenous(curve.control_points[i], curve.weights[i])));
    }

    int n   = controlPoints.size() - 1;
    int m   = n + degree + 1;
    int ph  = degree + times;
    int ph2 = floor(ph / 2);

    std::vector<std::vector<T>> bezalfs(degree + times + 1, std::vector<T>(degree + 1));
    bezalfs[0][0] = bezalfs[ph][degree] = 1.0;

    for (int i = 1; i <= ph2; i++)
    {
        T inv   = 1.0 / tinynurbs::util::binomial(ph, i);
        int mpi = std::min(degree, i);

        for (int j = std::max(0, i - times); j <= mpi; j++)
        {
            bezalfs[i][j] = inv * tinynurbs::util::binomial(degree, j) * tinynurbs::util::binomial(times, i - j);
        }
    }

    for (int i = ph2 + 1; i <= ph - 1; i++)
    {
        int mpi = std::min(degree, i);
        for (int j = std::max(0, i - times); j <= mpi; j++)
        {
            bezalfs[i][j] = bezalfs[ph - i][degree - j];
        }
    }

    int mh   = ph;
    int kind = ph + 1;
    int r    = -1;
    int a    = degree;
    int b    = degree + 1;
    int cind = 1;
    T ua     = knotVector[0];

    int moresize = controlPoints.size() + controlPoints.size() * times;
    std::vector<tvecnp1> updatedControlPoints(moresize, tvecnp1(MaxDistance, MaxDistance,
                                                                MaxDistance, 1));
    updatedControlPoints[0] = controlPoints[0];

    std::vector<T> updatedKnotVector(moresize + ph + 1, MaxDistance);
    for (int i = 0; i <= ph; i++)
    {
        updatedKnotVector[i] = ua;
    }

    std::vector<tvecnp1> bpts(degree + 1);
    for (int i = 0; i <= degree; i++)
    {
        bpts[i] = controlPoints[i];
    }

    std::vector<tvecnp1> nextbpts(degree - 1);

    while (b < m)
    {
        int i = b;
        while (b < m && equalTo(knotVector[b], knotVector[b + 1]))
        {
            b = b + 1;
        }
        int mul = b - i + 1;
        mh += mul + times;
        T ub = knotVector[b];

        int oldr = r;
        r        = degree - mul;

        int lbz = oldr > 0 ? floor((oldr + 2) / 2) : 1;
        int rbz = r > 0 ? floor(ph - (r + 1) / 2) : ph;

        if (r > 0)
        {
            T numer = ub - ua;
            std::vector<T> alfs(degree - 1);
            for (int k = degree; k > mul; k--)
            {
                alfs[k - mul - 1] = numer / (knotVector[a + k] - ua);
            }
            for (int j = 1; j <= r; j++)
            {
                int save = r - j;
                int s    = mul + j;

                for (int k = degree; k >= s; k--)
                {
                    bpts[k] = alfs[k - s] * bpts[k] + ((T) 1.0 - alfs[k - s]) * bpts[k - 1];
                }
                nextbpts[save] = bpts[degree];
            }
        }

        std::vector<tvecnp1> ebpts(degree + times + 1);
        for (int i = lbz; i <= ph; i++)
        {
            ebpts[i] = tvecnp1(0.0, 0.0, 0.0, 0.0);
            int mpi  = std::min(degree, i);
            for (int j = std::max(0, i - times); j <= mpi; j++)
            {
                ebpts[i] += bezalfs[i][j] * bpts[j];
            }
        }

        if (oldr > 1)
        {
            int first = kind - 2;
            int last  = kind;
            T den     = ub - ua;
            T bet     = (ub - updatedKnotVector[kind - 1]) / den;

            for (int tr = 1; tr < oldr; tr++)
            {
                int i  = first;
                int j  = last;
                int kj = j - kind + 1;

                while (j - i > tr)
                {
                    if (i < cind)
                    {
                        T alf = (ub - updatedKnotVector[i]) / (ua - updatedKnotVector[i]);
                        updatedControlPoints[i] =
                                alf * updatedControlPoints[i] + ((T) 1.0 - alf) * updatedControlPoints[i - 1];
                    }

                    if (j >= lbz)
                    {
                        if (j - tr <= kind - ph + oldr)
                        {
                            T gam     = (ub - updatedKnotVector[j - tr]) / den;
                            ebpts[kj] = gam * ebpts[kj] + ((T) 1.0 - gam) * ebpts[kj + 1];
                        }
                        else
                        {
                            ebpts[kj] = bet * ebpts[kj] + ((T) 1.0 - bet) * ebpts[kj + 1];
                        }
                    }

                    i  = i + 1;
                    j  = j - 1;
                    kj = kj - 1;
                }

                first -= 1;
                last += 1;
            }
        }

        if (a != degree)
        {
            for (int i = 0; i < ph - oldr; i++)
            {
                updatedKnotVector[kind++] = ua;
            }
        }

        for (int j = lbz; j <= rbz; j++)
        {
            updatedControlPoints[cind++] = ebpts[j];
        }

        if (b < m)
        {
            for (int j = 0; j < r; j++)
            {
                bpts[j] = nextbpts[j];
            }
            for (int j = r; j <= degree; j++)
            {
                bpts[j] = controlPoints[b - degree + j];
            }

            a  = b;
            b  = b + 1;
            ua = ub;
        }
        else
        {
            for (int i = 0; i <= ph; i++)
            {
                updatedKnotVector[kind + i] = ub;
            }
        }
    }

    for (int i = updatedControlPoints.size() - 1; i > 0; i--)
    {
        if (equalTo(updatedControlPoints[i].x, MaxDistance) &&
            equalTo(updatedControlPoints[i].y, MaxDistance) &&
            equalTo(updatedControlPoints[i].z, MaxDistance))
        {
            updatedControlPoints.pop_back();
            continue;
        }
        break;
    }
    for (int i = updatedKnotVector.size() - 1; i > 0; i--)
    {
        if (equalTo(updatedKnotVector[i], MaxDistance))
        {
            updatedKnotVector.pop_back();
            continue;
        }
        break;
    }

    std::vector<T> updatedWeights;
    std::vector<tvecn> updatedPoints;
    tinynurbs::util::homogenousToCartesian(updatedControlPoints, updatedPoints, updatedWeights);

    result.degree         = ph;
    result.knots          = updatedKnotVector;
    result.control_points = updatedPoints;
    result.weights        = updatedWeights;
}

template<typename T>
glm::vec<3, T> createRandomOrthogonal(const glm::vec<3, T> &xyz)
{
    glm::vec<3, T> current = xyz;
    glm::vec<3, T> normal  = glm::normalize(current);
    // glm::vec<3, T> tangent = normal.CrossProduct(XYZ(-normal.z, normal.x, normal.y));
    glm::vec<3, T> tangent = glm::cross(normal, glm::vec<3, T>(-normal.z, normal.x, normal.y));
    // glm::vec<3, T> bitangent = normal.CrossProduct(tangent);
    glm::vec<3, T> bitangent = glm::cross(normal, tangent);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(-glm::pi<T>(), glm::pi<T>());
    double angle = distrib(gen);
    return glm::normalize(tangent * (T) std::sin(angle) + bitangent * (T) std::cos(angle));
}

template<typename T>
bool IsAlmostEqualTo(const glm::vec<3, T> &a, const glm::vec<3, T> &b)
{
    return equalTo(a[0], b[0]) &&
           equalTo(a[1], b[1]) &&
           equalTo(a[2], b[2]);
}

template<typename T>
int getKnotMultiplicity(const std::vector<T> &knotVector, T knot)
{
    int size  = knotVector.size();
    int multi = 0;

    for (int index = 0; index < size; index++)
    {
        if (equalTo(knot, knotVector[index]))
        {
            multi++;
        }
    }

    return multi;
}

template<typename T>
std::map<T, int> getKnotMultiplicityMap(const std::vector<T> &knotVector)
{
    std::map<T, int> result;

    for (int i = 0; i < knotVector.size(); i++)
    {
        T knot   = knotVector[i];
        auto got = result.find(knot);
        if (got == result.end())
        {
            int multi = getKnotMultiplicity(knotVector, knotVector[i]);
            result.insert(std::make_pair(knotVector[i], multi));
        }
    }
    return result;
}

template<typename T>
map<T, int> getInternalKnotMultiplicityMap(const std::vector<T> &knotVector)
{
    auto result = getKnotMultiplicityMap(knotVector);
    if (!result.empty())
    {
        result.erase(result.begin());
        auto lastElementIterator = prev(result.end(), 1);
        result.erase(lastElementIterator);
    }
    return result;
}

template<typename T>
bool isClosed(const tinynurbs::RationalCurve<T> &curve)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    std::vector<T> knotVector = curve.knots;
    T first                   = knotVector[0];
    T end                     = knotVector[knotVector.size() - 1];

    tvecn startPoint = rationalCurvePoint(curve, first);
    tvecn endPoint   = rationalCurvePoint(curve, end);

    // double distance = startPoint.Distance(endPoint);
    T distance   = glm::distance(startPoint, endPoint);
    bool isMatch = equalTo(distance, (T) 0.0);
    if (isMatch) return true;

    // std::vector<tvecnp1> Cw = curve.ControlPoints;
    std::vector<tvecnp1> Cw;
    Cw.reserve(curve.control_points.size());
    for (size_t i = 0; i < curve.control_points.size(); i++)
    {
        Cw.push_back(tvecnp1(tinynurbs::util::cartesianToHomogenous(curve.control_points[i], curve.weights[i])));
    }

    int n = Cw.size() - 1;
    // tvecn last = Cw[n].ToXYZ(true);
    tvecn last = curve.control_points[n];

    bool flag = false;
    int index = 0;
    for (int i = 0; i < n; i++)
    {
        // tvecn current = Cw[i].ToXYZ(true);
        tvecn current = curve.control_points[i];
        // if (last.IsAlmostEqualTo(current))
        if (IsAlmostEqualTo(last, current))
        {
            index = i;
            flag  = true;
            break;
        }
    }

    if (!flag) return false;
    if (index == 0) return true;
    for (int i = index; i >= 0; i--)
    {
        // tvecn current = Cw[i].ToXYZ(true);
        tvecn current = curve.control_points[i];
        // tvecn another = Cw[n - index + i].ToXYZ(true);
        tvecn another = curve.control_points[n - index + i];
        // if (!another.IsAlmostEqualTo(current))
        if (!IsAlmostEqualTo(another, current))
        {
            return false;
        }
    }
    return true;
}

template<typename T>
bool isLinear(const tinynurbs::RationalCurve<T> &curve)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    int degree                       = curve.degree;
    std::vector<T> knotVector        = curve.knots;
    std::vector<tvecn> controlPoints = curve.control_points;

    int size = controlPoints.size();
    if (size == 2)
        return true;
    tvecn start = controlPoints[0];
    tvecn end   = controlPoints[size - 1];
    if (IsAlmostEqualTo(end - start, tvecn(0.0, 0.0, 0.0)))
        return false;
    tvecn dir = glm::normalize(end - start);
    for (int i = 0; i < size - 1; i++)
    {
        tvecn current = controlPoints[i];
        tvecn next    = controlPoints[i + 1];
        if (!IsAlmostEqualTo(glm::normalize(next - current), dir))
        {
            return false;
        }
    }

    auto map = getInternalKnotMultiplicityMap(knotVector);
    for (auto it = map.begin(); it != map.end(); it++)
    {
        T u        = it->first;
        tvecn cp   = rationalCurvePoint(curve, u);
        tvecn cp2s = glm::normalize(cp - start);
        tvecn cp2e = glm::normalize(cp - end);

        if (!IsAlmostEqualTo(cp2s, dir))
            return false;
        if (!IsAlmostEqualTo(cp2e, (T) -1 * dir))
            return false;
    }
    return true;
}

template<typename T>
std::vector<glm::vec<3, T>> projectNormal(const tinynurbs::RationalCurve<T> &curve)
{
    const std::vector<T> &knotVector = curve.knots;
    int size                         = knotVector.size();
    int m                            = size - 1;

    T v0 = knotVector[0];
    T vm = knotVector[m];

    std::vector<glm::vec<3, T>> Blist(size);
    glm::vec<3, T> T0 = glm::normalize(rationalCurveDerivatives(curve, 1, v0)[1]);

    bool flag = true;
    while (flag)
    {
        bool needReCal = false;
        Blist[0]       = createRandomOrthogonal(T0);

        for (int i = 1; i <= m; i++)
        {
            glm::vec<3, T> Ti = glm::normalize(rationalCurveDerivatives(curve, 1, knotVector[i])[1]);

            // needReCal = Ti.CrossProduct(Blist[i - 1]).IsZero();
            needReCal = IsZero(glm::cross(Ti, Blist[i - 1]));

            if (needReCal)
            {
                break;
            }

            glm::vec<3, T> bi = Blist[i - 1] - (glm::dot(Blist[i - 1], Ti)) * Ti;
            // Blist[i] = bi.Normalize();
            Blist[i] = glm::normalize(bi);
        }
        if (!needReCal)
        {
            flag = false;
        }
    }

    bool isClose = isClosed(curve);//
    if (isClose)
    {
        Blist[m] = Blist[0];

        std::vector<glm::vec<3, T>> Baver(m);
        Baver[m] = Blist[0];

        for (int i = m - 1; i >= 1; i--)
        {
            glm::vec<3, T> Ti1 = glm::normalize(rationalCurveDerivatives(curve, 1, knotVector[i + 1])[1]);
            glm::vec<3, T> bi  = Baver[i + 1] + Ti1;
            Baver[i]           = glm::normalize(bi);
        }

        for (int i = 1; i <= m - 1; i++)
        {
            Blist[i] = (T) 0.5 * (Blist[i] + Baver[i]);
        }
    }
    return Blist;
}

template<typename T>
void computeProjectionNormalFrame(const tinynurbs::RationalCurve<T> &curve,
                                  int num_samples,
                                  std::vector<glm::vec<3, T>> &tangents,
                                  std::vector<glm::vec<3, T>> &normals,
                                  std::vector<glm::vec<3, T>> &binormals)
{
    typedef glm::vec<3, T> tvecn;

    // tinynurbs::RationalCurve<T> path;
    // reparametrize(curve, T(0), T(1), path);

    std::vector<T> sample_params(num_samples);
    for (int i = 0; i < num_samples; ++i)
    {
        sample_params[i] = i / T(num_samples - 1);
    }

    tvecn T0 = glm::normalize(rationalCurveDerivatives(curve, 1, sample_params[0])[1]);
    tvecn N0 = createRandomOrthogonal(T0);
    tvecn B0 = glm::normalize(glm::cross(T0, N0));
    N0       = glm::normalize(glm::cross(B0, T0));

    tangents.push_back(T0);
    normals.push_back(N0);
    binormals.push_back(B0);

    for (int i = 1; i < num_samples; ++i)
    {
        tvecn Ti = glm::normalize(rationalCurveDerivatives(curve, 1, sample_params[i])[1]);

        tvecn Ni = normals.back() - glm::dot(normals.back(), Ti) * Ti;
        Ni       = glm::normalize(Ni);

        tvecn Bi = glm::normalize(glm::cross(Ti, Ni));

        tangents.push_back(Ti);
        normals.push_back(Ni);
        binormals.push_back(Bi);
    }
}

template<typename T>
void computeFrenetFrames(const tinynurbs::RationalCurve<T> &curve,
                         int num_samples,
                         std::vector<glm::vec<3, T>> &tangents,
                         std::vector<glm::vec<3, T>> &normals,
                         std::vector<glm::vec<3, T>> &binormals)
{
    typedef glm::vec<3, T> tvecn;

    const std::vector<T> &knots = curve.knots;
    T t_start                   = knots.front();
    T t_end                     = knots.back();

    std::vector<T> samples(num_samples);
    for (int i = 0; i < num_samples; ++i)
    {
        samples[i] = t_start + i * (t_end - t_start) / (num_samples - 1);
    }

    for (T t: samples)
    {
        tvecn tangent, normal, binormal;

        auto derivatives = rationalCurveDerivatives(curve, 2, t);

        tangent = glm::normalize(derivatives[1]);

        binormal = glm::normalize(glm::cross(derivatives[1], derivatives[2]));

        // normal = glm::normalize(derivatives[2] - tangent * glm::dot(derivatives[2], tangent));
        normal = glm::normalize(glm::cross(binormal, tangent));
        // binormal = glm::normalize(glm::cross(tangent, normal));

        tangents.push_back(tangent);
        normals.push_back(normal);
        binormals.push_back(binormal);
    }
}

template<typename T>
void computeRMF(const tinynurbs::RationalCurve<T> &curve,
                int num_samples,
                std::vector<glm::vec<3, T>> &tangents,
                std::vector<glm::vec<3, T>> &normals,
                std::vector<glm::vec<3, T>> &binormals)
{
    typedef glm::vec<3, T> tvecn;

    // tinynurbs::RationalCurve<T> path;
    // reparametrize(curve, T(0), T(1), path);

    std::vector<T> sample_params(num_samples);
    for (int i = 0; i < num_samples; ++i)
    {
        sample_params[i] = i / T(num_samples - 1);
    }

    tvecn T0 = glm::normalize(rationalCurveDerivatives(curve, 1, sample_params[0])[1]);
    tvecn N0 = glm::vec3(0, 1, 0);
    if (glm::dot(N0, T0) > 0.9)
    {
        N0 = glm::vec3(1, 0, 0);
    }
    tvecn B0 = glm::normalize(glm::cross(T0, N0));
    N0       = glm::normalize(glm::cross(B0, T0));

    tangents.push_back(T0);
    normals.push_back(N0);
    binormals.push_back(B0);

    for (int i = 1; i < num_samples; ++i)
    {
        tvecn Ti = glm::normalize(rationalCurveDerivatives(curve, 1, sample_params[i])[1]);

        tvecn Ni = normals.back() - glm::dot(Ti, normals.back()) * Ti;
        Ni       = glm::normalize(Ni);

        tvecn Bi = glm::normalize(glm::cross(Ti, Ni));

        tangents.push_back(Ti);
        normals.push_back(Ni);
        binormals.push_back(Bi);
    }
}

template<typename T>
T getTotalChordLength(const std::vector<glm::vec<3, T>> &throughPoints)
{
    int n    = throughPoints.size() - 1;
    T length = 0.0;
    for (int i = 1; i <= n; i++)
    {
        length += glm::distance(throughPoints[i], throughPoints[i - 1]);
    }
    return length;
}

template<typename T>
std::vector<T> getChordParameterization(const std::vector<glm::vec<3, T>> &throughPoints)
{
    int size = throughPoints.size();
    int n    = size - 1;

    std::vector<T> uk(size, 0.0);
    uk[n] = 1.0;

    T d = getTotalChordLength(throughPoints);
    for (int i = 1; i <= n - 1; i++)
    {
        uk[i] = uk[i - 1] + (glm::distance(throughPoints[i], throughPoints[i - 1])) / d;
    }
    return uk;
}

template<typename T>
std::vector<std::vector<T>> solveLinearSystem(const std::vector<std::vector<T>> &matrix,
                                              const std::vector<std::vector<T>> &right)
{
    std::vector<std::vector<T>> result(matrix.size(), std::vector<T>(right[0].size()));

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m(matrix.size(), matrix[0].size());
    for (int i = 0; i < matrix.size(); i++)
    {
        m.row(i) = Eigen::Matrix<T, Eigen::Dynamic, 1>::Map(&matrix[i][0], matrix[i].size());
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> r(right.size(), right[0].size());
    for (int i = 0; i < right.size(); i++)
    {
        r.row(i) = Eigen::Matrix<T, Eigen::Dynamic, 1>::Map(&right[i][0], right[i].size());
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> solve = m.lu().solve(r);
    for (int col = 0; col < solve.cols(); ++col)
    {
        for (int row = 0; row < solve.rows(); ++row)
        {
            result[row][col] = solve(row, col);
        }
    }
    return result;
}

template<typename T>
std::vector<T> averageKnotVector(int degree, const std::vector<T> &params)
{
    std::vector<T> uk = params;
    int size          = params.size();
    int n             = size - 1;
    int m             = n + degree + 1;

    std::vector<T> knotVector(m + 1, 0.0);
    for (int i = m - degree; i <= m; i++)
    {
        knotVector[i] = 1.0;
    }

    for (int j = 1; j <= n - degree; j++)
    {
        T sum = 0.0;
        for (int i = j; i <= j + degree - 1; i++)
        {
            sum += uk[i];
        }
        knotVector[j + degree] = (1.0 / degree) * sum;
    }
    return knotVector;
}

template<typename T>
void globalInterpolation(unsigned int degree, const std::vector<glm::vec<3, T>> &throughPoints,
                         tinynurbs::RationalCurve<T> &curve,
                         const std::vector<T> &params)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;
    const int NURBSMaxDegree = 7;
    // VALIDATE_ARGUMENT(degree >= 0 && degree <= Constants::NURBSMaxDegree, "degree",
    //                   "Degree must be greater than or equal zero and not exceed the maximun degree.");
    // VALIDATE_ARGUMENT(throughPoints.size() > degree, "throughPoints",
    //                   "ThroughPoints size must be greater than degree.");
    int size = throughPoints.size();
    int n    = size - 1;

    std::vector<T> uk(size);
    if (params.empty())
    {
        uk = getChordParameterization(throughPoints);
    }
    else
    {
        // VALIDATE_ARGUMENT(params.size() == size, "params", "Params size must be equal to throughPoints size.");
        uk = params;
    }
    std::vector<T> knotVector = averageKnotVector(degree, uk);

    std::vector<std::vector<T>> A(size, std::vector<T>(size));
    for (int i = 1; i < n; i++)
    {
        int spanIndex = findSpan(degree, knotVector, uk[i]);
        // T basis[NURBSMaxDegree + 1];
        std::vector<T> basis(NURBSMaxDegree + 1);
        // Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i], basis);
        basis = basisFuns(degree, spanIndex, knotVector, uk[i]);

        for (int j = 0; j <= degree; j++)
        {
            A[i][spanIndex - degree + j] = basis[j];
        }
    }
    A[0][0] = 1.0;
    A[n][n] = 1.0;

    std::vector<std::vector<T>> right(size, std::vector<T>(3));
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            right[i][j] = throughPoints[i][j];
        }
    }

    std::vector<tvecnp1> controlPoints(size);
    std::vector<std::vector<T>> result = solveLinearSystem(A, right);
    for (int i = 0; i < result.size(); i++)
    {
        tvecn temp = tvecn(0, 0, 0);
        for (int j = 0; j < 3; j++)
        {
            temp[j] = result[i][j];
        }
        controlPoints[i] = tvecnp1(temp, 1.0);
    }

    std::vector<T> updatedWeights;
    std::vector<tvecn> updatedPoints;
    tinynurbs::util::homogenousToCartesian(controlPoints, updatedPoints, updatedWeights);

    curve.degree         = degree;
    curve.knots          = knotVector;
    curve.control_points = updatedPoints;
    curve.weights        = updatedWeights;
}

template<typename T>
void CreateLoftSurface(const std::vector<tinynurbs::RationalCurve<T>> &sections,
                       tinynurbs::RationalSurface<T> &surface, int customTrajectoryDegree = 0,
                       const std::vector<T> &customTrajectoryKnotVector = {})
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    int degree_max = 0;
    for (int k = 0; k < sections.size(); k++)
    {
        tinynurbs::RationalCurve<T> current = sections[k];
        if (degree_max < current.degree)
        {
            degree_max = current.degree;
        }
    }

    int size = sections.size();

    std::vector<tinynurbs::RationalCurve<T>> internals(size);
    std::vector<std::vector<glm::vec<4, T>>> curvesControlPoints(size);
    for (int k = 0; k < size; k++)
    {
        tinynurbs::RationalCurve<T> current = sections[k];
        if (degree_max > current.degree)
        {
            int times = degree_max - current.degree;
            tinynurbs::RationalCurve<T> tc;
            elevateDegree(current, times, tc);//
            current.degree         = degree_max;
            current.knots          = tc.knots;
            current.control_points = tc.control_points;
            current.weights        = tc.weights;
        }
        // curvesControlPoints[k] = current.ControlPoints;
        curvesControlPoints[k] = tinynurbs::util::cartesianToHomogenous(current.control_points, current.weights);

        internals[k] = current;
    }

    int degreeU = degree_max;
    int degreeV = degree_max;

    std::vector<T> vl(size);
    if (customTrajectoryKnotVector.empty())
    {
        int K = size - 1;
        vl[0] = 0;
        vl[K] = 1;
        for (int k = 1; k <= K - 1; k++)
        {
            std::vector<tvecnp1> current = curvesControlPoints[k];
            int n                        = current.size() - 1;
            T sum                        = 0.0;
            for (int i = 0; i <= n; i++)
            {
                T delta = glm::distance(tinynurbs::util::homogenousToCartesian(curvesControlPoints[k][i]),
                                        (tinynurbs::util::homogenousToCartesian(curvesControlPoints[k - 1][i])));

                std::vector<tvecn> tempPoints(size);
                for (int j = 0; j <= K; j++)
                {
                    tempPoints[j] = tinynurbs::util::homogenousToCartesian(curvesControlPoints[j][i]);
                }

                T di = getTotalChordLength(tempPoints);
                sum += delta / di;
            }

            vl[k] = vl[k - 1] + (1.0 / (n + 1)) * sum;
        }
    }
    else
    {
        degreeV = customTrajectoryDegree;
        vl      = customTrajectoryKnotVector;
    }
    std::vector<T> knotVectorV = averageKnotVector(degreeV, vl);
    std::vector<std::vector<tvecn>> controlPoints;
    std::vector<std::vector<T>> weights;
    // tinynurbs::array2<glm::vec<3, T>> controlPoints;
    int column = curvesControlPoints[0].size();
    for (int c = 0; c < column; c++)
    {
        std::vector<tvecn> temp(size);
        std::vector<T> weightCache(size);
        for (int k = 0; k < size; k++)
        {
            tvecnp1 cp     = curvesControlPoints[k][c];
            T w            = cp.w;
            weightCache[k] = w;
            temp[k]        = tvecn(cp.x / w, cp.y / w, cp.z / w);
        }
        tinynurbs::RationalCurve<T> tc;
        globalInterpolation(degreeV, temp, tc, vl);
        for (int i = 0; i < tc.control_points.size(); i++)
        {
            // XYZW cp = tc.ControlPoints[i];
            // tvecnp1 cp = tinynurbs::util::cartesianToHomogenous(tc.control_points[i], tc.weights[i]);

            // T w = cp.w;
            // tc.ControlPoints[i] = XYZW(XYZ(cp.GetWX() / w, cp.GetWY() / w, cp.GetWz / w), weightCache[i]);
            tc.control_points[i] /= weightCache[i];
            tc.weights[i] = weightCache[i];
        }
        controlPoints.emplace_back(tc.control_points);
        weights.emplace_back(tc.weights);
    }
    std::vector<tvecn> controlPoints1D;
    for (const auto &row: controlPoints)
    {
        for (tvecn val: row)
        {
            controlPoints1D.push_back(val);
        }
    }

    std::vector<T> weights1D;
    for (const auto &row: weights)
    {
        for (T val: row)
        {
            weights1D.push_back(val);
        }
    }

    surface.degree_u       = degreeU;
    surface.degree_v       = degreeV;
    surface.knots_u        = internals[0].knots;
    surface.knots_v        = knotVectorV;
    surface.control_points = tinynurbs::array2<tvecn>(controlPoints.size(), controlPoints[0].size(), controlPoints1D);
    surface.weights        = tinynurbs::array2<T>(weights.size(), weights[0].size(), weights1D);
}

template<typename T>
glm::vec<3, T> worldToLocal(const glm::vec<3, T> &localOrigin, const glm::vec<3, T> &localXdir,
                            const glm::vec<3, T> &localYdir, const glm::vec<3, T> &localZdir,
                            const glm::vec<3, T> &worldPoint)
{
    glm::vec<3, T> traslation = (worldPoint - localOrigin);
    return {glm::dot(traslation, localXdir), glm::dot(traslation, localYdir), glm::dot(traslation, localZdir)};
}

template<typename T>
glm::vec<3, T> LocalToWorld(const glm::vec<3, T> &localOrigin, const glm::vec<3, T> &localXdir,
                            const glm::vec<3, T> &localYdir, const glm::vec<3, T> &localZdir,
                            const glm::vec<3, T> &worldPoint)
{
    return worldPoint.x * localXdir + worldPoint.y * localYdir + worldPoint.z * localZdir + localOrigin;
}

template<typename T>
void createSweepSurface1(const tinynurbs::RationalCurve<T> &profile, const tinynurbs::RationalCurve<T> &trajectory,
                         int minimumProfiles = 0, tinynurbs::RationalSurface<T> &surface = {},
                         std::vector<T> &v_bar                           = {},
                         std::vector<std::vector<glm::vec<3, T>>> &frame = {}, tinynurbs::RationalCurve<T> &sect = {})
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    // Determine number of sections for lofting.
    const int q                         = trajectory.degree;
    std::vector<T> trajectoryKnotVector = trajectory.knots;
    const int ktv                       = trajectoryKnotVector.size();
    const int K                         = minimumProfiles;
    int nsect                           = K + 1;

    tinynurbs::RationalCurve<T> trajectoryCopy = trajectory;
    if (ktv <= nsect + q)
    {
        // Insert knots to widest spans
        // to reach the required minimum number of sections.
        int m                           = nsect + q - ktv + 1;
        std::vector<T> insertedElements = getInsertedKnotElements(m, trajectoryKnotVector);
        refineKnotVector(trajectory, insertedElements, trajectoryCopy);
    }
    else
    {
        if (ktv > nsect + q + 1)
        {
            nsect = ktv - q - 1;
        }
    }

    // Determine lofting profile positions.
    std::vector<T> knotVectorV = trajectoryCopy.knots;
    auto minmax                = std::minmax_element(knotVectorV.begin(), knotVectorV.end());

    std::vector<T> vk(nsect);
    vk[0]         = *minmax.first;
    vk[nsect - 1] = *minmax.second;

    for (int k = 1; k < nsect - 1; k++)
    {
        T sum = 0.0;
        for (int i = 1; i <= q; i++)
        {
            sum += knotVectorV[k + i];
        }
        vk[k] = sum / (T) q;
    }

    // Compute trajectory normals.
    std::vector<tvecn> Bv = projectNormal(trajectoryCopy);

    // const std::vector<tvecnp1> &profileControlPoints = profile.ControlPoints;
    std::vector<tvecnp1> profileControlPoints;
    profileControlPoints.reserve(profile.control_points.size());
    for (size_t i = 0; i < profile.control_points.size(); i++)
    {
        profileControlPoints.push_back(
                tvecnp1(tinynurbs::util::cartesianToHomogenous(profile.control_points[i], profile.weights[i])));
    }

    int profileCpSize = profileControlPoints.size();

    // for each lofting section
    std::vector<tinynurbs::RationalCurve<T>> sections(nsect);
    std::vector<tvecnp1> transformedControlPoints(profileCpSize);
    for (int k = 0; k < nsect; k++)
    {
        // Build local coordinate system for the section.
        tvecn zdir             = glm::normalize(rationalCurveDerivatives(trajectoryCopy, 1, vk[k])[1]);
        unsigned int spanIndex = findSpan(trajectoryCopy.degree, trajectoryCopy.knots, vk[k]);
        const tvecn &xdir      = Bv[spanIndex];
        // tvecn ydir = zdir.CrossProduct(xdir).Normalize();
        tvecn ydir = glm::normalize(glm::cross(zdir, xdir));

        ///
        frame.push_back({xdir, ydir, zdir});

        // Transform the input profile to the lofting position.

        if (k == 0)
            for (int i = 0; i < profileCpSize; i++)
            {
                tvecnp1 wp = profileControlPoints[i];
                T w        = wp.w;
                // tvecn p = wp.ToXYZ(true);
                tvecn p = tinynurbs::util::homogenousToCartesian(wp);

                // tvecn transformed = worldToLocal(p, xdir, ydir, zdir, rationalCurvePoint(trajectoryCopy, vk[k]));
                tvecn transformed           = worldToLocal(rationalCurvePoint(trajectoryCopy, vk[k]), xdir, ydir, zdir, p);
                transformedControlPoints[i] = tvecnp1(transformed, w);
            }

        std::vector<tvecn> transformedPoints;
        std::vector<T> weights;
        tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformedPoints, weights);

        // Make section curve.
        // note: std::vector::swap relinks pointers, which is faster than copying.
        tinynurbs::RationalCurve<T> &section = sections[k];
        section.degree                       = profile.degree;
        section.knots                        = profile.knots;
        section.weights.swap(weights);
        section.control_points.swap(transformedPoints);

        for (int j = 0; j < profileCpSize; j++)
        {
            section.control_points[j] =
                    LocalToWorld(rationalCurvePoint(trajectoryCopy, vk[k]), xdir, ydir, zdir,
                                 section.control_points[j]);
        }
    }

    v_bar = vk;
    sect  = sections[2];

    // Do the lofting.
    CreateLoftSurface(sections, surface);
}

template<typename T>
std::vector<T> rescale(const std::vector<T> &knotVector, T min, T max)
{
    T origintMin = knotVector.front();
    T origintMax = knotVector.back();
    T k          = (max - min) / (origintMax - origintMin);

    int size = knotVector.size();
    std::vector<T> result(size);
    for (int i = 0; i < size; i++)
    {
        result[i] = (k * knotVector[i] - origintMin) + min;
    }
    return result;
}

template<typename T>
void reparametrize(const tinynurbs::RationalCurve<T> &curve, T min, T max, tinynurbs::RationalCurve<T> &result)
{
    std::vector<T> knotVector = curve.knots;

    if (equalTo(min, knotVector.front()) && equalTo(max, knotVector.back()))
    {
        result = curve;
        return;
    }

    std::vector<T> newKnotVector = rescale(knotVector, min, max);
    result                       = curve;
    result.knots                 = newKnotVector;
}

template<typename T>
std::vector<tinynurbs::RationalCurve<T>> decomposeToBeziers(const tinynurbs::RationalCurve<T> &curve)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    unsigned int degree              = curve.degree;
    std::vector<T> knotVector        = curve.knots;
    std::vector<tvecn> controlPoints = curve.control_points;
    std::vector<T> weights           = curve.weights;

    unsigned int knotSize = 2 * (degree + 1);
    std::vector<T> bezierKnots(knotSize);
    for (int i = 0; i < knotSize / 2; i++)
    {
        bezierKnots[i] = 0;
    }
    for (int i = knotSize / 2; i < knotSize; i++)
    {
        bezierKnots[i] = 1;
    }
    int bezierSize = controlPoints.size() - degree;
    std::vector<tinynurbs::RationalCurve<T>> beziers(bezierSize);
    for (int i = 0; i < bezierSize; i++)
    {
        beziers[i].degree         = degree;
        beziers[i].knots          = bezierKnots;
        beziers[i].control_points = std::vector<tvecn>(degree + 1);
        beziers[i].weights        = std::vector<T>(degree + 1);
    }

    int n = controlPoints.size() - 1;
    int m = n + degree + 1;

    unsigned int a = degree;
    unsigned int b = degree + 1;

    int nb = 0;
    for (int i = 0; i <= degree; i++)
    {
        beziers[nb].control_points[i] = controlPoints[i];
        beziers[nb].weights[i]        = weights[i];
    }

    while (b < m)
    {
        int i = b;
        while (b < m && equalTo(knotVector[b + 1], knotVector[b]))
        {
            b++;
        }
        int multi = b - i + 1;
        if (multi < degree)
        {
            T numerator = knotVector[b] - knotVector[a];
            std::vector<T> alphaVector(degree + 1);
            for (int j = degree; j > multi; j--)
            {
                alphaVector[j - multi - 1] = numerator / (knotVector[a + j] - knotVector[a]);
            }

            int r = degree - multi;
            for (int j = 1; j <= r; j++)
            {
                int save = r - j;
                int s    = multi + j;
                for (int k = degree; k >= s; k--)
                {
                    T alpha = alphaVector[k - s];
                    beziers[nb].control_points[k] =
                            (T) alpha * beziers[nb].control_points[k] + (T) (1.0 - alpha) * beziers[nb].control_points[k - 1];
                    beziers[nb].weights[k] = (T) alpha * beziers[nb].weights[k] + (T) (1.0 - alpha) *
                                                                                          beziers[nb].weights[k - 1];
                }

                if (b < m)
                {
                    beziers[nb + 1].control_points[save] = beziers[nb].control_points[degree];
                    beziers[nb + 1].weights[save]        = beziers[nb].weights[degree];
                }
            }

            nb++;
            if (b < m)
            {
                for (int i = degree - multi; i <= degree; i++)
                {
                    beziers[nb].control_points[i] = controlPoints[b - degree + i];
                    beziers[nb].weights[i]        = weights[b - degree + i];
                }

                a = b;
                b += 1;
            }
        }
    }
    return beziers;
}

template<typename T>
T approximateLength(const tinynurbs::RationalCurve<T> &curve)
{
    const std::vector<T> GaussLegendreAbscissae =
            {
                    -0.0640568928626056260850430826247450385909,
                    0.0640568928626056260850430826247450385909,
                    -0.1911188674736163091586398207570696318404,
                    0.1911188674736163091586398207570696318404,
                    -0.3150426796961633743867932913198102407864,
                    0.3150426796961633743867932913198102407864,
                    -0.4337935076260451384870842319133497124524,
                    0.4337935076260451384870842319133497124524,
                    -0.5454214713888395356583756172183723700107,
                    0.5454214713888395356583756172183723700107,
                    -0.6480936519369755692524957869107476266696,
                    0.6480936519369755692524957869107476266696,
                    -0.7401241915785543642438281030999784255232,
                    0.7401241915785543642438281030999784255232,
                    -0.8200019859739029219539498726697452080761,
                    0.8200019859739029219539498726697452080761,
                    -0.8864155270044010342131543419821967550873,
                    0.8864155270044010342131543419821967550873,
                    -0.9382745520027327585236490017087214496548,
                    0.9382745520027327585236490017087214496548,
                    -0.9747285559713094981983919930081690617411,
                    0.9747285559713094981983919930081690617411,
                    -0.9951872199970213601799974097007368118745,
                    0.9951872199970213601799974097007368118745,
            };

    const std::vector<T> GaussLegendreWeights =
            {
                    0.1279381953467521569740561652246953718517,
                    0.1279381953467521569740561652246953718517,
                    0.1258374563468282961213753825111836887264,
                    0.1258374563468282961213753825111836887264,
                    0.121670472927803391204463153476262425607,
                    0.121670472927803391204463153476262425607,
                    0.1155056680537256013533444839067835598622,
                    0.1155056680537256013533444839067835598622,
                    0.1074442701159656347825773424466062227946,
                    0.1074442701159656347825773424466062227946,
                    0.0976186521041138882698806644642471544279,
                    0.0976186521041138882698806644642471544279,
                    0.086190161531953275917185202983742667185,
                    0.086190161531953275917185202983742667185,
                    0.0733464814110803057340336152531165181193,
                    0.0733464814110803057340336152531165181193,
                    0.0592985849154367807463677585001085845412,
                    0.0592985849154367807463677585001085845412,
                    0.0442774388174198061686027482113382288593,
                    0.0442774388174198061686027482113382288593,
                    0.0285313886289336631813078159518782864491,
                    0.0285313886289336631813078159518782864491,
                    0.0123412297999871995468056670700372915759,
                    0.0123412297999871995468056670700372915759,
            };

    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    if (isLinear(curve))
    {
        std::vector<tvecn> controlPoints = curve.control_points;
        tvecn startPoint                 = controlPoints.front();
        tvecn endPoint                   = controlPoints.back();
        return glm::distance(startPoint, endPoint);
    }

    tinynurbs::RationalCurve<T> reCurve;
    reparametrize(curve, (T) 0.0, (T) 1.0, reCurve);

    unsigned int degree              = reCurve.degree;
    std::vector<T> knotVector        = reCurve.knots;
    std::vector<tvecn> controlPoints = reCurve.control_points;///4

    T length                                              = 0.0;
    std::vector<tinynurbs::RationalCurve<T>> bezierCurves = decomposeToBeziers(reCurve);
    for (int i = 0; i < bezierCurves.size(); i++)
    {
        tinynurbs::RationalCurve<T> bezierCurve = bezierCurves[i];

        std::vector<T> bKnots = bezierCurve.knots;
        T a                   = bKnots[0];
        T b                   = bKnots[bKnots.size() - 1];
        T coefficient         = (b - a) / 2.0;

        T bLength                = 0.0;
        std::vector<T> abscissae = GaussLegendreAbscissae;
        int size                 = abscissae.size();
        for (int i = 0; i < size; i++)
        {
            T t         = coefficient * abscissae[i] + (a + b) / 2.0;
            T derLength = glm::length(rationalCurveDerivatives(bezierCurve, 1, t)[1]);
            if (std::isnan(derLength))
                derLength = 0.0;
            bLength += GaussLegendreWeights[i] * derLength;
        }
        bLength = coefficient * bLength;
        length += bLength;
    }
    return length;
}

template<typename T>
bool splitAt(const tinynurbs::RationalCurve<T> &curve, T parameter, tinynurbs::RationalCurve<T> &left,
             tinynurbs::RationalCurve<T> &right)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    int degree                       = curve.degree;
    std::vector<T> knotVector        = curve.knots;
    std::vector<tvecn> controlPoints = curve.control_points;
    std::vector<T> weights           = curve.weights;

    if (IsLessThanOrEqual(parameter, knotVector[degree]) ||
        IsGreaterThanOrEqual(parameter, knotVector[knotVector.size() - degree - 1]))
    {
        return false;
    }

    int multi = getKnotMultiplicity(knotVector, parameter);
    std::vector<T> insert(degree + 1 - multi, parameter);
    left = curve;
    if (insert.size() > 0)
    {
        tinynurbs::RationalCurve<T> temp;
        refineKnotVector(left, insert, temp);
        left = temp;
    }

    int spanIndex      = findSpan(left.degree, left.knots, parameter) - degree;
    right.degree       = degree;
    int rControlPoints = left.control_points.size() - spanIndex;
    std::vector<tvecn> rightControlPoints(rControlPoints);
    std::vector<T> rightWeights(rControlPoints);
    std::vector<T> rightKnotVector(rControlPoints + degree + 1);
    right.knots          = rightKnotVector;
    right.control_points = rightControlPoints;
    right.weights        = rightWeights;

    for (int i = left.control_points.size() - 1, j = rControlPoints - 1; j >= 0; j--, i--)
    {
        right.control_points[j] = left.control_points[i];
        right.weights[j]        = left.weights[i];
    }

    for (int i = left.knots.size() - 1, j = rControlPoints + degree; j >= 0; j--, i--)
    {
        right.knots[j] = left.knots[i];
    }

    left.degree = degree;
    left.control_points.resize(spanIndex);
    left.weights.resize(spanIndex);
    left.knots.resize(spanIndex + degree + 1);

    return true;
}

template<typename T>
T getParamByLength(const tinynurbs::RationalCurve<T> &curve, T start, T end, T givenLength)
{
    const T DistanceEpsilon = 1E-4;

    T middle = (start + end) / 2.0;

    tinynurbs::RationalCurve<T> left;
    tinynurbs::RationalCurve<T> right;
    bool isSplited = splitAt(curve, middle, left, right);
    if (!isSplited) return start;
    T length = approximateLength(left);
    if (equalTo(length, givenLength, DistanceEpsilon))
    {
        return middle;
    }
    else if (IsGreaterThan(length, givenLength, DistanceEpsilon))
    {
        end = middle;
        return getParamByLength(curve, start, end, givenLength);
    }
    else if (IsLessThan(length, givenLength, DistanceEpsilon))
    {
        start = middle;
        return getParamByLength(curve, start, end, givenLength);
    }
    return middle;
}

template<typename T>
T GetParamOnCurve(const tinynurbs::RationalCurve<T> &curve, T givenLength)
{
    const T DistanceEpsilon = 1E-4;

    std::vector<T> knotVector = curve.knots;
    T start                   = knotVector.front();
    T end                     = knotVector.back();

    T totalLength = approximateLength(curve);
    if (IsLessThan(totalLength, givenLength, DistanceEpsilon))
    {
        return end;
    }

    for (int i = 0; i < knotVector.size(); i++)
    {
        tinynurbs::RationalCurve<T> left;
        tinynurbs::RationalCurve<T> right;
        T knot         = knotVector[i];
        bool isSplited = splitAt(curve, knot, left, right);
        if (!isSplited) continue;
        T length = approximateLength(left);
        if (equalTo(length, givenLength))
        {
            return knot;
        }
        if (IsGreaterThan(length, givenLength))
        {
            end = knot;
            break;
        }
    }
    T param = getParamByLength(curve, start, end, givenLength);
    return param;
}

template<typename T>
std::vector<T> getParamsOnCurve(const tinynurbs::RationalCurve<T> &curve, T givenLength)
{
    std::vector<T> result;

    std::vector<T> knotVector = curve.knots;
    T end                     = knotVector.back();

    T param = GetParamOnCurve(curve, givenLength);
    while (!equalTo(param, end))
    {
        result.emplace_back(param);

        tinynurbs::RationalCurve<T> left;
        tinynurbs::RationalCurve<T> right;
        bool isSplited = splitAt(curve, param, left, right);
        if (!isSplited) continue;
        param = GetParamOnCurve(right, givenLength);
    }
    return result;
}

template<typename T>
void createSweepSurface2(const tinynurbs::RationalCurve<T> &profile, const tinynurbs::RationalCurve<T> &trajectory,
                         int minimumProfiles = 0, int customTrajectoryDegree = 0,
                         tinynurbs::RationalSurface<T> &surface = {}, std::vector<T> &v_bar = {},
                         std::vector<std::vector<glm::vec<3, T>>> &frame = {})
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    int K              = minimumProfiles;
    int nsect          = K + 1;
    T trajectoryLength = approximateLength(trajectory);
    T stepLength       = trajectoryLength / (nsect - 1);
    tinynurbs::RationalCurve<T> path;
    reparametrize(trajectory, (T) 0, (T) 1, path);

    std::vector<T> segments(nsect);
    segments[0]           = 0.0;
    std::vector<T> params = getParamsOnCurve(path, stepLength);
    for (int i = 0; i < params.size(); i++)
    {
        segments[i + 1] = params[i];
    }

    v_bar = segments;

    std::vector<tvecn> Bv = projectNormal(path);
    // const std::vector<tvecn> &profileControlPoints = profile.ControlPoints;
    // const std::vector<T> &profileWeights = profile.weights;

    std::vector<tvecnp1> profileControlPoints = tinynurbs::util::cartesianToHomogenous(profile.control_points,
                                                                                       profile.weights);
    int profileCpSize                         = profileControlPoints.size();

    std::vector<tinynurbs::RationalCurve<T>> sections(nsect);
    std::vector<tvecnp1> transformedControlPoints(profileCpSize);

    for (int k = 0; k < nsect; k++)
    {
        tvecn zdir    = glm::normalize(rationalCurveDerivatives(path, 1, segments[k])[1]);
        int spanIndex = findSpan(path.degree, path.knots, segments[k]);
        tvecn xdir    = Bv[spanIndex];
        tvecn ydir    = glm::normalize(glm::cross(zdir, xdir));

        //frame
        frame.push_back({xdir, ydir, zdir});

        if (k == 0)
            for (int i = 0; i < profileCpSize; i++)
            {
                tvecnp1 wp = profileControlPoints[i];
                T w        = wp.w;
                // tvecn p = wp.ToXYZ(true);
                tvecn p = tinynurbs::util::homogenousToCartesian(wp);

                // tvecn transformed = worldToLocal(p, xdir, ydir, zdir, rationalCurvePoint(trajectoryCopy, vk[k]));
                tvecn transformed           = worldToLocal(rationalCurvePoint(path, segments[k]), xdir, ydir, zdir, p);
                transformedControlPoints[i] = tvecnp1(transformed, w);
            }

        std::vector<tvecn> transformedPoints;
        std::vector<T> weights;
        tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformedPoints, weights);

        // Make section curve.
        // note: std::vector::swap relinks pointers, which is faster than copying.
        tinynurbs::RationalCurve<T> &section = sections[k];
        section.degree                       = profile.degree;
        section.knots                        = profile.knots;
        section.weights.swap(weights);
        section.control_points.swap(transformedPoints);

        for (int j = 0; j < profileCpSize; j++)
        {
            section.control_points[j] =
                    LocalToWorld(rationalCurvePoint(path, segments[k]), xdir, ydir, zdir,
                                 section.control_points[j]);
        }
    }

    CreateLoftSurface(sections, surface, customTrajectoryDegree, segments);
}

template<typename T>
bool LeastSquaresApproximation(int degree, const std::vector<glm::vec<3, T>> &throughPoints,
                               int controlPointsCount, tinynurbs::RationalCurve<T> &curve)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    if (!(degree >= 0 && degree <= NURBSMaxDegree) || controlPointsCount <= 0)
    {
        return false;
    }

    int n             = controlPointsCount;
    int m             = throughPoints.size();
    std::vector<T> uk = getChordParameterization(throughPoints);
    std::vector<T> knotVector(n + degree + 1, 1.0);

    T d = (T) m / (T) n;
    for (int j = 0; j <= degree; j++)
    {
        knotVector[j] = 0.0;
    }
    for (int j = 1; j < n - degree; j++)
    {
        knotVector[degree + j] = 0.0;
        for (int k = j; k < j + degree; k++)
        {
            int i1 = k * d;
            T a    = k * d - i1;
            int i2 = ((k - 1) * d);
            knotVector[degree + j] += a * uk[i2] + (1 - a) * uk[i1];
        }
        knotVector[degree + j] /= degree;
    }

    std::vector<tvecn> controlPoints(n);
    std::vector<T> weights(n);
    std::vector<tvecn> R(n);
    std::vector<tvecn> rk(m);
    std::vector<T> funs(degree + 1);
    std::vector<std::vector<T>> N(m, std::vector<T>(n));
    R[0]            = throughPoints[0];
    R[n - 1]        = throughPoints[m - 1];
    N[0][0]         = 1.0;
    N[m - 1][n - 1] = 1.0;

    for (int i = 0; i < m; i++)
    {
        int spanIndex = findSpan(degree, knotVector, uk[i]);
        std::vector<T> basis;
        basis = basisFuns(degree, spanIndex, knotVector, uk[i]);
        for (int j = 0; j <= degree; j++)
        {
            N[i][spanIndex - degree + j] = basis[j];
        }
        rk[i] = throughPoints[i] - N[i][0] * throughPoints[0] - N[i][n - 1] * throughPoints[m - 1];
    }

    for (int i = 0; i < n; i++)
    {
        R[i] = tvecn(0.0, 0.0, 0.0);
        ;
        for (int j = 0; j < m; j++)
        {
            R[i] += N[j][i] * rk[j];
        }

        if (IsAlmostEqualTo(R[i], tvecn(0.0, 0.0, 0.0)))
        {
            return false;
        }
    }

    if (n - 2 > 0)
    {
        std::vector<std::vector<T>> X(n - 2, std::vector<T>(3));
        std::vector<std::vector<T>> B(n - 2, std::vector<T>(3));
        std::vector<std::vector<T>> Ns(m - 2, std::vector<T>(n - 2));
        for (int i = 0; i < n - 2; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                B[i][j] = R[i + 1][j];
            }
        }
        for (int i = 1; i <= m - 2; i++)
        {
            for (int j = 1; j <= n - 2; j++)
            {
                Ns[i - 1][j - 1] = N[i][j];
            }
        }
        std::vector<std::vector<T>> NsT;
        Transpose(Ns, NsT);
        auto NsTNs = MatrixMultiply(NsT, Ns);
        X          = solveLinearSystem(NsTNs, B);

        for (int i = 0; i < n - 2; i++)
        {
            T x                  = X[i][0];
            T y                  = X[i][1];
            T z                  = X[i][2];
            controlPoints[i + 1] = tvecn(x, y, z);
            weights[i + 1]       = (T) 1;
        }
    }
    controlPoints[0]     = throughPoints[0];
    weights[0]           = 1;
    controlPoints[n - 1] = throughPoints[m - 1];
    weights[n - 1]       = 1;

    curve.degree         = degree;
    curve.knots          = knotVector;
    curve.control_points = controlPoints;
    curve.weights        = weights;
    return true;
}

// approximate method
template<typename T>
bool globalApproximation(const std::vector<std::vector<glm::vec<3, T>>> &throughPoints, int degreeU, int degreeV,
                         int controlPointsRows, int controlPointsColumns, tinynurbs::RationalSurface<T> &surface)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    if (throughPoints.empty() || throughPoints[0].empty() || degreeU < 1 || degreeV < 1)
    {
        return false;
    }

    int rows    = controlPointsRows;
    int n       = rows - 1;
    int columns = controlPointsColumns;
    int m       = columns - 1;

    std::vector<T> knotVectorU;
    std::vector<T> knotVectorV;
    // std::vector<std::vector<tvecn> > controlPoints;
    std::vector<tvecn> controlPoints1D;

    std::vector<std::vector<tvecn>> tempControlPoints;
    for (int i = 0; i < rows; i++)
    {
        tinynurbs::RationalCurve<T> tc;
        bool result = LeastSquaresApproximation(degreeU, throughPoints[i], rows, tc);
        if (!result) return false;
        std::vector<tvecn> points = tc.control_points;
        tempControlPoints.emplace_back(points);
        knotVectorU = tc.knots;
    }

    std::vector<std::vector<tvecn>> preControlPoints;
    std::vector<std::vector<tvecn>> tPoints;
    for (int i = 0; i < columns; i++)
    {
        std::vector<tvecn> c = GetColumn(tempControlPoints, i);
        tinynurbs::RationalCurve<T> tc;
        bool result = LeastSquaresApproximation(degreeV, c, columns, tc);
        if (!result) return false;
        std::vector<tvecn> points = tc.control_points;
        tPoints.emplace_back(points);
        knotVectorV = tc.knots;
    }
    Transpose(tPoints, preControlPoints);
    // controlPoints = preControlPoints;

    for (const auto &row: preControlPoints)
        for (const auto &point: row)
            controlPoints1D.push_back(point);

    std::vector<T> weights(controlPoints1D.size(), (T) 1.0);

    surface.degree_u       = degreeU;
    surface.degree_v       = degreeV;
    surface.knots_u        = knotVectorU;
    surface.knots_v        = knotVectorV;
    surface.control_points = tinynurbs::array2<tvecn>(preControlPoints.size(), preControlPoints[0].size(),
                                                      controlPoints1D);
    surface.weights        = tinynurbs::array2<T>(preControlPoints.size(), preControlPoints[0].size(), weights);
    return true;
}

template<typename T>
void createSweepSurfaceWithApproximate(const tinynurbs::RationalCurve<T> &profile,
                                       const tinynurbs::RationalCurve<T> &trajectory,
                                       tinynurbs::RationalSurface<T> &surface = {})
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    // Determine number of sections for lofting.
    const int q                         = trajectory.degree;
    std::vector<T> trajectoryKnotVector = trajectory.knots;
    const int ktv                       = trajectoryKnotVector.size();
    // const int K = minimumProfiles;
    // int nsect = K + 1;

    // Compute trajectory normals.
    std::vector<tvecn> Bv = projectNormal(trajectory);

    // const std::vector<tvecnp1> &profileControlPoints = profile.ControlPoints;
    std::vector<tvecnp1> profileControlPoints;
    profileControlPoints.reserve(profile.control_points.size());
    for (size_t i = 0; i < profile.control_points.size(); i++)
    {
        profileControlPoints.push_back(
                tvecnp1(tinynurbs::util::cartesianToHomogenous(profile.control_points[i], profile.weights[i])));
    }

    int profileCpSize = profileControlPoints.size();

    // for each lofting section
    std::vector<tinynurbs::RationalCurve<T>> sections(ktv);
    std::vector<tvecnp1> transformedControlPoints(profileCpSize);
    for (int k = 0; k < ktv; k++)
    {
        // Build local coordinate system for the section.
        tvecn zdir             = glm::normalize(rationalCurveDerivatives(trajectory, 1, trajectoryKnotVector[k])[1]);
        unsigned int spanIndex = findSpan(trajectory.degree, trajectory.knots, trajectoryKnotVector[k]);
        const tvecn &xdir      = Bv[spanIndex];
        // tvecn ydir = zdir.CrossProduct(xdir).Normalize();
        tvecn ydir = glm::normalize(glm::cross(zdir, xdir));

        ///
        // frame.push_back({xdir, ydir, zdir});

        // Transform the input profile to the lofting position.

        if (k == 0)
            for (int i = 0; i < profileCpSize; i++)
            {
                tvecnp1 wp = profileControlPoints[i];
                T w        = wp.w;
                // tvecn p = wp.ToXYZ(true);
                tvecn p = tinynurbs::util::homogenousToCartesian(wp);

                // tvecn transformed = worldToLocal(p, xdir, ydir, zdir, rationalCurvePoint(trajectoryCopy, vk[k]));
                tvecn transformed           = worldToLocal(rationalCurvePoint(trajectory, trajectoryKnotVector[k]), xdir, ydir,
                                                           zdir, p);
                transformedControlPoints[i] = tvecnp1(transformed, w);
            }

        std::vector<tvecn> transformedPoints;
        std::vector<T> weights;
        tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformedPoints, weights);

        // Make section curve.
        // note: std::vector::swap relinks pointers, which is faster than copying.
        tinynurbs::RationalCurve<T> &section = sections[k];
        section.degree                       = profile.degree;
        section.knots                        = profile.knots;
        section.weights.swap(weights);
        section.control_points.swap(transformedPoints);

        for (int j = 0; j < profileCpSize; j++)
        {
            section.control_points[j] =
                    LocalToWorld(rationalCurvePoint(trajectory, trajectoryKnotVector[k]), xdir, ydir, zdir,
                                 section.control_points[j]);
        }
    }

    // v_bar = vk;
    // sect = sections[2];

    std::vector<std::vector<glm::vec<3, T>>> throughPoints;

    for (auto &section: sections)
    {
        std::vector<glm::vec<3, T>> points;
        for (auto &point: section.control_points)
        {
            points.emplace_back(point);
        }
        throughPoints.emplace_back(points);
    }

    // std::cout << "throughPoints.rows = " << throughPoints.size() << std::endl;
    // std::cout << "throughPoints.colums = " << throughPoints[0].size() << std::endl;

    // approximate method
    globalApproximation(throughPoints, 3, 3, 7, 7, surface);
}

template<typename T>
std::vector<glm::vec<3, T>> ComputeTangent(const std::vector<glm::vec<3, T>> &throughPoints)
{
    typedef glm::vec<3, T> tvecn;

    int size = throughPoints.size();
    std::vector<tvecn> tangents(size, tvecn(0, 0, 0));
    std::vector<tvecn> qq(size, tvecn(0, 0, 0));
    std::vector<T> delta(size, 0.0);

    auto params = getChordParameterization(throughPoints);
    for (int i = 1; i < size; i++)
    {
        delta[i] = params[i] - params[i - 1];
        qq[i]    = throughPoints[i] - throughPoints[i - 1];
    }
    for (int i = 1; i < size - 1; i++)
    {
        T a         = delta[i] / (delta[i] + delta[i + 1]);
        tangents[i] = glm::normalize((T) (1 - a) * qq[i] + (T) a * qq[i + 1]);
    }

    tangents[0]                   = glm::normalize((T) 2 * qq[1] / delta[1] - tangents[1]);
    tangents[tangents.size() - 1] = glm::normalize((T) 2 * qq[qq.size() - 1] / delta[delta.size() - 1] - tangents[tangents.size() - 2]);
    return tangents;
}

template<typename T>
std::vector<T> GetIndex(T size)
{
    std::vector<T> ind(2 * (size - 1) + 2);
    ind[0]              = 0;
    ind[ind.size() - 1] = 3 * size - 3;
    T ii                = 1;
    T jj                = 1;
    for (int i = 0; i < size - 1; i++)
    {
        ind[ii]     = jj;
        ind[ii + 1] = jj + 1;
        ii          = ii + 2;
        jj          = jj + 3;
    }
    return ind;
}

template<typename T>
bool BicubicLocalInterpolation(const std::vector<std::vector<glm::vec<3, T>>> &throughPoints,
                               tinynurbs::RationalSurface<T> &surface)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    if (throughPoints.size() <= 0 || throughPoints[0].size() <= 0)
        return false;

    int degreeU = 3;
    int degreeV = 3;

    std::vector<T> knotVectorU;
    std::vector<T> knotVectorV;

    int row    = throughPoints.size();
    int n      = row - 1;
    int column = throughPoints[0].size();
    int m      = column - 1;

    std::vector<std::vector<std::vector<tvecn>>> td(
            n + 1, std::vector<std::vector<tvecn>>(m + 1, std::vector<tvecn>(3, tvecn(0.0, 0.0, 0.0))));

    std::vector<T> ub(n + 1, 0.0);
    std::vector<T> vb(m + 1, 0.0);

    std::vector<T> r(m + 1);
    std::vector<T> s(n + 1);

    T total = 0.0;
    for (int l = 0; l <= m; l++)
    {
        std::vector<tvecn> columnData = GetColumn(throughPoints, l);
        std::vector<tvecn> tvkl       = ComputeTangent(columnData);

        r[l] = 0.0;
        for (int k = 0; k <= n; k++)
        {
            td[k][l][1] = tvkl[k];

            if (k > 0)
            {
                T d = glm::distance(throughPoints[k][l], throughPoints[k - 1][l]);
                ub[k] += d;
                r[l] += d;
            }
        }
        total += r[l];
    }
    for (int k = 1; k < n; k++)
    {
        ub[k] = ub[k - 1] + ub[k] / total;
    }
    ub[n] = 1.0;
    total = 0.0;

    for (int k = 0; k <= n; k++)
    {
        std::vector<tvecn> tukl = ComputeTangent(throughPoints[k]);
        s[k]                    = 0.0;
        for (int l = 0; l <= m; l++)
        {
            td[k][l][0] = tukl[l];

            if (l > 0)
            {
                T d = glm::distance(throughPoints[k][l], throughPoints[k][l - 1]);
                vb[l] += d;
                s[k] += d;
            }
        }
        total += s[k];
    }
    for (int l = 1; l < m; l++)
    {
        vb[l] = vb[l - 1] + vb[l] / total;
    }
    vb[m] = 1.0;
    total = 0.0;

    int kuSize = 2 * ub.size() + 2 + 2;
    knotVectorU.resize(kuSize);
    for (int i = 0; i < 4; i++)
    {
        knotVectorU[i] = knotVectorU[i + 1] = knotVectorU[i + 2] = knotVectorU[i + 3] = 0.0;
        knotVectorU[kuSize - 1] = knotVectorU[kuSize - 2] = knotVectorU[kuSize - 3] = knotVectorU[kuSize - 4] = 1.0;
    }
    int ii = 4;
    for (int i = 1; i < ub.size() - 1; i++)
    {
        knotVectorU[ii]     = ub[i];
        knotVectorU[ii + 1] = ub[i];
        ii += 2;
    }

    int kvSize = 2 * vb.size() + 2 + 2;
    knotVectorV.resize(kvSize);
    for (int i = 0; i < 4; i++)
    {
        knotVectorV[i] = knotVectorV[i + 1] = knotVectorV[i + 2] = knotVectorV[i + 3] = 0.0;
        knotVectorV[kvSize - 1] = knotVectorV[kvSize - 2] = knotVectorV[kvSize - 3] = knotVectorV[kvSize - 4] = 1.0;
    }

    ii = 4;
    for (int i = 1; i < vb.size() - 1; i++)
    {
        knotVectorV[ii]     = vb[i];
        knotVectorV[ii + 1] = vb[i];
        ii += 2;
    }

    std::vector<std::vector<tvecn>> bezierControlPoints(3 * row - 2,
                                                        std::vector<tvecn>(3 * column - 2, tvecn(0, 0, 0)));
    for (int i = 0; i < row; i++)
    {
        int n = throughPoints[i].size();
        std::vector<tvecn> Tmp;
        for (int c = 0; c < td[0].size(); c++)
        {
            Tmp.emplace_back(td[i][c][1]);
        }

        std::vector<tvecn> temp(3 * n - 2, tvecn(0, 0, 0));
        for (int j = 0; j < n; j++)
        {
            temp[3 * i] = throughPoints[i][j];
        }
        for (int j = 0; j < n - 1; j++)
        {
            T a             = (ub[j + 1] - ub[j]) * getTotalChordLength(throughPoints[i]);
            temp[3 * j + 1] = throughPoints[i][j] + a / (T) 3.0 * Tmp[j];
            temp[3 * j + 2] = throughPoints[i][j + 1] - a / (T) 3.0 * Tmp[j];
        }
        bezierControlPoints[3 * i] = temp;
    }
    for (int i = 0; i < column; i++)
    {
        auto columnData = GetColumn(throughPoints, i);
        int n           = columnData.size();
        std::vector<tvecn> Tmp;
        for (int r = 0; r < td.size(); r++)
        {
            Tmp.emplace_back(td[r][i][0]);
        }

        std::vector<tvecn> temp(3 * n - 2, tvecn(0, 0, 0));
        for (int j = 0; j < n; j++)
        {
            temp[3 * i] = columnData[j];
        }
        for (int j = 0; j < n - 1; j++)
        {
            T a             = (vb[j + 1] - vb[j]) * getTotalChordLength(columnData);
            temp[3 * j + 1] = columnData[j] + a / (T) 3.0 * Tmp[j];
            temp[3 * j + 2] = columnData[j + 1] - a / (T) 3.0 * Tmp[j];
        }
        for (int r = 0; r < bezierControlPoints.size(); r++)
        {
            bezierControlPoints[r][3 * i] = temp[r];
        }
    }


    for (int k = 1; k < n; k++)
    {
        T ak = (ub[k] - ub[k - 1]) / ((ub[k] - ub[k - 1]) + (ub[k + 1] - ub[k]));
        for (int l = 1; l < m; l++)
        {
            T bl = bl = (vb[l] - vb[l - 1]) / ((vb[l] - vb[l - 1]) + (vb[l + 1] - vb[l]));

            tvecn dvukl = (T) (1 - ak) * (td[k][l][1] - td[k - 1][l][1]) / (ub[k] - ub[k - 1]) + ak * (td[k + 1][l][1] - td[k][l][1]) / (ub[k + 1] - ub[k]);
            tvecn duvkl = (T) (1 - bl) * (td[k][l][0] - td[k][l - 1][0]) / (vb[l] - vb[l - 1]) + bl * (td[k][l + 1][0] - td[k][l][0]) / (vb[l + 1] - vb[l]);

            td[k][l][2] = (ak * duvkl + bl * dvukl) / (ak + bl);
        }
    }

    for (int k = 0; k < n; k++)
    {
        for (int l = 0; l < m; l++)
        {
            T gamma = (ub[k + 1] - ub[k]) * (vb[l + 1] - vb[l]) / 9.0;
            int ii  = 3 * k;
            int jj  = 3 * l;
            bezierControlPoints[ii + 1][jj + 1] =
                    gamma * td[k][l][2] + bezierControlPoints[ii][jj + 1] + bezierControlPoints[ii + 1][jj] -
                    bezierControlPoints[ii][jj];
            bezierControlPoints[ii + 2][jj + 1] =
                    -gamma * td[k + 1][l][2] + bezierControlPoints[ii + 3][jj + 1] - bezierControlPoints[ii + 3][jj] +
                    bezierControlPoints[ii + 2][jj];
            bezierControlPoints[ii + 1][jj + 2] =
                    -gamma * td[k][l + 1][2] + bezierControlPoints[ii + 1][jj + 3] - bezierControlPoints[ii][jj + 3] +
                    bezierControlPoints[ii][jj + 2];
            bezierControlPoints[ii + 2][jj + 2] =
                    gamma * td[k + 1][l + 1][2] + bezierControlPoints[ii + 2][jj + 3] + bezierControlPoints[ii + 3][jj + 2] - bezierControlPoints[ii + 3][jj + 3];
        }
    }

    std::vector<std::vector<tvecn>> columnFilter;
    auto ind = GetIndex(column);
    for (int c = 0; c < ind.size(); c++)
    {
        auto columnData = GetColumn(bezierControlPoints, ind[c]);
        columnFilter.emplace_back(columnData);
    }
    std::vector<std::vector<tvecn>> Tcf;
    Transpose(columnFilter, Tcf);
    ind = GetIndex(row);
    std::vector<std::vector<tvecn>> rowFilter;
    for (int r = 0; r < ind.size(); r++)
    {
        rowFilter.emplace_back(Tcf[ind[r]]);
    }

    std::vector<tvecn> controlPoints1D;

    for (const auto &row: rowFilter)
        for (const auto &point: row)
            controlPoints1D.push_back(point);

    std::vector<T> weights(controlPoints1D.size(), (T) 1.0);

    surface.degree_u       = degreeU;
    surface.degree_v       = degreeV;
    surface.knots_u        = knotVectorU;
    surface.knots_v        = knotVectorV;
    surface.control_points = tinynurbs::array2<tvecn>(rowFilter.size(), rowFilter[0].size(), controlPoints1D);
    surface.weights        = tinynurbs::array2<T>(rowFilter.size(), rowFilter[0].size(), weights);
    return true;
}

// Frenet Frame
// template<typename T>
// void computeFrenetFrame(const tinynurbs::RationalCurve<T> &curve, T t,
//                         glm::vec<3, T> &tangent,
//                         glm::vec<3, T> &normal,
//                         glm::vec<3, T> &binormal)
// {
//     typedef glm::vec<3, T> tvecn;
//
//
//     auto derivatives = rationalCurveDerivatives(curve, 2, t);
//     tangent = glm::normalize(derivatives[1]);
//
//
//     normal = glm::normalize(derivatives[2] - tangent * glm::dot(derivatives[2], tangent));
//
//
//     binormal = glm::normalize(glm::cross(tangent, normal));
// }

template<typename T>
glm::mat4 maxtrixT1(T t)
{
    T angle = glm::pi<T>() * t;
    return glm::mat4(
            glm::vec4(cos(angle), -sin(angle), 0.0f, 0.0f),
            glm::vec4(sin(angle), cos(angle), 0.0f, 0.0f),
            glm::vec4(0.0f, 0.0f, 1.0f, 0.0f),
            glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
}

template<typename T>
glm::mat4 maxtrixT2(T t)
{
    T angle = glm::pi<T>() * t;
    return glm::mat4(
            glm::vec4(2 * cos(angle), -sin(angle), 0.0f, 0.0f),
            glm::vec4(sin(angle), 3 * cos(angle), 0.0f, 0.0f),
            glm::vec4(0.0f, 0.0f, 4.0f, 0.0f),
            glm::vec4(0.0f, 0.0f, 0.0f, 1.0f));
}

template<typename T>
tinynurbs::RationalCurve<T> offsetProfile3D(const tinynurbs::RationalCurve<T> &profile,
                                            const glm::vec<3, T> &origin,
                                            const glm::vec<3, T> &tangent,
                                            const glm::vec<3, T> &normal,
                                            const glm::vec<3, T> &binormal)
{
    typedef glm::vec<3, T> tvecn;

    std::vector<tvecn> offset_points;
    for (const auto &point: profile.control_points)
    {
        tvecn transformed = origin + point.x * tangent + point.y * normal + point.z * binormal;
        offset_points.push_back(transformed);
    }


    tinynurbs::RationalCurve<T> offset_curve = profile;
    offset_curve.control_points              = offset_points;
    return offset_curve;
}

template<typename T>
tinynurbs::RationalCurve<T> offsetProfile3DWithShape(const tinynurbs::RationalCurve<T> &profile,
                                                     const glm::vec<3, T> &origin,
                                                     const glm::vec<3, T> &tangent,
                                                     const glm::vec<3, T> &normal,
                                                     const glm::vec<3, T> &binormal,
                                                     const glm::mat4 &shapeMatrix)
{
    tinynurbs::RationalCurve<T> offset_curve = profile;
    for (size_t i = 0; i < profile.control_points.size(); ++i)
    {
        glm::vec<3, T> p               = profile.control_points[i];
        glm::vec4 transformed          = shapeMatrix * glm::vec4(p, profile.weights[i]);
        p                              = tinynurbs::util::homogenousToCartesian(transformed);
        offset_curve.control_points[i] = origin +
                                         p.x * tangent +
                                         p.y * normal +
                                         p.z * binormal;
    }
    return offset_curve;
}

template<typename T>
void createSweepSurfaceWithInterpolation(const tinynurbs::RationalCurve<T> &profile,
                                         const tinynurbs::RationalCurve<T> &trajectory, int num_samples = 0,
                                         tinynurbs::RationalSurface<T> &surface = {}, std::vector<T> &v_bar = {},
                                         std::vector<std::vector<glm::vec<3, T>>> &frame = {}, int frameType = 0,
                                         int matrixType = 0)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    const int q                         = trajectory.degree;
    std::vector<T> trajectoryKnotVector = trajectory.knots;
    const int ktv                       = trajectoryKnotVector.size();
    const int K                         = num_samples;
    int nsect                           = K + 1;

    tinynurbs::RationalCurve<T> trajectoryCopy = trajectory;
    if (ktv <= nsect + q)
    {
        // Insert knots to widest spans
        // to reach the required minimum number of sections.
        int m                           = nsect + q - ktv + 1;
        std::vector<T> insertedElements = getInsertedKnotElements(m, trajectoryKnotVector);
        refineKnotVector(trajectory, insertedElements, trajectoryCopy);
    }
    else
    {
        if (ktv > nsect + q + 1)
        {
            nsect = ktv - q - 1;
        }
    }

    // Determine lofting profile positions.
    std::vector<T> knotVectorV = trajectoryCopy.knots;
    auto minmax                = std::minmax_element(knotVectorV.begin(), knotVectorV.end());

    std::vector<T> vk(nsect);
    vk[0]         = *minmax.first;
    vk[nsect - 1] = *minmax.second;

    for (int k = 1; k < nsect - 1; k++)
    {
        T sum = 0.0;
        for (int i = 1; i <= q; i++)
        {
            sum += knotVectorV[k + i];
        }
        vk[k] = sum / (T) q;
    }

    // path.knots = segments;

    v_bar = vk;

    // std::vector<tvecn> Bv = projectNormal(trajectoryCopy);
    // const std::vector<tvecn> &profileControlPoints = profile.ControlPoints;
    // const std::vector<T> &profileWeights = profile.weights;

    std::vector<tvecnp1> profileControlPoints = tinynurbs::util::cartesianToHomogenous(profile.control_points,
                                                                                       profile.weights);
    int profileCpSize                         = profileControlPoints.size();

    std::vector<tinynurbs::RationalCurve<T>> sections(nsect);
    std::vector<tvecnp1> transformedControlPoints(profileCpSize);

    std::vector<tvecn> tangents, normals, binormals;

    switch (frameType)
    {
        case 0:
            computeFrenetFrames(trajectoryCopy, nsect, tangents, normals, binormals);
            break;
        case 1:
            computeProjectionNormalFrame(trajectoryCopy, nsect, tangents, normals, binormals);
            break;
        case 2:
            computeRMF(trajectoryCopy, nsect, tangents, normals, binormals);
        default:
            break;
    }

    // frame
    for (int i = 0; i < tangents.size(); i++)
    {
        frame.push_back({tangents[i], normals[i], binormals[i]});
    }

    for (int k = 0; k < nsect; k++)
    {
        // tvecn zdir = glm::normalize(rationalCurveDerivatives(path, 1, segments[k])[1]);
        // int spanIndex = findSpan(path.degree, path.knots, segments[k]);
        // tvecn xdir = Bv[spanIndex];
        // tvecn ydir = glm::normalize(glm::cross(zdir, xdir));

        //frame
        // frame.push_back({xdir, ydir, zdir});

        if (k == 0)
            for (int i = 0; i < profileCpSize; i++)
            {
                tvecnp1 wp = profileControlPoints[i];
                T w        = wp.w;
                // tvecn p = wp.ToXYZ(true);
                tvecn p = tinynurbs::util::homogenousToCartesian(wp);

                // tvecn transformed = worldToLocal(p, xdir, ydir, zdir, rationalCurvePoint(trajectoryCopy, vk[k]));
                tvecn transformed           = worldToLocal(rationalCurvePoint(trajectoryCopy, vk[k]), tangents[k], normals[k],
                                                           binormals[k], p);
                transformedControlPoints[i] = tvecnp1(transformed, w);
            }

        std::vector<tvecn> transformedPoints;
        std::vector<T> weights;
        tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformedPoints, weights);

        // Make section curve.
        // note: std::vector::swap relinks pointers, which is faster than copying.
        tinynurbs::RationalCurve<T> &section = sections[k];
        section.degree                       = profile.degree;
        section.knots                        = profile.knots;
        section.weights.swap(weights);
        section.control_points.swap(transformedPoints);

        tvecn path_point = rationalCurvePoint(trajectoryCopy, vk[k]);

        // for (int j = 0; j < profileCpSize; j++)
        // {
        //     section.control_points[j] =
        //             LocalToWorld(rationalCurvePoint(trajectoryCopy, vk[k]), tangents[k], normals[k],
        //                          binormals[k], section.control_points[j]);

        if (matrixType == 0)
        {
            section = offsetProfile3D(
                    section, path_point, tangents[k], normals[k],
                    binormals[k]);
        }
        else if (matrixType == 1)
        {
            section = offsetProfile3DWithShape(
                    section, path_point, tangents[k], normals[k],
                    binormals[k], maxtrixT1(vk[k]));
        }
        else
        {
            section = offsetProfile3DWithShape(
                    section, path_point, tangents[k], normals[k],
                    binormals[k], maxtrixT2(vk[k]));
        }

        // tvecn p = section.control_points[j];
        // glm::vec4 transformed = maxtrixT1(vk[k]) * glm::vec4(p, section.weights[j]);
        // section.control_points[j] = tinynurbs::util::homogenousToCartesian(transformed);
        // }
    }

    // v_bar = vk;
    // sect = sections[2];

    std::vector<std::vector<tvecn>> throughPoints(sections.size());
    for (int i = 0; i < sections.size(); i++)
    {
        for (int j = sections[i].degree; j < sections[i].knots.size() - 1 - sections[i].degree; ++j)
        {
            throughPoints[i].emplace_back(rationalCurvePoint(sections[i], sections[i].knots[j]));
        }
    }


    // for (auto &section: sections)
    // {
    //     std::vector<glm::vec<3, T> > points;
    //     for (auto &point: section.control_points)
    //     {
    //         points.emplace_back(point);
    //     }
    //     throughPoints.emplace_back(points);
    // }

    transpose(throughPoints);

    // std::cout << "throughPoints.rows = " << throughPoints.size() << std::endl;
    // std::cout << "throughPoints.colums = " << throughPoints[0].size() << std::endl;

    // approximate method
    BicubicLocalInterpolation(throughPoints, surface);
}

template<typename T>
void createSweepSurfaceWithSkinning(const tinynurbs::RationalCurve<T> &profile,
                                    const tinynurbs::RationalCurve<T> &trajectory,
                                    int num_samples                        = 0,
                                    tinynurbs::RationalSurface<T> &surface = {}, std::vector<T> &v_bar = {},
                                    std::vector<std::vector<glm::vec<3, T>>> &frame = {}, int frameType = 0,
                                    int matrixType = 0)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    const int q                         = trajectory.degree;
    std::vector<T> trajectoryKnotVector = trajectory.knots;
    const int ktv                       = trajectoryKnotVector.size();
    const int K                         = num_samples;
    int nsect                           = K + 1;

    tinynurbs::RationalCurve<T> trajectoryCopy = trajectory;
    if (ktv <= nsect + q)
    {
        // Insert knots to widest spans
        // to reach the required minimum number of sections.
        int m                           = nsect + q - ktv + 1;
        std::vector<T> insertedElements = getInsertedKnotElements(m, trajectoryKnotVector);
        refineKnotVector(trajectory, insertedElements, trajectoryCopy);
    }
    else
    {
        if (ktv > nsect + q + 1)
        {
            nsect = ktv - q - 1;
        }
    }

    // Determine lofting profile positions.
    std::vector<T> knotVectorV = trajectoryCopy.knots;
    auto minmax                = std::minmax_element(knotVectorV.begin(), knotVectorV.end());

    std::vector<T> vk(nsect);
    vk[0]         = *minmax.first;
    vk[nsect - 1] = *minmax.second;

    for (int k = 1; k < nsect - 1; k++)
    {
        T sum = 0.0;
        for (int i = 1; i <= q; i++)
        {
            sum += knotVectorV[k + i];
        }
        vk[k] = sum / (T) q;
    }

    // path.knots = segments;

    v_bar = vk;

    // std::vector<tvecn> Bv = projectNormal(trajectoryCopy);
    // const std::vector<tvecn> &profileControlPoints = profile.ControlPoints;
    // const std::vector<T> &profileWeights = profile.weights;

    std::vector<tvecnp1> profileControlPoints = tinynurbs::util::cartesianToHomogenous(profile.control_points,
                                                                                       profile.weights);
    int profileCpSize                         = profileControlPoints.size();

    std::vector<tinynurbs::RationalCurve<T>> sections(nsect);
    std::vector<tvecnp1> transformedControlPoints(profileCpSize);

    std::vector<tvecn> tangents, normals, binormals;

    switch (frameType)
    {
        case 0:
            computeFrenetFrames(trajectoryCopy, nsect, tangents, normals, binormals);
            break;
        case 1:
            computeProjectionNormalFrame(trajectoryCopy, nsect, tangents, normals, binormals);
            break;
        case 2:
            computeRMF(trajectoryCopy, nsect, tangents, normals, binormals);
        default:
            break;
    }

    // frame
    for (int i = 0; i < tangents.size(); i++)
    {
        frame.push_back({tangents[i], normals[i], binormals[i]});
    }

    for (int k = 0; k < nsect; k++)
    {
        // tvecn zdir = glm::normalize(rationalCurveDerivatives(path, 1, segments[k])[1]);
        // int spanIndex = findSpan(path.degree, path.knots, segments[k]);
        // tvecn xdir = Bv[spanIndex];
        // tvecn ydir = glm::normalize(glm::cross(zdir, xdir));

        //frame
        // frame.push_back({xdir, ydir, zdir});

        if (k == 0)
            for (int i = 0; i < profileCpSize; i++)
            {
                tvecnp1 wp = profileControlPoints[i];
                T w        = wp.w;
                // tvecn p = wp.ToXYZ(true);
                tvecn p = tinynurbs::util::homogenousToCartesian(wp);

                // tvecn transformed = worldToLocal(p, xdir, ydir, zdir, rationalCurvePoint(trajectoryCopy, vk[k]));
                tvecn transformed           = worldToLocal(rationalCurvePoint(trajectoryCopy, vk[k]), tangents[k], normals[k],
                                                           binormals[k], p);
                transformedControlPoints[i] = tvecnp1(transformed, w);
            }

        std::vector<tvecn> transformedPoints;
        std::vector<T> weights;
        tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformedPoints, weights);

        // Make section curve.
        // note: std::vector::swap relinks pointers, which is faster than copying.
        tinynurbs::RationalCurve<T> &section = sections[k];
        section.degree                       = profile.degree;
        section.knots                        = profile.knots;
        section.weights.swap(weights);
        section.control_points.swap(transformedPoints);

        tvecn path_point = rationalCurvePoint(trajectoryCopy, vk[k]);

        // for (int j = 0; j < profileCpSize; j++)
        // {
        //     section.control_points[j] =
        //             LocalToWorld(rationalCurvePoint(trajectoryCopy, vk[k]), tangents[k], normals[k],
        //                          binormals[k], section.control_points[j]);

        if (matrixType == 0)
        {
            section = offsetProfile3D(
                    section, path_point, tangents[k], normals[k],
                    binormals[k]);
        }
        else if (matrixType == 1)
        {
            section = offsetProfile3DWithShape(
                    section, path_point, tangents[k], normals[k],
                    binormals[k], maxtrixT1(vk[k]));
        }
        else
        {
            section = offsetProfile3DWithShape(
                    section, path_point, tangents[k], normals[k],
                    binormals[k], maxtrixT2(vk[k]));
        }

        // tvecn p = section.control_points[j];
        // glm::vec4 transformed = maxtrixT1(vk[k]) * glm::vec4(p, section.weights[j]);
        // section.control_points[j] = tinynurbs::util::homogenousToCartesian(transformed);
        // }
    }

    CreateLoftSurface(sections, surface);
}

template<typename T>
void createSweepSurfaceWithOffset(const tinynurbs::RationalCurve<T> &profile,
                                  const tinynurbs::RationalCurve<T> &trajectory,
                                  T offset_distance, int num_samples,
                                  tinynurbs::RationalSurface<T> &surface, std::vector<T> &v_bar = {},
                                  std::vector<std::vector<glm::vec<3, T>>> &frame = {}, int frameType = 0,
                                  int matrixType = 0)
{
    typedef glm::vec<3, T> tvecn;
    typedef glm::vec<4, T> tvecnp1;

    // Step 1:
    int trajectory_degree           = trajectory.degree;
    std::vector<T> trajectory_knots = trajectory.knots;
    int num_sections                = num_samples + 1;

    // Step 2:
    tinynurbs::RationalCurve<T> refined_trajectory = trajectory;
    if (trajectory_knots.size() <= num_sections + trajectory_degree)
    {
        int additional_knots          = num_sections + trajectory_degree - trajectory_knots.size() + 1;
        std::vector<T> inserted_knots = getInsertedKnotElements(additional_knots, trajectory_knots);
        refineKnotVector(trajectory, inserted_knots, refined_trajectory);
    }

    std::vector<T> knots_v = refined_trajectory.knots;
    auto minmax            = std::minmax_element(knots_v.begin(), knots_v.end());

    // Step 3:
    std::vector<T> vk(num_sections);
    vk[0]                = *minmax.first;
    vk[num_sections - 1] = *minmax.second;
    for (int k = 1; k < num_sections - 1; k++)
    {
        T sum = 0.0;
        for (int i = 1; i <= trajectory_degree; i++)
        {
            sum += knots_v[k + i];
        }
        vk[k] = sum / T(trajectory_degree);
    }

    v_bar = vk;

    // Step 4:
    std::vector<tvecn> tangents, normals, binormals;
    // computeRMF(refined_trajectory, num_sections, tangents, normals, binormals);
    switch (frameType)
    {
        case 0:
            computeFrenetFrames(refined_trajectory, num_sections, tangents, normals, binormals);
            break;
        case 1:
            computeProjectionNormalFrame(refined_trajectory, num_sections, tangents, normals, binormals);
            break;
        case 2:
            computeRMF(refined_trajectory, num_sections, tangents, normals, binormals);
        default:
            break;
    }

    // frame
    for (int i = 0; i < tangents.size(); i++)
    {
        frame.push_back({tangents[i], normals[i], binormals[i]});
    }
    // Step 5:
    std::vector<tinynurbs::RationalCurve<T>> offset_sections(num_sections);
    std::vector<tvecnp1> transformedControlPoints(profile.control_points.size());
    std::vector<tvecnp1> transformedControlPoints3D(profile.control_points.size());
    // offset_sections[0] = profile;
    for (int k = 0; k < num_sections; k++)
    {
        std::vector<tvecn> transformed_control_points;
        std::vector<T> weights = profile.weights;

        if (k == 0)
        {
            for (int i = 0; i < profile.control_points.size(); i++)
            {
                tvecn p = profile.control_points[i];
                T w     = profile.weights[i];
                // tvecn p = wp.ToXYZ(true);
                // tvecn p = tinynurbs::util::homogenousToCartesian(wp);

                // tvecn transformed = worldToLocal(p, xdir, ydir, zdir, rationalCurvePoint(trajectoryCopy, vk[k]));
                tvecn transformed           = worldToLocal(rationalCurvePoint(refined_trajectory, vk[k]), tangents[k], normals[k],
                                                           binormals[k], p);
                transformedControlPoints[i] = tvecnp1(transformed, w);
            }

            tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformed_control_points, weights);
            tinynurbs::RationalCurve<T> section;
            section.degree         = profile.degree;
            section.knots          = profile.knots;
            section.control_points = transformed_control_points;
            section.weights        = weights;

            offset_sections[k] = section;

            continue;
        }

        tinynurbs::util::homogenousToCartesian(transformedControlPoints, transformed_control_points, weights);

        // rationalCurveDerivatives(refined_trajectory, 2, vk[k])[1];
        //
        for (auto &point: transformed_control_points)
        {
            // tvecn offset_point = point;
            if (k < num_sections - 1)
                point = point + offset_distance * normals[k - 1];
            // tvecn world_point = LocalToWorld(rationalCurvePoint(refined_trajectory, vk[k]),
            //                                  tangents[k], normals[k], binormals[k], offset_point);
            // point = world_point;
        }

        //
        tinynurbs::RationalCurve<T> section;
        section.degree         = profile.degree;
        section.knots          = profile.knots;
        section.control_points = transformed_control_points;
        section.weights        = weights;

        offset_sections[k] = section;
    }

    for (int k = 0; k < num_sections; k++)
    {
        tvecn path_point = rationalCurvePoint(refined_trajectory, vk[k]);

        if (matrixType == 0)
        {
            offset_sections[k] = offsetProfile3D(
                    offset_sections[k], path_point, tangents[k], normals[k],
                    binormals[k]);
        }
        else if (matrixType == 1)
        {
            offset_sections[k] = offsetProfile3DWithShape(
                    offset_sections[k], path_point, tangents[k], normals[k],
                    binormals[k], maxtrixT1(vk[k]));
        }
        else
        {
            offset_sections[k] = offsetProfile3DWithShape(
                    offset_sections[k], path_point, tangents[k], normals[k],
                    binormals[k], maxtrixT2(vk[k]));
        }
    }

    // Step 6:
    CreateLoftSurface(offset_sections, surface);
}

template<typename T>
bool HasSelfIntersection(const tinynurbs::RationalSurface<T> &surface, T tolerance = T(1e-6))
{
    using Vec3 = glm::vec<3, T>;

    // Number of sampling points along u and v directions
    const int numSamplesU = 50;// Adjust as necessary
    const int numSamplesV = 50;// Adjust as necessary

    // Generate parametric sampling points
    std::vector<T> uSamples(numSamplesU), vSamples(numSamplesV);
    T uStart = surface.knots_u.front();
    T uEnd   = surface.knots_u.back();
    T vStart = surface.knots_v.front();
    T vEnd   = surface.knots_v.back();

    for (int i = 0; i < numSamplesU; ++i)
    {
        uSamples[i] = uStart + (uEnd - uStart) * i / (numSamplesU - 1);
    }
    for (int i = 0; i < numSamplesV; ++i)
    {
        vSamples[i] = vStart + (vEnd - vStart) * i / (numSamplesV - 1);
    }

    // Check for self-intersection
    for (size_t i = 0; i < numSamplesU; ++i)
    {
        for (size_t j = 0; j < numSamplesV; ++j)
        {
            Vec3 point1 = rationalSurfacePoint(surface, uSamples[i], vSamples[j]);

            for (size_t k = i; k < numSamplesU; ++k)
            {
                for (size_t l = (k == i ? j + 1 : 0); l < numSamplesV; ++l)
                {
                    Vec3 point2 = rationalSurfacePoint(surface, uSamples[k], vSamples[l]);

                    // Check distance between the two points
                    if (glm::distance(point1, point2) < tolerance)
                    {
                        return true;
                    }
                }
            }
        }
    }

    return false;// No self-intersections found
}

template<typename T>
struct Triangle
{
    glm::vec<3, T> v0, v1, v2;
};

template<typename T>
bool rayIntersectsTriangle(
        const glm::vec<3, T> &rayOrigin,
        const glm::vec<3, T> &rayVector,
        const Triangle<T> &triangle,
        T epsilon = std::numeric_limits<T>::epsilon())
{
    glm::vec<3, T> edge1 = triangle.v1 - triangle.v0;
    glm::vec<3, T> edge2 = triangle.v2 - triangle.v0;
    glm::vec<3, T> h     = glm::cross(rayVector, edge2);
    T a                  = glm::dot(edge1, h);

    if (a > -epsilon && a < epsilon)
    {
        return false;//
    }

    T f              = 1.0 / a;
    glm::vec<3, T> s = rayOrigin - triangle.v0;
    T u              = f * glm::dot(s, h);

    if (u < 0.0 || u > 1.0)
    {
        return false;
    }

    glm::vec<3, T> q = glm::cross(s, edge1);
    T v              = f * glm::dot(rayVector, q);

    if (v < 0.0 || u + v > 1.0)
    {
        return false;
    }

    //
    T t = f * glm::dot(edge2, q);

    if (t > epsilon)
    {
        //
        // return rayOrigin + rayVector * t;
        return true;
    }
    else
    {
        //
        return false;
    }
}

//
template<typename T>
bool trianglesIntersect(
        const Triangle<T> &tri1,
        const Triangle<T> &tri2,
        T epsilon = std::numeric_limits<T>::epsilon())
{
    //
    if (rayIntersectsTriangle(tri1.v0, tri1.v1 - tri1.v0, tri2, epsilon) ||
        rayIntersectsTriangle(tri1.v1, tri1.v2 - tri1.v1, tri2, epsilon) ||
        rayIntersectsTriangle(tri1.v2, tri1.v0 - tri1.v2, tri2, epsilon))
    {
        return true;
    }

    //
    if (rayIntersectsTriangle(tri2.v0, tri2.v1 - tri2.v0, tri1, epsilon) ||
        rayIntersectsTriangle(tri2.v1, tri2.v2 - tri2.v1, tri1, epsilon) ||
        rayIntersectsTriangle(tri2.v2, tri2.v0 - tri2.v2, tri1, epsilon))
    {
        return true;
    }

    return false;
}

//
template<typename T>
std::vector<Triangle<T>> discretizeSurface(const tinynurbs::RationalSurface<T> &surface, int resolution_u,
                                           int resolution_v)
{
    std::vector<Triangle<T>> triangles;
    T u_start = surface.knots_u.front();
    T u_end   = surface.knots_u.back();
    T v_start = surface.knots_v.front();
    T v_end   = surface.knots_v.back();

    for (int i = 0; i < resolution_u; ++i)
    {
        for (int j = 0; j < resolution_v; ++j)
        {
            T u0 = u_start + (u_end - u_start) * i / resolution_u;
            T u1 = u_start + (u_end - u_start) * (i + 1) / resolution_u;
            T v0 = v_start + (v_end - v_start) * j / resolution_v;
            T v1 = v_start + (v_end - v_start) * (j + 1) / resolution_v;

            glm::vec<3, T> p00 = rationalSurfacePoint(surface, u0, v0);
            glm::vec<3, T> p10 = rationalSurfacePoint(surface, u1, v0);
            glm::vec<3, T> p01 = rationalSurfacePoint(surface, u0, v1);
            glm::vec<3, T> p11 = rationalSurfacePoint(surface, u1, v1);

            //
            triangles.push_back({p00, p10, p11});
            triangles.push_back({p00, p11, p01});
        }
    }
    return triangles;
}

//
template<typename T>
bool hasSelfIntersection(const tinynurbs::RationalSurface<T> &surface, int resolution_u = 50, int resolution_v = 50,
                         T epsilon = std::numeric_limits<T>::epsilon())
{
    auto triangles = discretizeSurface(surface, resolution_u, resolution_v);

    //
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        for (size_t j = i + 1; j < triangles.size(); ++j)
        {
            if (trianglesIntersect(triangles[i], triangles[j], epsilon))
            {
                return true;
            }
        }
    }
    return false;
}

#endif// ALGORITHM_H
