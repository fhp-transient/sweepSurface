/**
 * Functions for modifying NURBS curves and surfaces.
 *
 * Use of this source code is governed by a BSD-style license that can be found
 * in the LICENSE file.
 */

#ifndef TINYNURBS_MODIFY_H
#define TINYNURBS_MODIFY_H

#include <glm/glm.hpp>
#include "tinynurbs/core/check.h"
#include "tinynurbs/core/curve.h"
#include "tinynurbs/core/surface.h"
#include "tinynurbs/util/util.h"
#include <tuple>
#include <vector>

namespace tinynurbs {

    /////////////////////////////////////////////////////////////////////

    namespace internal {

        /**
         * Insert knots in the curve
         * @param[in] deg Degree of the curve
         * @param[in] knots Knot vector of the curve
         * @param[in] cp Control points of the curve
         * @param[in] u Parameter to insert knot(s) at
         * @param[in] r Number of times to insert knot
         * @param[out] new_knots Updated knot vector
         * @param[out] new_cp Updated control points
         */
        template <int dim, typename T>
        void curveKnotInsert(unsigned int deg, const std::vector<T>& knots,
                             const std::vector<glm::vec<dim, T>>& cp, T u,
                             unsigned int r, std::vector<T>& new_knots,
                             std::vector<glm::vec<dim, T>>& new_cp) {
            int k = findSpan(deg, knots, u);
            unsigned int s = knotMultiplicity(knots, u);
            assert(s <= deg); // Multiplicity cannot be greater than degree
            if (s == deg) {
                new_knots = knots;
                new_cp = cp;
                return;
            }
            if ((r + s) > deg) {
                r = deg - s;
            }

            // Insert new knots between k and (k + 1)
            new_knots.resize(knots.size() + r);
            for (int i = 0; i < k + 1; ++i) {
                new_knots[i] = knots[i];
            }
            for (unsigned int i = 1; i < r + 1; ++i) {
                new_knots[k + i] = u;
            }
            for (int i = k + 1; i < knots.size(); ++i) {
                new_knots[i + r] = knots[i];
            }
            // Copy unaffected control points
            new_cp.resize(cp.size() + r);
            for (int i = 0; i < k - deg + 1; ++i) {
                new_cp[i] = cp[i];
            }
            for (int i = k - s; i < cp.size(); ++i) {
                new_cp[i + r] = cp[i];
            }
            // Copy affected control points
            std::vector<glm::vec<dim, T>> tmp;
            tmp.resize(deg - s + 1);
            for (int i = 0; i < deg - s + 1; ++i) {
                tmp[i] = cp[k - deg + i];
            }
            // Modify affected control points
            for (int j = 1; j <= r; ++j) {
                int L = k - deg + j;
                for (int i = 0; i < deg - j - s + 1; ++i) {
                    T a = (u - knots[L + i]) /
                          (knots[i + k + 1] - knots[L + i]);
                    tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                }
                new_cp[L] = tmp[0];
                new_cp[k + r - j - s] = tmp[deg - j - s];
            }
            int L = k - deg + r;
            for (int i = L + 1; i < k - s; ++i) {
                new_cp[i] = tmp[i - L];
            }
        }

        /**
         * Insert knots in the surface along one direction
         * @param[in] degree Degree of the surface along which to insert knot
         * @param[in] knots Knot vector
         * @param[in] cp 2D array of control points
         * @param[in] knot Knot value to insert
         * @param[in] r Number of times to insert
         * @param[in] along_u Whether inserting along u-direction
         * @param[out] new_knots Updated knot vector
         * @param[out] new_cp Updated control points
         */
        template <int dim, typename T>
        void surfaceKnotInsert(unsigned int degree, const std::vector<T>& knots,
                               const array2<glm::vec<dim, T>>& cp, T knot,
                               unsigned int r, bool along_u,
                               std::vector<T>& new_knots,
                               array2<glm::vec<dim, T>>& new_cp) {
            int span = findSpan(degree, knots, knot);
            unsigned int s = knotMultiplicity(knots, knot);
            assert(s <=
                   degree); // Knot multiplicity cannot be greater than degree
            if (s == degree) {
                new_cp = cp;
                new_knots = knots;
                return;
            }
            if ((r + s) > degree) {
                r = degree - s;
            }

            // Create a new knot vector
            new_knots.resize(knots.size() + r);
            for (int i = 0; i <= span; ++i) {
                new_knots[i] = knots[i];
            }
            for (int i = 1; i <= r; ++i) {
                new_knots[span + i] = knot;
            }
            for (int i = span + 1; i < knots.size(); ++i) {
                new_knots[i + r] = knots[i];
            }
            // Compute alpha
            array2<T> alpha(degree - s, r + 1, T(0));
            for (int j = 1; j <= r; ++j) {
                int L = span - degree + j;
                for (int i = 0; i <= degree - j - s; ++i) {
                    alpha(i, j) = (knot - knots[L + i]) /
                                  (knots[i + span + 1] - knots[L + i]);
                }
            }

            // Create a temporary container for affected control points per
            // row/column
            std::vector<glm::vec<dim, T>> tmp(degree + 1);

            if (along_u) {
                // Create new control points with additional rows
                new_cp.resize(cp.rows() + r, cp.cols());

                // Update control points
                // Each row is a u-isocurve, each col is a v-isocurve
                for (int col = 0; col < cp.cols(); ++col) {
                    // Copy unaffected control points
                    for (int i = 0; i <= span - degree; ++i) {
                        new_cp(i, col) = cp(i, col);
                    }
                    for (int i = span - s; i < cp.rows(); ++i) {
                        new_cp(i + r, col) = cp(i, col);
                    }
                    // Copy affected control points to temp array
                    for (int i = 0; i < degree - s + 1; ++i) {
                        tmp[i] = cp(span - degree + i, col);
                    }
                    // Insert knot
                    for (int j = 1; j <= r; ++j) {
                        int L = span - degree + j;
                        for (int i = 0; i <= degree - j - s; ++i) {
                            T a = alpha(i, j);
                            tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                        }
                        new_cp(L, col) = tmp[0];
                        new_cp(span + r - j - s, col) = tmp[degree - j - s];
                    }
                    int L = span - degree + r;
                    for (int i = L + 1; i < span - s; ++i) {
                        new_cp(i, col) = tmp[i - L];
                    }
                }
            } else {
                // Create new control points with additional columns
                new_cp.resize(cp.rows(), cp.cols() + r);

                // Update control points
                // Each row is a u-isocurve, each col is a v-isocurve
                for (int row = 0; row < cp.rows(); ++row) {
                    // Copy unaffected control points
                    for (int i = 0; i <= span - degree; ++i) {
                        new_cp(row, i) = cp(row, i);
                    }
                    for (int i = span - s; i < cp.cols(); ++i) {
                        new_cp(row, i + r) = cp(row, i);
                    }
                    // Copy affected control points to temp array
                    for (int i = 0; i < degree - s + 1; ++i) {
                        tmp[i] = cp(row, span - degree + i);
                    }
                    // Insert knot
                    for (int j = 1; j <= r; ++j) {
                        int L = span - degree + j;
                        for (int i = 0; i <= degree - j - s; ++i) {
                            T a = alpha(i, j);
                            tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                        }
                        new_cp(row, L) = tmp[0];
                        new_cp(row, span + r - j - s) = tmp[degree - j - s];
                    }
                    int L = span - degree + r;
                    for (int i = L + 1; i < span - s; ++i) {
                        new_cp(row, i) = tmp[i - L];
                    }
                }
            }
        }

        /**
         * Split the curve into two
         * @param[in] degree Degree of curve
         * @param[in] knots Knot vector
         * @param[in] control_points Array of control points
         * @param[in] u Parameter to split curve
         * @param[out] left_knots Knots of the left part of the curve
         * @param[out] left_control_points Control points of the left part of
         * the curve
         * @param[out] right_knots Knots of the right part of the curve
         * @param[out] right_control_points Control points of the right part of
         * the curve
         */
        template <int dim, typename T>
        void curveSplit(unsigned int degree, const std::vector<T>& knots,
                        const std::vector<glm::vec<dim, T>>& control_points,
                        T u, std::vector<T>& left_knots,
                        std::vector<glm::vec<dim, T>>& left_control_points,
                        std::vector<T>& right_knots,
                        std::vector<glm::vec<dim, T>>& right_control_points) {
            std::vector<T> tmp_knots;
            std::vector<glm::vec<dim, T>> tmp_cp;

            int span = findSpan(degree, knots, u);
            int r = degree - knotMultiplicity(knots, u);

            internal::curveKnotInsert(degree, knots, control_points, u, r,
                                      tmp_knots, tmp_cp);

            left_knots.clear();
            right_knots.clear();
            left_control_points.clear();
            right_control_points.clear();

            int span_l = findSpan(degree, tmp_knots, u) + 1;
            for (int i = 0; i < span_l; ++i) {
                left_knots.push_back(tmp_knots[i]);
            }
            left_knots.push_back(u);

            for (int i = 0; i < degree + 1; ++i) {
                right_knots.push_back(u);
            }
            for (int i = span_l; i < tmp_knots.size(); ++i) {
                right_knots.push_back(tmp_knots[i]);
            }

            int ks = span - degree + 1;
            for (int i = 0; i < ks + r; ++i) {
                left_control_points.push_back(tmp_cp[i]);
            }
            for (int i = ks + r - 1; i < tmp_cp.size(); ++i) {
                right_control_points.push_back(tmp_cp[i]);
            }
        }

        /**
         * Split the surface into two along given parameter direction
         * @param[in] degree Degree of surface along given direction
         * @param[in] knots Knot vector of surface along given direction
         * @param[in] control_points Array of control points
         * @param[in] param Parameter to split curve
         * @param[in] along_u Whether the direction to split along is the
         * u-direction
         * @param[out] left_knots Knots of the left part of the curve
         * @param[out] left_control_points Control points of the left part of
         * the curve
         * @param[out] right_knots Knots of the right part of the curve
         * @param[out] right_control_points Control points of the right part of
         * the curve
         */
        template <int dim, typename T>
        void surfaceSplit(unsigned int degree, const std::vector<T>& knots,
                          const array2<glm::vec<dim, T>>& control_points,
                          T param, bool along_u, std::vector<T>& left_knots,
                          array2<glm::vec<dim, T>>& left_control_points,
                          std::vector<T>& right_knots,
                          array2<glm::vec<dim, T>>& right_control_points) {
            std::vector<T> tmp_knots;
            array2<glm::vec<dim, T>> tmp_cp;

            int span = findSpan(degree, knots, param);
            unsigned int r = degree - knotMultiplicity(knots, param);
            internal::surfaceKnotInsert(degree, knots, control_points, param, r,
                                        along_u, tmp_knots, tmp_cp);

            left_knots.clear();
            right_knots.clear();

            int span_l = findSpan(degree, tmp_knots, param) + 1;
            for (int i = 0; i < span_l; ++i) {
                left_knots.push_back(tmp_knots[i]);
            }
            left_knots.push_back(param);

            for (int i = 0; i < degree + 1; ++i) {
                right_knots.push_back(param);
            }
            for (int i = span_l; i < tmp_knots.size(); ++i) {
                right_knots.push_back(tmp_knots[i]);
            }

            int ks = span - degree + 1;
            if (along_u) {
                size_t ii = 0;
                left_control_points.resize(ks + r, tmp_cp.cols());
                for (int i = 0; i < ks + r; ++i) {
                    for (int j = 0; j < tmp_cp.cols(); ++j) {
                        left_control_points[ii++] = tmp_cp(i, j);
                    }
                }
                ii = 0;
                right_control_points.resize(tmp_cp.rows() - ks - r + 1,
                                            tmp_cp.cols());
                for (int i = ks + r - 1; i < tmp_cp.rows(); ++i) {
                    for (int j = 0; j < tmp_cp.cols(); ++j) {
                        right_control_points[ii++] = tmp_cp(i, j);
                    }
                }
            } else {
                size_t ii = 0;
                left_control_points.resize(tmp_cp.rows(), ks + r);
                for (int i = 0; i < tmp_cp.rows(); ++i) {
                    for (int j = 0; j < ks + r; ++j) {
                        left_control_points[ii++] = tmp_cp(i, j);
                    }
                }
                ii = 0;
                right_control_points.resize(tmp_cp.rows(),
                                            tmp_cp.cols() - ks - r + 1);
                for (int i = 0; i < tmp_cp.rows(); ++i) {
                    for (int j = ks + r - 1; j < tmp_cp.cols(); ++j) {
                        right_control_points[ii++] = tmp_cp(i, j);
                    }
                }
            }
        }

        /**
         * Knot vector refinement for a nonrational NURBS curve.
         * @param[in] degree Degree of the curve.
         * @param[in] knots Original knot vector.
         * @param[in] control_points Original control points.
         * @param[in] knots_diff Knots to be inserted.
         * @param[out] new_knots Refined knot vector.
         * @param[out] new_control_points Control points of the refined curve.
         */
        template <int dim, typename T>
        void
            RefineKnotVectCurve(int degree, const std::vector<T>& knots,
                const std::vector<glm::vec<dim, T>>& control_points,
                const std::vector<T>& knots_diff,
                std::vector<T>& new_knots,
                std::vector<glm::vec<dim, T>>& new_control_points) {
            if (knots_diff.size() == 0) {
                new_control_points = control_points;
                new_knots = knots;
                return;
            }
            assert(std::is_sorted(knots_diff.begin(), knots_diff.end()));
            int a = findSpan(degree, knots, knots_diff[0]),
                b = findSpan(degree, knots, knots_diff[knots_diff.size() - 1]) +
                1;
            new_knots.resize(knots.size() + knots_diff.size());
            new_control_points.resize(new_knots.size() - degree - 1);
            // Save unaltered control points and knots.
            std::copy(control_points.begin(),
                control_points.begin() + a - degree + 1,
                new_control_points.begin());
            std::copy(control_points.begin() + b - 1, control_points.end(),
                new_control_points.begin() + b + knots_diff.size() - 1);
            std::copy(knots.begin(), knots.begin() + a + 1, new_knots.begin());
            std::copy(knots.begin() + b + degree, knots.end(),
                new_knots.begin() + b + degree + knots_diff.size());
            int i = b + degree - 1, k = b + degree + knots_diff.size() - 1;
            for (int j = knots_diff.size() - 1; j >= 0; j--) {
                while (knots_diff[j] <
                    knots[i] + std::numeric_limits<T>::epsilon() &&
                    i > a) {
                    new_control_points[k - degree - 1] =
                        control_points[i - degree - 1];
                    new_knots[k] = knots[i];
                    k--, i--;
                }
                new_control_points[k - degree - 1] =
                    new_control_points[k - degree];
                for (int l = 1; l <= degree; l++) {
                    int ind = k - degree + l;
                    if (util::close(new_knots[k + l], knots_diff[j])) {
                        new_control_points[ind - 1] = new_control_points[ind];
                    }
                    else {
                        // Compute the new control points.
                        T alfa = (new_knots[k + l] - knots_diff[j]) /
                            (new_knots[k + l] - knots[i - degree + l]);
                        new_control_points[ind - 1] =
                            alfa * new_control_points[ind - 1] +
                            (1 - alfa) * new_control_points[ind];
                    }
                }
                new_knots[k] = knots_diff[j];
                k--;
            }
        }

        /**
         * Degree elevation for a nonrational NURBS curve.
         * @param[in] degree Degree of the original curve.
         * @param[in] knots Original knot vector.
         * @param[in] control_points Original control points.
         * @param[in] degree_inc Degree increment.
         * @param[out] new_knots Knot vector of the resulting curve.
         * @param[out] new_control_points Control points of the resulting curve.
         */
        template <int dim, typename T>
        void DegreeElevateCurve(int degree, const std::vector<T>& knots,
            const std::vector<glm::vec<dim, T>>& control_points,
            int degree_inc, std::vector<T>& new_knots,
            std::vector<glm::vec<dim, T>>& new_control_points) {
            assert(degree_inc >= 0);
            if (degree_inc == 0) {
                new_knots = knots;
                new_control_points = control_points;
                return;
            }
            new_knots.resize(control_points.size() +
                (control_points.size() + 4 - degree) * degree_inc);
            new_control_points.resize(control_points.size() + degree +
                (control_points.size() + 5 - degree) *
                degree_inc);
            int degree_h = degree + degree_inc;
            int degree_h_2 = degree_h / 2;
            // Coefficients for degree elevating the Bezier segments
            array2<T> bezalfs(degree_h + 1, degree + 1);
            // Bezier control points of the current segment
            typedef glm::vec<dim, T> tvecn;
            std::vector<tvecn> bpts(degree + 1);
            // Elevated Bezier control points of the current segment
            std::vector<tvecn> ebpts(degree_h + 1);
            // Left most control points of the next Bezier segment
            std::vector<tvecn> Nextbpts(degree - 1);
            // Knot insertion alphas
            std::vector<T> alphas(degree - 1);
            // Compute Bezier degree elevation coefficients
            bezalfs(0, 0) = bezalfs(degree_h, degree) = 1;
            for (int i = 1; i <= degree_h_2; i++) {
                T inv = 1.0 / util::binomial(degree_h, i);
                int mpi = std::min(degree, i);
                for (int j = std::max(0, i - degree_inc); j <= mpi; j++) {
                    bezalfs(i, j) = inv * util::binomial(degree, j) *
                        util::binomial(degree_inc, i - j);
                }
            }
            for (int i = degree_h_2 + 1; i <= degree_h - 1; i++) {
                int mpi = std::min(degree, i);
                for (int j = std::max(0, i - degree_inc); j <= mpi; j++) {
                    bezalfs(i, j) = bezalfs(degree_h - i, degree - j);
                }
            }
            int mh = degree_h;
            T ua = knots[0];
            std::fill(new_knots.begin(), new_knots.begin() + degree_h + 1, ua);
            // Initialize first Bezier seg
            int kind = degree_h + 1, cind = 1, lbz, rbz, r = -1;
            new_control_points[0] = control_points[0];
            std::copy(control_points.begin(),
                control_points.begin() + degree + 1, bpts.begin());
            int a = degree, b = degree + 1;
            // Big loop through knot vector
            while (b < knots.size() - 1) {
                int i = b;
                while (b < knots.size() - 1 &&
                    util::close(knots[b], knots[b + 1])) {
                    b++;
                }
                int mul = b - i + 1, oldr = r;
                T ub = knots[b];
                mh = mh + mul + degree_inc;
                // Insert `knots[b]` `r` times
                r = degree - mul;
                lbz = oldr > 0 ? (oldr + 2) / 2 : 1;
                rbz = r > 0 ? degree_h - (r + 1) / 2 : degree_h;
                if (r > 0) {
                    // Insert knot to get Bezier segment
                    T numer = ub - ua;
                    for (int k = degree; k > mul; k--) {
                        alphas[k - mul - 1] = numer / (knots[a + k] - ua);
                    }
                    for (int j = 1; j <= r; j++) {
                        int save = r - j, s = mul + j;
                        for (int k = degree; k >= s; k--) {
                            T alf = alphas[k - s];
                            bpts[k] = alf * bpts[k] + ((T)1.0 - alf) * bpts[k - 1];
                        }
                        Nextbpts[save] = bpts[degree];
                    }
                } // End of "insert knot"
                // Degree elevate Bezier
                int RBZ = r > 0 ? rbz + 1 : rbz;
                for (int i = lbz; i <= RBZ; i++) {
                    ebpts[i] = tvecn(0);
                    int mpi = std::min(degree, i);
                    for (int j = std::max(0, i - degree_inc); j <= mpi; j++) {
                        ebpts[i] = ebpts[i] + bezalfs(i, j) * bpts[j];
                    }
                } // End of degree elevating Bezier
                if (oldr > 1) {
                    // Remove `knot[a]` `oldr` times
                    T alfj;
                    if (oldr > 2) {
                        alfj = (ua - new_knots[kind - 1]) /
                            (ub - new_knots[kind - 1]);
                    }
                    int first = kind - 2, last = kind;
                    // Knot removal loop
                    for (int tr = 1; tr < oldr; tr++) {
                        int i = first, j = last;
                        int kj = j - kind + 1;
                        // Loop and compute the new control points for one
                        // removal step
                        while (j - i > tr) {
                            if (i < cind) {
                                T alf = (ua - new_knots[i]) /
                                    (ub - new_knots[i]);
                                new_control_points[i] =
                                    (new_control_points[i] -
                                        ((T)1.0 - alf) * new_control_points[i - 1]) /
                                    alf;
                            }
                            if (kj >= lbz) {
                                ebpts[kj] = ebpts[kj] -
                                    alfj * ebpts[kj + 1] / ((T)1.0 - alfj);
                            }
                            i++, j--, kj--;
                        }
                        first--, last++;
                    }
                } // End of removing `knots[a]`
                if (a != degree) {
                    // Load the knot `ua`
                    std::fill(new_knots.begin() + kind,
                        new_knots.begin() + kind + degree_h - oldr, ua);
                    kind = kind + degree_h - oldr;
                }
                // Load control points into `new_control_points`
                std::copy(ebpts.begin() + lbz, ebpts.begin() + rbz + 1,
                    new_control_points.begin() + cind);
                cind = cind + rbz - lbz + 1;
                if (b < knots.size() - 1) {
                    // Set up for next pass through loop
                    lbz = 1;
                    std::copy(Nextbpts.begin(), Nextbpts.begin() + r,
                        bpts.begin());
                    std::copy(control_points.begin() + (b - degree + r),
                        control_points.begin() + b + 1, bpts.begin() + r);
                    a = b;
                    b++;
                    ua = ub;
                }
                else {
                    // End knot
                    std::fill(new_knots.begin() + kind,
                        new_knots.begin() + kind + degree_h + 1, ub);
                }
            } // End of while-loop (b < knots.size() - 1)
            new_knots.resize(mh + 1);
            new_control_points.resize(mh - degree_h);
        }


    } // namespace internal

    /////////////////////////////////////////////////////////////////////

    /**
     * Insert knots in the curve
     * @param[in] crv Curve object
     * @param[in] u Parameter to insert knot at
     * @param[in] repeat Number of times to insert
     * @return New curve with #repeat knots inserted at u
     */
    template <typename T>
    Curve<T> curveKnotInsert(const Curve<T>& crv, T u,
                             unsigned int repeat = 1) {
        Curve<T> new_crv;
        new_crv.degree = crv.degree;
        internal::curveKnotInsert(crv.degree, crv.knots, crv.control_points, u,
                                  repeat, new_crv.knots,
                                  new_crv.control_points);
        return new_crv;
    }

    /**
     * Insert knots in the rational curve
     * @param[in] crv RationalCurve object
     * @param[in] u Parameter to insert knot at
     * @param[in] repeat Number of times to insert
     * @return New RationalCurve object with #repeat knots inserted at u
     */
    template <typename T>
    RationalCurve<T> curveKnotInsert(const RationalCurve<T>& crv, T u,
                                     unsigned int repeat = 1) {
        RationalCurve<T> new_crv;
        new_crv.degree = crv.degree;

        // Convert to homogenous coordinates
        std::vector<glm::vec<4, T>> Cw;
        Cw.reserve(crv.control_points.size());
        for (int i = 0; i < crv.control_points.size(); ++i) {
            Cw.push_back(util::cartesianToHomogenous(crv.control_points[i],
                                                     crv.weights[i]));
        }

        // Perform knot insertion and get new knots and control points
        std::vector<glm::vec<4, T>> new_Cw;
        std::vector<T> new_knots;
        internal::curveKnotInsert(crv.degree, crv.knots, Cw, u, repeat,
                                  new_crv.knots, new_Cw);

        // Convert back to cartesian coordinates
        new_crv.control_points.reserve(new_Cw.size());
        new_crv.weights.reserve(new_Cw.size());
        for (int i = 0; i < new_Cw.size(); ++i) {
            new_crv.control_points.push_back(
                util::homogenousToCartesian(new_Cw[i]));
            new_crv.weights.push_back(new_Cw[i].w);
        }
        return new_crv;
    }

    /**
     * Insert knots in the surface along u-direction
     * @param[in] srf Surface object
     * @param[in] u Knot value to insert
     * @param[in] repeat Number of times to insert
     * @return New Surface object after knot insertion
     */
    template <typename T>
    Surface<T> surfaceKnotInsertU(const Surface<T>& srf, T u,
                                  unsigned int repeat = 1) {
        Surface<T> new_srf;
        new_srf.degree_u = srf.degree_u;
        new_srf.degree_v = srf.degree_v;
        new_srf.knots_v = srf.knots_v;
        internal::surfaceKnotInsert(new_srf.degree_u, srf.knots_u,
                                    srf.control_points, u, repeat, true,
                                    new_srf.knots_u, new_srf.control_points);
        return new_srf;
    }

    /**
     * Insert knots in the rational surface along u-direction
     * @param[in] srf RationalSurface object
     * @param[in] u Knot value to insert
     * @param[in] repeat Number of times to insert
     * @return New RationalSurface object after knot insertion
     */
    template <typename T>
    RationalSurface<T> surfaceKnotInsertU(const RationalSurface<T>& srf, T u,
                                          unsigned int repeat = 1) {
        RationalSurface<T> new_srf;
        new_srf.degree_u = srf.degree_u;
        new_srf.degree_v = srf.degree_v;
        new_srf.knots_v = srf.knots_v;

        // Original control points in homogenous coordinates
        array2<glm::vec<4, T>> Cw(srf.control_points.rows(),
                                  srf.control_points.cols());
        for (int i = 0; i < srf.control_points.rows(); ++i) {
            for (int j = 0; j < srf.control_points.cols(); ++j) {
                Cw(i, j) = util::cartesianToHomogenous(srf.control_points(i, j),
                                                       srf.weights(i, j));
            }
        }

        // New knots and new homogenous control points after knot insertion
        std::vector<T> new_knots_u;
        array2<glm::vec<4, T>> new_Cw;
        internal::surfaceKnotInsert(srf.degree_u, srf.knots_u, Cw, u, repeat,
                                    true, new_srf.knots_u, new_Cw);

        // Convert back to cartesian coordinates
        new_srf.control_points.resize(new_Cw.rows(), new_Cw.cols());
        new_srf.weights.resize(new_Cw.rows(), new_Cw.cols());
        for (int i = 0; i < new_Cw.rows(); ++i) {
            for (int j = 0; j < new_Cw.cols(); ++j) {
                new_srf.control_points(i, j) = util::homogenousToCartesian(
                    new_Cw(i, j));
                new_srf.weights(i, j) = new_Cw(i, j).w;
            }
        }
        return new_srf;
    }

    /**
     * Insert knots in the surface along v-direction
     * @param[in] srf Surface object
     * @param[in] v Knot value to insert
     * @param[in] repeat Number of times to insert
     * @return New Surface object after knot insertion
     */
    template <typename T>
    Surface<T> surfaceKnotInsertV(const Surface<T>& srf, T v,
                                  unsigned int repeat = 1) {
        Surface<T> new_srf;
        new_srf.degree_u = srf.degree_u;
        new_srf.degree_v = srf.degree_v;
        new_srf.knots_u = srf.knots_u;
        // New knots and new control points after knot insertion
        internal::surfaceKnotInsert(srf.degree_v, srf.knots_v,
                                    srf.control_points, v, repeat, false,
                                    new_srf.knots_v, new_srf.control_points);
        return new_srf;
    }

    /**
     * Insert knots in the rational surface along v-direction
     * @param[in] srf RationalSurface object
     * @param[in] v Knot value to insert
     * @param[in] repeat Number of times to insert
     * @return New RationalSurface object after knot insertion
     */
    template <typename T>
    RationalSurface<T> surfaceKnotInsertV(const RationalSurface<T>& srf, T v,
                                          unsigned int repeat = 1) {
        RationalSurface<T> new_srf;
        new_srf.degree_u = srf.degree_u;
        new_srf.degree_v = srf.degree_v;
        new_srf.knots_u = srf.knots_u;
        // Original control points in homogenous coordinates
        array2<glm::vec<4, T>> Cw(srf.control_points.rows(),
                                  srf.control_points.cols());
        for (int i = 0; i < srf.control_points.rows(); ++i) {
            for (int j = 0; j < srf.control_points.cols(); ++j) {
                Cw(i, j) = util::cartesianToHomogenous(srf.control_points(i, j),
                                                       srf.weights(i, j));
            }
        }

        // New knots and new homogenous control points after knot insertion
        std::vector<T> new_knots_v;
        array2<glm::vec<4, T>> new_Cw;
        internal::surfaceKnotInsert(srf.degree_v, srf.knots_v, Cw, v, repeat,
                                    false, new_srf.knots_v, new_Cw);

        // Convert back to cartesian coordinates
        new_srf.control_points.resize(new_Cw.rows(), new_Cw.cols());
        new_srf.weights.resize(new_Cw.rows(), new_Cw.cols());
        for (int i = 0; i < new_Cw.rows(); ++i) {
            for (int j = 0; j < new_Cw.cols(); ++j) {
                new_srf.control_points(i, j) = util::homogenousToCartesian(
                    new_Cw(i, j));
                new_srf.weights(i, j) = new_Cw(i, j).w;
            }
        }
        return new_srf;
    }

    /**
     * Split a curve into two
     * @param[in] crv Curve object
     * @param[in] u Parameter to split at
     * @return Tuple with first half and second half of the curve
     */
    template <typename T>
    std::tuple<Curve<T>, Curve<T>> curveSplit(const Curve<T>& crv, T u) {
        Curve<T> left, right;
        left.degree = crv.degree;
        right.degree = crv.degree;
        internal::curveSplit(crv.degree, crv.knots, crv.control_points, u,
                             left.knots, left.control_points, right.knots,
                             right.control_points);
        return std::make_tuple(std::move(left), std::move(right));
    }

    /**
     * Split a rational curve into two
     * @param[in] crv RationalCurve object
     * @param[in] u Parameter to split at
     * @return Tuple with first half and second half of the curve
     */
    template <typename T>
    std::tuple<RationalCurve<T>, RationalCurve<T>>
    curveSplit(const RationalCurve<T>& crv, T u) {
        RationalCurve<T> left, right;
        left.degree = crv.degree;
        right.degree = crv.degree;

        std::vector<glm::vec<4, T>> Cw, left_Cw, right_Cw;
        Cw.reserve(crv.control_points.size());
        for (int i = 0; i < crv.control_points.size(); ++i) {
            Cw.push_back(util::cartesianToHomogenous(crv.control_points[i],
                                                     crv.weights[i]));
        }

        internal::curveSplit(crv.degree, crv.knots, Cw, u, left.knots, left_Cw,
                             right.knots, right_Cw);

        left.control_points.reserve(left_Cw.size());
        left.weights.reserve(left_Cw.size());
        right.control_points.reserve(right_Cw.size());
        right.weights.reserve(right_Cw.size());
        for (int i = 0; i < left_Cw.size(); ++i) {
            left.control_points.push_back(
                util::homogenousToCartesian(left_Cw[i]));
            left.weights.push_back(left_Cw[i].w);
        }
        for (int i = 0; i < right_Cw.size(); ++i) {
            right.control_points.push_back(
                util::homogenousToCartesian(right_Cw[i]));
            right.weights.push_back(right_Cw[i].w);
        }
        return std::make_tuple(std::move(left), std::move(right));
    }

    /**
     * Split a surface into two along u-direction
     * @param[in] srf Surface object
     * @param[in] u Parameter along u-direction to split the surface
     * @return Tuple with first and second half of the surfaces
     */
    template <typename T>
    std::tuple<Surface<T>, Surface<T>> surfaceSplitU(const Surface<T>& srf,
                                                     T u) {
        Surface<T> left, right;
        left.degree_u = srf.degree_u;
        left.degree_v = srf.degree_v;
        left.knots_v = srf.knots_v;
        right.degree_u = srf.degree_u;
        right.degree_v = srf.degree_v;
        right.knots_v = srf.knots_v;
        internal::surfaceSplit(srf.degree_u, srf.knots_u, srf.control_points, u,
                               true, left.knots_u, left.control_points,
                               right.knots_u, right.control_points);
        return std::make_tuple(std::move(left), std::move(right));
    }

    /**
     * Split a rational surface into two along u-direction
     * @param[in] srf RationalSurface object
     * @param[in] u Parameter along u-direction to split the surface
     * @return Tuple with first and second half of the surfaces
     */
    template <typename T>
    std::tuple<RationalSurface<T>, RationalSurface<T>>
    surfaceSplitU(const RationalSurface<T>& srf, T u) {
        RationalSurface<T> left, right;
        left.degree_u = srf.degree_u;
        left.degree_v = srf.degree_v;
        left.knots_v = srf.knots_v;
        right.degree_u = srf.degree_u;
        right.degree_v = srf.degree_v;
        right.knots_v = srf.knots_v;

        // Compute homogenous coordinates of control points and weights
        array2<glm::vec<4, T>> Cw = util::cartesianToHomogenous(
            srf.control_points, srf.weights);

        // Split surface with homogenous coordinates
        array2<glm::vec<4, T>> left_Cw, right_Cw;
        internal::surfaceSplit(srf.degree_u, srf.knots_u, Cw, u, true,
                               left.knots_u, left_Cw, right.knots_u, right_Cw);

        // Convert back to cartesian coordinates
        util::homogenousToCartesian(left_Cw, left.control_points, left.weights);
        util::homogenousToCartesian(right_Cw, right.control_points,
                                    right.weights);

        return std::make_tuple(std::move(left), std::move(right));
    }

    /**
     * Split a surface into two along v-direction
     * @param[in] srf Surface object
     * @param[in] v Parameter along v-direction to split the surface
     * @return Tuple with first and second half of the surfaces
     */
    template <typename T>
    std::tuple<Surface<T>, Surface<T>> surfaceSplitV(const Surface<T>& srf,
                                                     T v) {
        Surface<T> left, right;
        left.degree_u = srf.degree_u;
        left.degree_v = srf.degree_v;
        left.knots_u = srf.knots_u;
        right.degree_u = srf.degree_u;
        right.degree_v = srf.degree_v;
        right.knots_u = srf.knots_u;
        internal::surfaceSplit(srf.degree_v, srf.knots_v, srf.control_points, v,
                               false, left.knots_v, left.control_points,
                               right.knots_v, right.control_points);
        return std::make_tuple(std::move(left), std::move(right));
    }

    /**
     * Split a rational surface into two along v-direction
     * @param[in] srf RationalSurface object
     * @param[in] v Parameter along v-direction to split the surface
     * @return Tuple with first and second half of the surfaces
     */
    template <typename T>
    std::tuple<RationalSurface<T>, RationalSurface<T>>
    surfaceSplitV(const RationalSurface<T>& srf, T v) {
        RationalSurface<T> left, right;
        left.degree_u = srf.degree_u;
        left.degree_v = srf.degree_v;
        left.knots_u = srf.knots_u;
        right.degree_u = srf.degree_u;
        right.degree_v = srf.degree_v;
        right.knots_u = srf.knots_u;

        // Compute homogenous coordinates of control points and weights
        array2<glm::vec<4, T>> Cw = util::cartesianToHomogenous(
            srf.control_points, srf.weights);

        // Split surface with homogenous coordinates
        array2<glm::vec<4, T>> left_Cw, right_Cw;
        internal::surfaceSplit(srf.degree_v, srf.knots_v, Cw, v, false,
                               left.knots_v, left_Cw, right.knots_v, right_Cw);

        // Convert back to cartesian coordinates
        util::homogenousToCartesian(left_Cw, left.control_points, left.weights);
        util::homogenousToCartesian(right_Cw, right.control_points,
                                    right.weights);

        return std::make_tuple(std::move(left), std::move(right));
    }

    /**
     * Knot vector refinement for a polynomial NURBS curve.
     * @param[in] crv Original curve.
     * @param[in] knots_diff Knots to be inserted.
     * @return new_crv Curve after knot refinement.
     */
    template <typename T>
    Curve<T> RefineKnotVectCurve(const Curve<T>& crv,
        const std::vector<T>& knots_diff) {
        Curve<T> new_crv;
        new_crv.degree = crv.degree;
        internal::RefineKnotVectCurve(crv.degree, crv.knots, crv.control_points,
            knots_diff, new_crv.knots,
            new_crv.control_points);
        return new_crv;
    }

    /**
     * Knot vector refinement for a rational NURBS curve.
     * @param[in] crv Original curve.
     * @param[in] knots_diff Knots to be inserted.
     * @return new_crv Curve after knot refinement.
     */
    template <typename T>
    RationalCurve<T> RefineKnotVectCurve(const RationalCurve<T>& crv,
        const std::vector<T>& knots_diff) {
        std::vector<glm::vec<4, T>> Cw;
        Cw.reserve(crv.control_points.size());
        for (int i = 0; i < crv.control_points.size(); i++) {
            // Compute homogenous coordinates of control points
            Cw.push_back(util::cartesianToHomogenous(crv.control_points[i],
                crv.weights[i]));
        }
        RationalCurve<T> new_crv;
        new_crv.degree = crv.degree;
        std::vector<glm::vec<4, T>> new_Cw;
        // Knot insertion using homogenous coordinates
        internal::RefineKnotVectCurve(crv.degree, crv.knots, Cw, knots_diff,
            new_crv.knots, new_Cw);
        for (const auto& val : new_Cw) {
            // Convert back to Cartesian coordinates
            new_crv.control_points.push_back(util::homogenousToCartesian(val));
            new_crv.weights.push_back(val.w);
        }
        return new_crv;
    }

    /**
     * Degree elevation for a nonrational NURBS curve.
     * @param[in] crv Original curve.
     * @param[in] degree_inc Degree increment.
     * @return new_crv Curve after degree elevation
     */
    template <typename T>
    Curve<T> DegreeElevateCurve(const Curve<T>& crv, int degree_inc) {
        Curve<T> new_crv;
        new_crv.degree = crv.degree + degree_inc;
        internal::DegreeElevateCurve(crv.degree, crv.knots, crv.control_points,
            degree_inc, new_crv.knots,
            new_crv.control_points);
        return new_crv;
    }

    /**
     * Degree elevation for a rational NURBS curve.
     * @param[in] crv Original curve.
     * @param[in] degree_inc Degree increment.
     * @return new_crv Curve after degree elevation
     */
    template <typename T>
    RationalCurve<T> DegreeElevateCurve(const RationalCurve<T>& crv,
        int degree_inc) {
        std::vector<glm::vec<4, T>> Cw;
        Cw.reserve(crv.control_points.size());
        for (int i = 0; i < crv.control_points.size(); i++) {
            // Compute homogenous coordinates of control points
            Cw.push_back(util::cartesianToHomogenous(crv.control_points[i],
                crv.weights[i]));
        }
        RationalCurve<T> new_crv;
        new_crv.degree = crv.degree + degree_inc;
        std::vector<glm::vec<4, T>> new_Cw;
        // Knot insertion using homogenous coordinates
        internal::DegreeElevateCurve(crv.degree, crv.knots, Cw, degree_inc,
            new_crv.knots, new_Cw);
        for (const auto& val : new_Cw) {
            // Convert back to Cartesian coordinates
            new_crv.control_points.push_back(util::homogenousToCartesian(val));
            new_crv.weights.push_back(val.w);
        }
        return new_crv;
    }
} // namespace tinynurbs

#endif // TINYNURBS_MODIFY_H
