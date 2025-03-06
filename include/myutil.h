//
// Created by fhp on 24-12-11.
//

#ifndef MYUTIL_H
#define MYUTIL_H

#include <vector>
#include <glm/ext/scalar_constants.hpp>

template<typename T>
bool equalTo(T a, T b, T eps = std::numeric_limits<T>::epsilon())
{
    return (std::abs(a - b) < eps);
}

template<typename T>
bool IsGreaterThan(T value1, T value2, T eps = std::numeric_limits<T>::epsilon())
{
    return value1 > value2 && !equalTo(value1, value2, eps);
}

template<typename T>
bool IsGreaterThanOrEqual(T value1, T value2, T tolerance = std::numeric_limits<T>::epsilon())
{
    return (value1 - value2 > tolerance) || equalTo(value1, value2, tolerance);
}

template<typename T>
bool IsLessThan(T value1, T value2, T tolerance = std::numeric_limits<T>::epsilon())
{
    return value1 < value2 && !equalTo(value1, value2, tolerance);
}

template<typename T>
bool IsLessThanOrEqual(T value1, T value2, T tolerance = std::numeric_limits<T>::epsilon())
{
    return value1 < value2 || equalTo(value1, value2, tolerance);
}

template<typename T>
void InsertMidKnotCore(std::vector<T> &unqiueKnotVector, std::vector<T> &insert, int limitNumber)
{
    if (insert.size() == limitNumber)
    {
        return;
    }
    else
    {
        T standard = std::numeric_limits<T>::epsilon();
        int index = -1;
        for (int i = 0; i < unqiueKnotVector.size() - 1; i++)
        {
            T delta = unqiueKnotVector[i + 1] - unqiueKnotVector[i];
            if (IsGreaterThan(delta, standard))
            {
                standard = delta;
                index = i;
            }
        }
        T current = unqiueKnotVector[index] + standard / 2.0;
        unqiueKnotVector.emplace_back(current);
        std::sort(unqiueKnotVector.begin(), unqiueKnotVector.end());
        insert.emplace_back(current);

        InsertMidKnotCore(unqiueKnotVector, insert, limitNumber);
    }
}

template<typename T>
std::vector<T> getInsertedKnotElements(int insertKnotsNumber, const std::vector<T> &knotVector)
{
    std::vector<T> unqiueKnotVector = knotVector;
    unqiueKnotVector.erase(std::unique(unqiueKnotVector.begin(), unqiueKnotVector.end()), unqiueKnotVector.end());

    std::vector<T> insert;
    InsertMidKnotCore(unqiueKnotVector, insert, insertKnotsNumber);
    std::sort(insert.begin(), insert.end());
    return insert;
}

template<typename T>
bool IsZero(const glm::vec<3, T> &vec, T epsilon = glm::epsilon<T>())
{
    return glm::length(vec) <= epsilon;
}

enum class IntegratorType :int
{
    Simpson = 0,
    GaussLegendre = 1,
    Chebyshev = 2,
};

static const int NURBSMaxDegree = 7;

template<typename T>
static void Transpose(const std::vector<std::vector<T> > &matrix, std::vector<std::vector<T> > &transposed)
{
    std::vector<T> temp;

    for (int i = 0; i < matrix[0].size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
        {
            temp.emplace_back(matrix[j][i]);
        }
        transposed.emplace_back(temp);
        temp.erase(temp.begin(), temp.end());
    }
}

template<typename T>
std::vector<std::vector<T> > MatrixMultiply(const std::vector<std::vector<T> > &left,
                                            const std::vector<std::vector<T> > &right)
{
    int m = left.size();
    int n = left[0].size();
    int p = right[0].size();

    std::vector<std::vector<T> > result(m, std::vector<T>(p, 0.0));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            for (int k = 0; k < n; k++)
            {
                result[i][j] += left[i][k] * right[k][j];
            }
        }
    }
    return result;
}

template<typename T>
static std::vector<T> GetColumn(const std::vector<std::vector<T> > &matrix, int columnIndex)
{
    int size = matrix.size();
    std::vector<T> result(size);
    for (int i = 0; i < size; i++)
    {
        result[i] = matrix[i][columnIndex];
    }
    return result;
}

template<typename T>
void transpose(std::vector<std::vector<T>>& matrix) {
    if (matrix.empty() || matrix[0].empty()) return;

    int rows = matrix.size();
    int cols = matrix[0].size();

    // 创建一个新的转置矩阵
    std::vector<std::vector<T>> transposed(cols, std::vector<T>(rows));

    // 进行转置
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }

    // 替换原矩阵
    matrix = transposed;
}

#endif //MYUTIL_H
