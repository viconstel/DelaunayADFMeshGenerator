#include "CoxDeBoorAlgorithm.h"

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <cstring>


CoxDeBoorAlgorithm::CoxDeBoorAlgorithm()
    : m_dim(0)
    , m_spline_order(0)
    , m_knots_count(0)
    , m_knot_vector(nullptr)
    , m_control_points_count(0)
    , m_control_points(nullptr)
    , m_threads_count(0)
{}


CoxDeBoorAlgorithm::CoxDeBoorAlgorithm(unsigned      p_dim,
                                       unsigned      p_spline_order,
                                       unsigned      p_knots_count,
                                       const double *p_knot_vector,
                                       unsigned      p_control_points_count,
                                       const double *p_control_points)
    : m_dim(p_dim)
    , m_spline_order(p_spline_order)
    , m_knots_count(p_knots_count)
    , m_knot_vector(p_knot_vector)
    , m_control_points_count(p_control_points_count)
    , m_control_points(p_control_points)
    , m_threads_count(0)
{
    initializeDependentData();
}


void CoxDeBoorAlgorithm::init(unsigned      p_dim,
                              unsigned      p_spline_order,
                              unsigned      p_knots_count,
                              const double *p_knot_vector,
                              unsigned      p_control_points_count,
                              const double *p_control_points)
{
    m_dim                  = p_dim;
    m_spline_order         = p_spline_order;
    m_knots_count          = p_knots_count;
    m_knot_vector          = p_knot_vector;
    m_control_points_count = p_control_points_count;
    m_control_points       = p_control_points;

    initializeDependentData();
}


void CoxDeBoorAlgorithm::initializeDependentData()
{
    m_max_der_order = 20;
    m_max_dim = 4;

    if ((m_dim > m_max_dim) || // Требуется увеличить значение m_max_dim
        (m_spline_order > m_max_der_order)) // Требуется увеличить значение m_max_der_order
        abort();

    if (m_threads_count == 0) {

        m_threads_count = omp_get_max_threads();

        m_compressed_control_points.resize(m_threads_count * m_max_dim * m_max_der_order);
        // Выделение памяти под структуру данных для рекурсивного вычисления контрольных точек
        m_recursion_control_points.resize(m_threads_count * (m_max_der_order + 1) * m_max_dim * (m_max_der_order + 1));

        m_cox_de_boor_algorithm.resize(m_threads_count);

        m_compressed_control_points_array = m_compressed_control_points.data();
        m_recursion_control_points_array = m_recursion_control_points.data();
        m_cox_de_boor_algorithm_array = m_cox_de_boor_algorithm.data();
    }
}


//void CoxDeBoorAlgorithm::clear()
//{
//    m_dim                  = 0;
//    m_spline_order         = 0;
//    m_knots_count          = 0;
//    m_knot_vector          = nullptr;
//    m_control_points_count = 0;
//    m_control_points       = nullptr;
//}


int CoxDeBoorAlgorithm::getPoint(double p_param, double *o_point) const
{
#   define _(point, coord) (m_dim * (point) + (coord))

    unsigned  param_ind = getParamInd(p_param);
    double   *compressed_control_points = &m_compressed_control_points_array[omp_get_thread_num() * m_max_dim * m_max_der_order];

//    for (unsigned p = 0; p < m_spline_order; ++p)
//        for (unsigned d = 0; d < m_dim; ++d)
//            compressed_control_points[_(p, d)] =
//                m_control_points[_(param_ind - m_spline_order + 1 + p, d)];
    memcpy( compressed_control_points, &m_control_points[_(param_ind - m_spline_order + 1, 0)], m_dim * m_spline_order * sizeof(double) );
    
    unsigned i;
    double val_1, val_2;
    for (unsigned r = 1; r < m_spline_order; ++r)
        for (unsigned j = m_spline_order - 1; j >= r; --j)
        {
            i = param_ind - m_spline_order + 1 + j;
            val_1 = p_param - m_knot_vector[i];
            val_2 = m_knot_vector[i + m_spline_order - r] - p_param;

            for (unsigned d = 0; d < m_dim; ++d)
                compressed_control_points[_(j, d)] =
                    (val_1 * compressed_control_points[_(j, d)] +
                     val_2 * compressed_control_points[_(j - 1, d)]) /
                    (val_1 + val_2);
        }

//    std::copy(compressed_control_points + m_dim * (m_spline_order - 1),
//              compressed_control_points + m_dim * m_spline_order,
//              o_point);
    memcpy( o_point, compressed_control_points + m_dim * (m_spline_order - 1), m_dim * sizeof(double) );

#   undef _
    return 0;
}


int CoxDeBoorAlgorithm::getDerPoint(unsigned p_der_order, double p_param, double *o_der_point) const
{
    if (p_der_order >= m_spline_order) {
        for (unsigned d = 0; d < m_dim; ++d)
            o_der_point[d] = 0.0;
        return 0;
    }

#   define _(point, coord) (m_dim * (point) + (coord))
#   define _RECURSION_CONTROL_POINTS_INDEX(l, coord) ((coord) + (l) * (m_max_dim * (m_max_der_order + 1)))

    unsigned param_ind = getParamInd(p_param);
    unsigned i         = param_ind - m_spline_order + p_der_order + 1;
    unsigned res_ind   = 0;
    double   delta     = 0.;

//    // Выделение памяти под структуру данных для рекурсивного вычисления контрольных точек
//    double **recursion_control_points = new double*[p_der_order + 1];
//    for (unsigned j = 0; j < p_der_order; ++j)
//        recursion_control_points[j] = new double[m_dim * (p_der_order - j + 1)];

//    // Следующий элемент хранит все (m_SplineOrder - p_DerOrder) вычисленных контрольных точек
//    recursion_control_points[p_der_order] = new double[m_dim * (m_spline_order - p_der_order)];

    double *recursion_control_points = &m_recursion_control_points_array[omp_get_thread_num() * (m_max_der_order + 1) * m_max_dim * (m_max_der_order + 1)];

    if ((p_der_order > m_max_der_order)) {//требуется увеличить значение m_max_der_order
        abort();
    }

    for (; i < param_ind + 1; ++i, ++res_ind)
    {
//        std::copy(&m_control_points[m_dim * (i + 1)] - m_dim * (p_der_order + 1),
//                  &m_control_points[m_dim * (i + 1)],
//                  recursion_control_points[0]);
        memcpy(&recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(0,0)], &m_control_points[m_dim * (i + 1)] - m_dim * (p_der_order + 1), m_dim * (p_der_order + 1) * sizeof(double) );

        for (unsigned j = 1; j < p_der_order + 1; ++j)
        {
            if (j != p_der_order)
            {
                for (unsigned k = 0; k < p_der_order - j + 1; ++k)
                {
                    delta = m_knot_vector[i + m_spline_order - (p_der_order - k - j + 1)] -
                            m_knot_vector[i - (p_der_order - k - j)];

                    for (unsigned d = 0; d < m_dim; ++d)
                        recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(j,_(k, d))] =
                            (m_spline_order - j) *
                            ((recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(j - 1,_(k + 1, d))] -
                                recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(j - 1,_(k, d))]) / delta);
                }
            }
            else
            {
                delta = m_knot_vector[i + m_spline_order - (p_der_order - j + 1)] -
                        m_knot_vector[i - (p_der_order - j)];

                for (unsigned d = 0; d < m_dim; ++d)
                    recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(j,_(res_ind, d))] =
                        (m_spline_order - j) *
                        ((recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(j - 1,_(1, d))] -
                          recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(j - 1,_(0, d))]) / delta);
            }
        }
    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
                m_dim,
                m_spline_order - p_der_order,
                2 * (m_spline_order - p_der_order),
                m_knot_vector + param_ind - m_spline_order + p_der_order + 1,
                m_spline_order - p_der_order,
                &recursion_control_points[_RECURSION_CONTROL_POINTS_INDEX(p_der_order,0)]);
    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param, o_der_point);

//    // Освобождение памяти из-под структуры данных для рекурсивного вычисления контрольных точек
//    for (unsigned j = 0; j < p_der_order + 1; ++j)
//        delete[] recursion_control_points[j];
//    delete[] recursion_control_points;

#   undef _
#   undef _RECURSION_CONTROL_POINTS_INDEX
    return 0;
}


unsigned CoxDeBoorAlgorithm::getParamInd(double p_param) const
{
//    // При собирании compressed_control_points значение (param_ind - m_spline_order + 1 + p) должно меняться от 0 до m_control_points_count-1
//    // Это значение меняется
//    // от param_ind - m_spline_order + 1
//    // до param_ind.

//    // Вычисляем минимальный допустимый param_ind, чтобы номер точки не получился меньше нуля:
//    // param_ind - m_spline_order + 1 >= 0 =>
//    // param_ind >= m_spline_order - 1

//    // Вычисляем максимальный допустимый param_ind, чтобы номер точки не получился больше m_control_points_count-1:
//    // param_ind == m_control_points_count-1
//    // т.к. известно, что m_knots_count = m_control_points_count + m_spline_order, то
//    // m_control_points_count = m_knots_count - m_spline_order
//    // значит
//    // param_ind <= m_knots_count - m_spline_order - 1

////    for (unsigned i = m_spline_order-1; i < m_knots_count - m_spline_order; ++i)
////    {
////        if (!(m_knot_vector[i] > p_param) && p_param < m_knot_vector[i + 1])
////            return i;
////    }
//    for (unsigned i = m_spline_order-1; i < m_knots_count - m_spline_order; ++i)
//    {
//        if (p_param < m_knot_vector[i+1])
//            return i;
//    }

//    // Считать p_param == m_knot_vector[m_knots_count - m_spline_order])
//    // не началом отрезка (m_knots_count - m_spline_order), а концом предыдущего (последнего допустимого).
//    return m_knots_count - m_spline_order - 1;

    ///< ALGORITHM A2.1 из книги
    ///  Les Piegl, Wayne Tiller. The NURBS Book. 2nd Edition.

    unsigned n = m_knots_count - m_spline_order - 1;

    // Считать p_param == m_knot_vector[m_knots_count - m_spline_order])
    // не началом отрезка (m_knots_count - m_spline_order), а концом предыдущего (последнего допустимого).
    if (p_param >= m_knot_vector[n + 1])
        return n; /* Special case */
    else if (p_param < m_knot_vector[m_spline_order])
        return m_spline_order - 1;

    unsigned low, high, mid, p;
    p = m_spline_order - 1;
    /* Do binary search */
    low = p;
    high = n + 1;
    mid = (low + high) / 2;
    while ((p_param < m_knot_vector[mid]) || (p_param >= m_knot_vector[mid + 1])) {
        if (p_param < m_knot_vector[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }

    return mid;
}


NFunction::NFunction()
    : m_dim(0)
    , m_order(0)
    , m_knots_count(0)
    , m_knot_vector(nullptr)
    , m_N_array(nullptr)
    , m_left_array(nullptr)
    , m_right_array(nullptr)
{}


NFunction::NFunction(unsigned      p_dim,
                     unsigned      p_spline_order,
                     unsigned      p_knots_count,
                     const double *p_knot_vector)
    : m_dim(p_dim)
    , m_order(p_spline_order)
    , m_knots_count(p_knots_count)
    , m_knot_vector(p_knot_vector)
{
    initializeDependentData();
}


NFunction::~NFunction()
{

}


void NFunction::init(unsigned      p_dim,
                     unsigned      p_order,
                     unsigned      p_knots_count,
                     const double *p_knot_vector)
{
    m_dim                  = p_dim;
    m_order         = p_order;
    m_knots_count          = p_knots_count;
    m_knot_vector          = p_knot_vector;

    initializeDependentData();
}


void NFunction::initializeDependentData()
{
    m_N.resize(m_order);
    m_left.resize(m_order);
    m_right.resize(m_order);
    m_ndu.resize(m_order * m_order);
    m_a.resize(2 * m_order);
    m_ders.resize(2 * m_order);

    cacheArrays();
}


void NFunction::cacheArrays()
{
    m_N_array = m_N.data();
    m_left_array = m_left.data();
    m_right_array = m_right.data();
    m_ndu_array = m_ndu.data();
    m_a_array = m_a.data();
}


void NFunction::clear()
{
    m_dim                  = 0;
    m_order         = 0;
    m_knots_count          = 0;
    m_knot_vector          = nullptr;

    m_N.clear();
    m_left.clear();
    m_right.clear();
    m_ndu.clear();
    m_a.clear();
}


//double NFunction::operator()(unsigned p_k, double p_param) const
//{
////#   define _(point, coord) (m_dim * (point) + (coord))

//    unsigned  param_ind = getParamInd(p_param);

//    unsigned displacement = param_ind - m_spline_order + 1;
//    if ((p_k < displacement) || (p_k > displacement + m_spline_order - 1))
//        return 0;

//    for (unsigned i = 0; i < m_spline_order-1; i++) {
//        m_N[i] = 0;
//    }
//    m_N[m_spline_order-1] = 1;

//    for (unsigned r = 1; r < m_spline_order; ++r) {
//        for (unsigned j = 0; j <= r; ++j) {
//            unsigned i = param_ind + j - r;

//            double new_basis = 0;
//            if (j > 0) {
//                double val_1 = (p_param - m_knot_vector[i]) / (m_knot_vector[i+r] - m_knot_vector[i]);
//                new_basis += val_1*m_N[m_spline_order - 1 + j - r];
//            }
//            if (j < r) {
//                double val_2 = (m_knot_vector[i+r+1] - p_param) / (m_knot_vector[i+r+1] - m_knot_vector[i+1]);
//                new_basis += val_2*m_N[m_spline_order + j - r];
//            }
//            m_N[m_spline_order - 1 + j - r] = new_basis;
//        }
//    }

//    //if ((param_ind + 1 < p_k) || (param_ind > p_k + m_spline_order))
//    //    return 0;

//    //for (unsigned q = 1; q <= m_spline_order; ++q) {
//    //    for (unsigned k = p_k; k <= p_k + m_spline_order - q; ++k) {
//    //        if (q == 1) {
//    //            if ((m_knot_vector[k] <= p_param) && (p_param < m_knot_vector[k+1]) /*&& (m_knot_vector[k] != m_knot_vector[k + 1])*/)
//    //                m_N[k - p_k] = 1.0;
//    //            else if ((k == p_k + m_spline_order - q) &&
//    //                     (p_param == m_knot_vector[k+1]))
//    //                m_N[k - p_k] = 1.0;
//    //            else
//    //                m_N[k - p_k] = 0.0;
//    //        }
//    //        else {
//    //            if (m_knot_vector[k+q-1] != m_knot_vector[k])
//    //                m_N[k - p_k]  = (p_param - m_knot_vector[k]) * m_N[k - p_k]
//    //                        / (m_knot_vector[k+q-1] - m_knot_vector[k]);
//                //else
//                //    m_N[k - p_k] = 0;
//    //            if (m_knot_vector[k+q] != m_knot_vector[k+1])
//    //                m_N[k - p_k] += (m_knot_vector[k+q] - p_param) * m_N[k + 1 - p_k]
//    //                        / (m_knot_vector[k+q] - m_knot_vector[k+1]);
//    //        }
//    //    }
//    //}

////#   undef _
//    return m_N[p_k - displacement];
//    //return m_N[0];
//}


int NFunction::getValues(const double p_param, const double* &o_N, unsigned &o_param_ind) const
{
//    unsigned low, high, mid, p;
//    p = m_spline_order - 1;

    unsigned i = getParamInd(p_param);
    double saved, temp;

    ///< ALGORITHM A2.2 из книги
    ///  Les Piegl, Wayne Tiller. The NURBS Book. 2nd Edition.

    m_N_array[0] = 1.0;

    for (unsigned j = 1; j < m_order; j++) {
        m_left_array[j] = p_param - m_knot_vector[i + 1 - j];
        m_right_array[j] = m_knot_vector[i + j] - p_param;
        saved = 0.0;
        for (unsigned  r = 0; r < j;  r++) {
            temp  = m_N_array[r]/(m_right_array[r+1]+m_left_array[j-r]);
            m_N_array[r] = saved + m_right_array[r + 1] * temp;
            saved = m_left_array[j-r] * temp;
        }
        m_N_array[j] = saved;
    }

    o_N = m_N_array;
    o_param_ind = i;

    return 0;
}

int NFunction::getDerValues(const double p_param, const double *&o_ders, unsigned &o_param_ind) const
{
    unsigned i = getParamInd(p_param);
    unsigned n = m_knots_count - m_order - 1;
    double saved, temp, d;
    unsigned p, s1, s2, rk, pk, j1, j2, j, r;
    p = m_order - 1;

#   define _(point, coord) (m_order * (point) + (coord))

    ///< ALGORITHM A2.3 из книги
    ///  Les Piegl, Wayne Tiller. The NURBS Book. 2nd Edition.

    m_ndu_array[_(0,0)] = 1.0;

    for (j = 1; j < m_order; j++) {
        m_left_array[j] = p_param - m_knot_vector[i + 1 - j];
        m_right_array[j] = m_knot_vector[i + j] - p_param;
        saved = 0.0;
        for (r = 0; r < j;  r++) {
            /* Lower triangle */
            m_ndu_array[_(j,r)] = m_right_array[r+1]+m_left_array[j-r];
            temp = m_ndu_array[_(r,j-1)]/m_ndu_array[_(j,r)];
            /* Upper triangle */
            m_ndu_array[_(r,j)] = saved + m_right_array[r + 1] * temp;
            saved = m_left_array[j-r] * temp;
        }
        m_ndu_array[_(j,j)] = saved;
    }

    for (j = 0; j < m_order; j++) /* Load the basis functioms */
        m_ders_array[_(0,j)] = m_ndu_array[_(j,p)];
    /* This section computes the derivatives (Eq. [2.9]) */
    for (r = 0; r <= p;  r++) { /* Loop over function index */
        s1=0; s2=1; /* Alternate rows in array a */
        m_a_array[_(0,0)] = 1.0;
        /* Loop to compute kth derivative */
        for (unsigned k=1; k<=n; k++) {
            d = 0.0;
            rk = r-k;   pk = p-k;
            if (r >= k) {
                m_a_array[_(s2,0)] = m_a_array[_(s1,0)]/m_ndu_array[_(pk+1,rk)];
                d = m_a_array[_(s2,0)]*m_ndu_array[_(rk,pk)];
            }
            if (rk >= -1)   j1 = 1;
            else            j1 = rk;// Убран Минус критично ли?
            if (r-1 <= pk)  j2 = k-1;
            else            j2 = p-r;
            for (j=j1; j<=j2; j++) {
                m_a_array[_(s2,j)] = (m_a_array[_(s1,j)]-m_a_array[_(s1,j-1)])/m_ndu_array[_(pk+1,rk+j)];
                d += m_a_array[_(s2,j)]*m_ndu_array[_(rk+j,pk)];
            }
            if (r <= pk) {
                m_a_array[_(s2,k)] = -m_a_array[_(s1,k-1)]/m_ndu_array[_(pk+1,r)];
                d += m_a_array[_(s2,k)]*m_ndu_array[_(r,pk)];
            }
            m_ders_array[_(k,r)] = d;
            j=s1;   s1=s2;  s2=j;   /* Switch rows */
        }
    }
    /* Multiply through by the correct factors */
    /* (Eq. [2.9]) */
    r = p;
    for (unsigned k=1; k<=n; k++) {
        for (j=0; j<=p; j++)    m_ders_array[_(k,j)] *= r;
        r *= (p-k);
    }

#   undef _

    o_ders = m_ders_array;
    o_param_ind = i;

    return 0;
}


unsigned NFunction::getParamInd(double p_param) const
{
    ///< ALGORITHM A2.1 из книги

    unsigned n = m_knots_count - m_order - 1;

    if (p_param >= m_knot_vector[n + 1])
        return n; /* Special case */
    else if (p_param < m_knot_vector[m_order])
        return m_order - 1;

    unsigned low, high, mid, p;
    p = m_order - 1;
    /* Do binary search */
    low = p;
    high = n + 1;
    mid = (low + high) / 2;
    while ((p_param < m_knot_vector[mid]) || (p_param >= m_knot_vector[mid + 1])) {
        if (p_param < m_knot_vector[mid])
            high = mid;
        else
            low = mid;
        mid = (low + high) / 2;
    }

    return mid;
}

//unsigned NFunction::getParamInd(double p_param) const
//{
//    for (unsigned i = m_spline_order-1; i < m_knots_count - m_spline_order; ++i)
//    {
//        if (p_param < m_knot_vector[i+1])
//            return i;
//    }

//    return m_knots_count - 1 - m_spline_order;

//}
