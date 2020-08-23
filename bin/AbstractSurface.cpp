#include "AbstractSurface.h"
//#include <BasicMathFunctions>
//#include <unt_StdDefines>
#include <omp.h>


///////////////////////////////////////////////////////////////////////////////
// AbstractSurface
///////////////////////////////////////////////////////////////////////////////

#define POINT_COMPARISON_PRECISION 0.01

AbstractSurface::AbstractSurface(
    unsigned p_id,
    unsigned p_dim,
    double   p_begin_1,
    double   p_end_1,
    double   p_begin_2,
    double   p_end_2,
    bool     p_is_closed_1,
    bool     p_is_closed_2,
    std::map<PoleKeyType, PoleValType> p_poles)
    : m_id(p_id)
    , m_dim(p_dim)
    , m_begin_1(p_begin_1)
    , m_end_1(p_end_1)
    , m_begin_2(p_begin_2)
    , m_end_2(p_end_2)
    , m_is_closed_1(p_is_closed_1)
    , m_is_closed_2(p_is_closed_2)
    , m_poles(p_poles)
{
    initializeDependentData();
}


AbstractSurface::AbstractSurface(
    unsigned p_dim,
    double   p_begin_1,
    double   p_end_1,
    double   p_begin_2,
    double   p_end_2,
    bool     p_is_closed_1,
    bool     p_is_closed_2,
    std::map<PoleKeyType, PoleValType> p_poles)
    : AbstractSurface(
        0,
        p_dim,
        p_begin_1,
        p_end_1,
        p_begin_2,
        p_end_2,
        p_is_closed_1,
        p_is_closed_2,
        p_poles)
{ }


/*AbstractSurface::AbstractSurface(const geo::PrimitiveParameters &p_parameters, unsigned &p_n)
{
    m_id          = std::any_cast<unsigned>(p_parameters[p_n++]);
    m_dim         = std::any_cast<unsigned>(p_parameters[p_n++]);
    m_begin_1     = std::any_cast<double  >(p_parameters[p_n++]);
    m_end_1       = std::any_cast<double  >(p_parameters[p_n++]);
    m_begin_2     = std::any_cast<double  >(p_parameters[p_n++]);
    m_end_2       = std::any_cast<double  >(p_parameters[p_n++]);
    m_is_closed_1 = std::any_cast<bool    >(p_parameters[p_n++]);
    m_is_closed_2 = std::any_cast<bool    >(p_parameters[p_n++]);
    m_poles       = std::any_cast<std::map<PoleKeyType, PoleValType>>(p_parameters[p_n++]);

    initializeDependentData();
}*/


AbstractSurface::~AbstractSurface()
{ }

void AbstractSurface::initializeDependentData()
{
    int threads_count = omp_get_max_threads();
    m_max_der_order = 20;

    m_r1.resize(threads_count * m_dim);
    m_r2.resize(threads_count * m_dim);
    m_counts.resize(threads_count * (m_max_der_order + 1));
    m_counts_base.resize(threads_count * (m_max_der_order + 2));
    m_der_point1.resize(threads_count * m_dim * (m_max_der_order + 3) * (m_max_der_order + 2) / 2);
    m_der_point2.resize(threads_count * m_dim * (m_max_der_order + 3) * (m_max_der_order + 2) / 2);

    m_r1_array = m_r1.data();
    m_r2_array = m_r2.data();
    m_counts_array = m_counts.data();
    m_counts_base_array = m_counts_base.data();
    m_der_point1_array = m_der_point1.data();
    m_der_point2_array = m_der_point2.data();
}


std::pair<double, double> AbstractSurface::getBegin() const
{
    return std::make_pair(m_begin_1, m_begin_2);
}


std::pair<double, double> AbstractSurface::getEnd() const
{
    return std::make_pair(m_end_1, m_end_2);
}

unsigned AbstractSurface::getId() const
{
    return m_id;
}

/*
void AbstractSurface::getNormalPoint(double p_param_1, double p_param_2, double *o_norm_point) const
{
    double *r1 = &m_r1_array[omp_get_thread_num() * m_dim];
    double *r2 = &m_r2_array[omp_get_thread_num() * m_dim];
    double length;

    if (m_poles.find(std::make_pair(1, p_param_1)) != m_poles.end()) {
        getDer1Point(p_param_1, m_poles.at(std::make_pair(1, p_param_1)).first,  r1);
        getDer1Point(p_param_1, m_poles.at(std::make_pair(1, p_param_1)).second, r2);

    } else if (m_poles.find(std::make_pair(2, p_param_2)) != m_poles.end()) {
        getDer2Point(m_poles.at(std::make_pair(1, p_param_2)).first,  p_param_2, r1);
        getDer2Point(m_poles.at(std::make_pair(1, p_param_2)).second, p_param_2, r2);
    } else {
        getDer1Point(p_param_1, p_param_2, r1);
        getDer2Point(p_param_1, p_param_2, r2);
    }

    basicMathFunctions::crossProduct3D(r1, r2, o_norm_point);
    length = basicMathFunctions::length3D(o_norm_point);

    o_norm_point[0] /= length;
    o_norm_point[1] /= length;
    o_norm_point[2] /= length;
}
*/

/*
void AbstractSurface::getDerNormalPoint(unsigned p_der_order,
                                        double   p_param_1,
                                        double   p_param_2,
                                        double  *o_der_norm_point) const
{
    if ((p_der_order > m_max_der_order)) // Требуется увеличить значение m_max_der_order
        abort();

    std::vector<std::vector<unsigned>> bin; // Binomial coefficients C^k_n

    //basicMathFunctions::bci(bin, p_der_order);

    unsigned m;
    double u1;
    double u2;
    double u3;
    double v1;
    double v2;
    double v3;
    double g11;
    double g12;
    double g22;
    double g;
    double b1;
    double b2;
    double b3;

//    double *r1 = &m_r1_array[omp_get_thread_num() * m_dim];
//    double *r2 = &m_r2_array[omp_get_thread_num() * m_dim];

    unsigned *counts = &m_counts_array[omp_get_thread_num() * (m_max_der_order + 1)];
    unsigned *counts_base = &m_counts_base_array[omp_get_thread_num() * (m_max_der_order + 2)];
    
    counts[0] = 0;
    for (unsigned k = 1; k <= p_der_order; k++)
        counts[k] = counts[k - 1] + p_der_order + 2 - k;

    counts_base[0] = 0;
    for (unsigned k = 1; k <= p_der_order+1; k++)
        counts_base[k] = counts_base[k - 1] + p_der_order + 3 - k;

    double *der_point1 = &m_der_point1_array[
        omp_get_thread_num() * m_dim * (m_max_der_order + 3) * (m_max_der_order + 2) / 2];
    double *der_point2 = &m_der_point2_array[
        omp_get_thread_num() * m_dim * (m_max_der_order + 3) * (m_max_der_order + 2) / 2];

    /// Идея вывода этих формул родилась после знакомства с выводом их частных случаев:
    /// см. параграф 1.7 (Вторая квадратичная фомула поверхности,
    ///                   Деривационные формулы  Вейнгартена)
    ///     параграф 1.8 (Третья квадратичная формула поверхности)
    /// Голованов Н.Н. Геометрическое моделирование
    /// М., Изд-во Физматлит, 2002. 472 с.

#   define  _(k, l, coord) ((coord) + m_dim * (counts[k]      + l))
#   define _b(k, l, coord) ((coord) + m_dim * (counts_base[k] + l))
*/
    /*
    if (m_poles.find(std::make_pair(1, p_param_1)) != m_poles.end()) {

        getPointAndDerivs(p_der_order + 1, p_param_1, m_poles.at(std::make_pair(1, p_param_1)).first, der_point1);
        getPointAndDerivs(p_der_order + 1, p_param_1, m_poles.at(std::make_pair(1, p_param_1)).second, der_point2);
        g11 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(1, 0, 0)], &der_point1[_b(1, 0, 0)]);
        g12 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(1, 0, 0)], &der_point2[_b(1, 0, 0)]);
        g22 = basicMathFunctions::dotProduct(m_dim, &der_point2[_b(1, 0, 0)], &der_point2[_b(1, 0, 0)]);
        g = g11 * g22 - g12 * g12;
        getNormalPoint(p_param_1, p_param_2, &o_der_norm_point[_(0, 0, 0)]);
        for (unsigned k = 0; k <= p_der_order; k++)
            for (unsigned l = 0; l <= p_der_order - k; l++)
            {
                u1 = 0;
                u2 = 0;
                u3 = 0;

                for (unsigned i = 0; i <= k; i++)
                {
                    for (unsigned j = 0; j <= l; j++)
                    {
                        if ((i == 0) && (j == 0))
                            continue;

                        u1 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                            &der_point1[_b(i + 1, j, 0)],
                            &o_der_norm_point[_(k - i, l - j, 0)]);
                        u2 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                            &der_point2[_b(i + 1, j, 0)],
                            &o_der_norm_point[_(k - i, l - j, 0)]);

                        if ((i != k) || (j != l))
                            u3 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                                &o_der_norm_point[_(i, j, 0)],
                                &o_der_norm_point[_(k - i, l - j, 0)]);
                    }
                    u1 *= -bin[k][i];
                    u2 *= -bin[k][i];
                    u3 *= -bin[k][i];
                }

                if ((k == 0) && (l == 0))
                    u3 = 1;
                else
                    u3 /= 2.0;

                b1 = (u1 * g22 - u2 * g12) / g;
                b2 = (u2 * g11 - u1 * g12) / g;
                b3 = u3;
                for (m = 0; m < m_dim; ++m)
                    o_der_norm_point[_(k, l, m)] =
                        b1 * der_point1[_b(1, 0, m)] +
                        b2 * der_point2[_b(1, 0, m)] +
                        b3 * o_der_norm_point[_(0, 0, m)];
            }

    } else if (m_poles.find(std::make_pair(2, p_param_2)) != m_poles.end()) {

        getPointAndDerivs(p_der_order + 1, m_poles.at(std::make_pair(1, p_param_2)).first, p_param_2, der_point1);
        getPointAndDerivs(p_der_order + 1, m_poles.at(std::make_pair(1, p_param_2)).second, p_param_2, der_point2);
        g11 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(0, 1, 0)], &der_point1[_b(0, 1, 0)]);
        g12 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(0, 1, 0)], &der_point2[_b(0, 1, 0)]);
        g22 = basicMathFunctions::dotProduct(m_dim, &der_point2[_b(0, 1, 0)], &der_point2[_b(0, 1, 0)]);
        g = g11 * g22 - g12 * g12;
        getNormalPoint(p_param_1, p_param_2, &o_der_norm_point[_(0, 0, 0)]);
        for (unsigned k = 0; k <= p_der_order; k++)
            for (unsigned l = 0; l <= p_der_order - k; l++)
            {
                u1 = 0;
                u2 = 0;
                u3 = 0;

                for (unsigned i = 0; i <= k; i++)
                {
                    for (unsigned j = 0; j <= l; j++)
                    {
                        if ((i == 0) && (j == 0))
                            continue;

                        u1 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                            &der_point1[_b(i, j + 1, 0)],
                            &o_der_norm_point[_(k - i, l - j, 0)]);
                        u2 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                            &der_point2[_b(i, j + 1, 0)],
                            &o_der_norm_point[_(k - i, l - j, 0)]);

                        if ((i != k) || (j != l))
                            u3 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                                &o_der_norm_point[_(i, j, 0)],
                                &o_der_norm_point[_(k - i, l - j, 0)]);
                    }
                    u1 *= -bin[k][i];
                    u2 *= -bin[k][i];
                    u3 *= -bin[k][i];
                }

                if ((k == 0) && (l == 0))
                    u3 = 1;
                else
                    u3 /= 2.0;

                b1 = (u1 * g22 - u2 * g12) / g;
                b2 = (u2 * g11 - u1 * g12) / g;
                b3 = u3;
                for (m = 0; m < m_dim; ++m)
                    o_der_norm_point[_(k, l, m)] =
                        b1 * der_point1[_b(0, 1, m)] +
                        b2 * der_point2[_b(0, 1, m)] +
                        b3 * o_der_norm_point[_(0, 0, m)];
            }

    } else {

        getPointAndDerivs(p_der_order + 1, p_param_1, p_param_2, der_point1);
    //    getDer1Point(p_param_1, p_param_2, r1);
    //    getDer2Point(p_param_1, p_param_2, r2);
        g11 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(1, 0, 0)], &der_point1[_b(1, 0, 0)]);
        g12 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(1, 0, 0)], &der_point1[_b(0, 1, 0)]);
        g22 = basicMathFunctions::dotProduct(m_dim, &der_point1[_b(0, 1, 0)], &der_point1[_b(0, 1, 0)]);
        g = g11 * g22 - g12 * g12;
        getNormalPoint(p_param_1, p_param_2, &o_der_norm_point[_(0, 0, 0)]);
        for (unsigned k = 0; k <= p_der_order; k++)
            for (unsigned l = 0; l <= p_der_order - k; l++)
            {
                u1 = 0;
                u2 = 0;
                u3 = 0;

                for (unsigned i = 0; i <= k; i++)
                {
                    for (unsigned j = 0; j <= l; j++)
                    {
                        if ((i == 0) && (j == 0))
                            continue;

                        u1 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                            &der_point1[_b(i + 1, j, 0)],
                            &o_der_norm_point[_(k - i, l - j, 0)]);
                        u2 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                            &der_point1[_b(i, j + 1, 0)],
                            &o_der_norm_point[_(k - i, l - j, 0)]);

                        if ((i != k) || (j != l))
                            u3 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                                &o_der_norm_point[_(i, j, 0)],
                                &o_der_norm_point[_(k - i, l - j, 0)]);
                    }
                    u1 *= -bin[k][i];
                    u2 *= -bin[k][i];
                    u3 *= -bin[k][i];
                }

                if ((k == 0) && (l == 0))
                    u3 = 1;
                else
                    u3 /= 2.0;

                //for (unsigned j=1; j<=l; j++) {
                //    u1 -= bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                //                                                    &der_point[_b(1, j, 0)],
                //                                                    &o_der_norm_point[_(k, l - j, 0)]);
                //    u2 -= bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                //                                                    &der_point[_b(0, j + 1, 0)],
                //                                                    &o_der_norm_point[_(k, l - j, 0)]);
                //    u3 -= bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                //                                                    &o_der_norm_point[_(0, j, 0)],
                //                                                    &o_der_norm_point[_(k, l-j, 0)]);
                //}

                //for (unsigned i=1; i<=k; i++) {
                //    u1 -= bin[k][i] * basicMathFunctions::dotProduct(m_dim,
                //                                                    &der_point[_b(i + 1, 0, 0)],
                //                                                    &o_der_norm_point[_(k - i, l, 0)]);
                //    u2 -= bin[k][i] * basicMathFunctions::dotProduct(m_dim,
                //                                                    &der_point[_b(i, 1, 0)],
                //                                                    &o_der_norm_point[_(k - i, l, 0)]);
                //    u3 -= bin[k][i] * basicMathFunctions::dotProduct(m_dim,
                //                                                    &o_der_norm_point[_(i ,0, 0)],
                //                                                    &o_der_norm_point[_(k - i, l, 0)]);
                //    v1 = 0;
                //    v2 = 0;
                //    v3 = 0;
                //    for (unsigned j=1; j<=l; j++) {
                //        v1 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                //                                                        &der_point[_b(i+1, j, 0)],
                //                                                        &o_der_norm_point[_(k-i, l-j, 0)]);
                //        v2 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                //                                                        &der_point[_b(i, j+1, 0)],
                //                                                        &o_der_norm_point[_(k-i, 1-j, 0)]);
                //        v3 += bin[l][j] * basicMathFunctions::dotProduct(m_dim,
                //                                                        &o_der_norm_point[_(i,j,0)],
                //                                                        &o_der_norm_point[_(k-i,l-j,0)]);
                //    }
                //    u1 -= bin[k][i] * v1;
                //    u2 -= bin[k][i] * v2;
                //    u3 -= bin[k][i] * v3;
                //}
                b1 = (u1 * g22 - u2 * g12) / g;
                b2 = (u2 * g11 - u1 * g12) / g;
                b3 = u3;
                for (m = 0; m < m_dim; ++m)
                    o_der_norm_point[_(k, l, m)] =
                        b1 * der_point1[_b(1, 0, m)] +
                        b2 * der_point1[_b(0, 1, m)] +
                        b3 * o_der_norm_point[_(0, 0, m)];
            }
    }


    
#   undef _
#   undef _b
}
*/

bool AbstractSurface::isClosed1() const
{
    return  m_is_closed_1;
}


bool AbstractSurface::isClosed2() const
{
    return  m_is_closed_2;
}

/*
bool AbstractSurface::checkClosed1() const
{
#   define POINT_PARAM(number, coord) ((number) * 2 + (coord))

    double period_2 = m_end_2 - m_begin_2;
    double step = period_2 / 4.;

    double test_param_points_1[10] {
        m_begin_1, m_begin_2,
        m_begin_1, m_begin_2 + step,
        m_begin_1, m_begin_2 + step * 2,
        m_begin_1, m_begin_2 + step * 3,
        m_begin_1, m_begin_2 + step * 4
    };

    double test_param_points_2[10] {
        m_end_1, m_begin_2,
        m_end_1, m_begin_2 + step,
        m_end_1, m_begin_2 + step * 2,
        m_end_1, m_begin_2 + step * 3,
        m_end_1, m_begin_2 + step * 4
    };

    double test_3d_points_1[15];
    double test_3d_points_2[15];

    getPoint(test_param_points_1[POINT_PARAM(0, 0)],
             test_param_points_1[POINT_PARAM(0, 1)],
             test_3d_points_1);

    getPoint(test_param_points_1[POINT_PARAM(1, 0)],
             test_param_points_1[POINT_PARAM(1, 1)],
             test_3d_points_1 + 3);

    getPoint(test_param_points_1[POINT_PARAM(2, 0)],
             test_param_points_1[POINT_PARAM(2, 1)],
             test_3d_points_1 + 6);

    getPoint(test_param_points_1[POINT_PARAM(3, 0)],
             test_param_points_1[POINT_PARAM(3, 1)],
             test_3d_points_1 + 9);

    getPoint(test_param_points_1[POINT_PARAM(4, 0)],
             test_param_points_1[POINT_PARAM(4, 1)],
             test_3d_points_1 + 12);

    getPoint(test_param_points_2[POINT_PARAM(0, 0)],
             test_param_points_2[POINT_PARAM(0, 1)],
             test_3d_points_2);

    getPoint(test_param_points_2[POINT_PARAM(1, 0)],
             test_param_points_2[POINT_PARAM(1, 1)],
             test_3d_points_2 + 3);

    getPoint(test_param_points_2[POINT_PARAM(2, 0)],
             test_param_points_2[POINT_PARAM(2, 1)],
             test_3d_points_2 + 6);

    getPoint(test_param_points_2[POINT_PARAM(3, 0)],
             test_param_points_2[POINT_PARAM(3, 1)],
             test_3d_points_2 + 9);

    getPoint(test_param_points_2[POINT_PARAM(4, 0)],
             test_param_points_2[POINT_PARAM(4, 1)],
             test_3d_points_2 + 12);

    double diff_vec[15] {
        (test_3d_points_2     )[0] - (test_3d_points_1     )[0],
        (test_3d_points_2     )[1] - (test_3d_points_1     )[1],
        (test_3d_points_2     )[2] - (test_3d_points_1     )[2],

        (test_3d_points_2 + 3 )[0] - (test_3d_points_1 + 3 )[0],
        (test_3d_points_2 + 3 )[1] - (test_3d_points_1 + 3 )[1],
        (test_3d_points_2 + 3 )[2] - (test_3d_points_1 + 3 )[2],

        (test_3d_points_2 + 6 )[0] - (test_3d_points_1 + 6 )[0],
        (test_3d_points_2 + 6 )[1] - (test_3d_points_1 + 6 )[1],
        (test_3d_points_2 + 6 )[2] - (test_3d_points_1 + 6 )[2],

        (test_3d_points_2 + 9 )[0] - (test_3d_points_1 + 9 )[0],
        (test_3d_points_2 + 9 )[1] - (test_3d_points_1 + 9 )[1],
        (test_3d_points_2 + 9 )[2] - (test_3d_points_1 + 9 )[2],

        (test_3d_points_2 + 12)[0] - (test_3d_points_1 + 12)[0],
        (test_3d_points_2 + 12)[1] - (test_3d_points_1 + 12)[1],
        (test_3d_points_2 + 12)[2] - (test_3d_points_1 + 12)[2]
    };

    double norms_vec[5] {
        basicMathFunctions::length3D(diff_vec     ),
        basicMathFunctions::length3D(diff_vec + 3 ),
        basicMathFunctions::length3D(diff_vec + 6 ),
        basicMathFunctions::length3D(diff_vec + 9 ),
        basicMathFunctions::length3D(diff_vec + 12)
    };

#   undef POINT_PARAM

    return norms_vec[0] < POINT_COMPARISON_PRECISION &&
           norms_vec[1] < POINT_COMPARISON_PRECISION &&
           norms_vec[2] < POINT_COMPARISON_PRECISION &&
           norms_vec[3] < POINT_COMPARISON_PRECISION &&
           norms_vec[4] < POINT_COMPARISON_PRECISION;
}*/

/*
bool AbstractSurface::checkClosed2() const
{
#   define POINT_PARAM(number, coord) ((number) * 2 + (coord))

    double period_1 = m_end_1 - m_begin_1;
    double step = period_1 / 4.;

    double test_param_points_1[10] {
        m_begin_1           , m_begin_2,
        m_begin_1 + step    , m_begin_2,
        m_begin_1 + step * 2, m_begin_2,
        m_begin_1 + step * 3, m_begin_2,
        m_begin_1 + step * 4, m_begin_2
    };

    double test_param_points_2[10] {
        m_begin_1           , m_end_2,
        m_begin_1 + step    , m_end_2,
        m_begin_1 + step * 2, m_end_2,
        m_begin_1 + step * 3, m_end_2,
        m_begin_1 + step * 4, m_end_2
    };

    double test_3d_points_1[15];
    double test_3d_points_2[15];

    getPoint(test_param_points_1[POINT_PARAM(0, 0)],
             test_param_points_1[POINT_PARAM(0, 1)],
             test_3d_points_1);

    getPoint(test_param_points_1[POINT_PARAM(1, 0)],
             test_param_points_1[POINT_PARAM(1, 1)],
             test_3d_points_1 + 3);

    getPoint(test_param_points_1[POINT_PARAM(2, 0)],
             test_param_points_1[POINT_PARAM(2, 1)],
             test_3d_points_1 + 6);

    getPoint(test_param_points_1[POINT_PARAM(3, 0)],
             test_param_points_1[POINT_PARAM(3, 1)],
             test_3d_points_1 + 9);

    getPoint(test_param_points_1[POINT_PARAM(4, 0)],
             test_param_points_1[POINT_PARAM(4, 1)],
             test_3d_points_1 + 12);

    getPoint(test_param_points_2[POINT_PARAM(0, 0)],
             test_param_points_2[POINT_PARAM(0, 1)],
             test_3d_points_2);

    getPoint(test_param_points_2[POINT_PARAM(1, 0)],
             test_param_points_2[POINT_PARAM(1, 1)],
             test_3d_points_2 + 3);

    getPoint(test_param_points_2[POINT_PARAM(2, 0)],
             test_param_points_2[POINT_PARAM(2, 1)],
             test_3d_points_2 + 6);

    getPoint(test_param_points_2[POINT_PARAM(3, 0)],
             test_param_points_2[POINT_PARAM(3, 1)],
             test_3d_points_2 + 9);

    getPoint(test_param_points_2[POINT_PARAM(4, 0)],
             test_param_points_2[POINT_PARAM(4, 1)],
             test_3d_points_2 + 12);

    double diff_vec[15] {
        (test_3d_points_2     )[0] - (test_3d_points_1     )[0],
        (test_3d_points_2     )[1] - (test_3d_points_1     )[1],
        (test_3d_points_2     )[2] - (test_3d_points_1     )[2],

        (test_3d_points_2 + 3 )[0] - (test_3d_points_1 + 3 )[0],
        (test_3d_points_2 + 3 )[1] - (test_3d_points_1 + 3 )[1],
        (test_3d_points_2 + 3 )[2] - (test_3d_points_1 + 3 )[2],

        (test_3d_points_2 + 6 )[0] - (test_3d_points_1 + 6 )[0],
        (test_3d_points_2 + 6 )[1] - (test_3d_points_1 + 6 )[1],
        (test_3d_points_2 + 6 )[2] - (test_3d_points_1 + 6 )[2],

        (test_3d_points_2 + 9 )[0] - (test_3d_points_1 + 9 )[0],
        (test_3d_points_2 + 9 )[1] - (test_3d_points_1 + 9 )[1],
        (test_3d_points_2 + 9 )[2] - (test_3d_points_1 + 9 )[2],

        (test_3d_points_2 + 12)[0] - (test_3d_points_1 + 12)[0],
        (test_3d_points_2 + 12)[1] - (test_3d_points_1 + 12)[1],
        (test_3d_points_2 + 12)[2] - (test_3d_points_1 + 12)[2]
    };

    double norms_vec[5] {
        basicMathFunctions::length3D(diff_vec     ),
        basicMathFunctions::length3D(diff_vec + 3 ),
        basicMathFunctions::length3D(diff_vec + 6 ),
        basicMathFunctions::length3D(diff_vec + 9 ),
        basicMathFunctions::length3D(diff_vec + 12)
    };

#   undef POINT_PARAM

    return norms_vec[0] < POINT_COMPARISON_PRECISION &&
           norms_vec[1] < POINT_COMPARISON_PRECISION &&
           norms_vec[2] < POINT_COMPARISON_PRECISION &&
           norms_vec[3] < POINT_COMPARISON_PRECISION &&
           norms_vec[4] < POINT_COMPARISON_PRECISION;
}
*/

bool AbstractSurface::checkDegenerateBoundary1Begin() const
{
#   define POINT_PARAM(number, coord) ((number) * 2 + (coord))

    double period_2 = m_end_2 - m_begin_2;
    double step = period_2 / 4.;

    double test_param_points_2[10] {
        m_begin_1, m_begin_2,
        m_begin_1, m_begin_2 + step,
        m_begin_1, m_begin_2 + step * 2,
        m_begin_1, m_begin_2 + step * 3,
        m_begin_1, m_begin_2 + step * 4
    };

    double test_3d_points_2[15];

    getPoint(test_param_points_2[POINT_PARAM(0, 0)],
             test_param_points_2[POINT_PARAM(0, 1)],
             test_3d_points_2);

    getPoint(test_param_points_2[POINT_PARAM(1, 0)],
             test_param_points_2[POINT_PARAM(1, 1)],
             test_3d_points_2 + 3);

    getPoint(test_param_points_2[POINT_PARAM(2, 0)],
             test_param_points_2[POINT_PARAM(2, 1)],
             test_3d_points_2 + 6);

    getPoint(test_param_points_2[POINT_PARAM(3, 0)],
             test_param_points_2[POINT_PARAM(3, 1)],
             test_3d_points_2 + 9);

    getPoint(test_param_points_2[POINT_PARAM(4, 0)],
             test_param_points_2[POINT_PARAM(4, 1)],
             test_3d_points_2 + 12);

    double diff_vec[12] {
        test_3d_points_2[0] - (test_3d_points_2 + 3 )[0],
        test_3d_points_2[1] - (test_3d_points_2 + 3 )[1],
        test_3d_points_2[2] - (test_3d_points_2 + 3 )[2],

        test_3d_points_2[0] - (test_3d_points_2 + 6 )[0],
        test_3d_points_2[1] - (test_3d_points_2 + 6 )[1],
        test_3d_points_2[2] - (test_3d_points_2 + 6 )[2],

        test_3d_points_2[0] - (test_3d_points_2 + 9 )[0],
        test_3d_points_2[1] - (test_3d_points_2 + 9 )[1],
        test_3d_points_2[2] - (test_3d_points_2 + 9 )[2],

        test_3d_points_2[0] - (test_3d_points_2 + 12)[0],
        test_3d_points_2[1] - (test_3d_points_2 + 12)[1],
        test_3d_points_2[2] - (test_3d_points_2 + 12)[2]
    };

    double coord_sum = diff_vec[0] + diff_vec[1]  + diff_vec[2] +
                       diff_vec[3] + diff_vec[4]  + diff_vec[5] +
                       diff_vec[6] + diff_vec[7]  + diff_vec[8] +
                       diff_vec[9] + diff_vec[10] + diff_vec[11];

#   undef POINT_PARAM



    return coord_sum < POINT_COMPARISON_PRECISION;
}


bool AbstractSurface::checkDegenerateBoundary1End() const
{
#   define POINT_PARAM(number, coord) ((number) * 2 + (coord))

    double period_2 = m_end_2 - m_begin_2;
    double step = period_2 / 4.;

    double test_param_points_2[10] {
        m_end_1, m_begin_2,
        m_end_1, m_begin_2 + step,
        m_end_1, m_begin_2 + step * 2,
        m_end_1, m_begin_2 + step * 3,
        m_end_1, m_begin_2 + step * 4
    };

    double test_3d_points_2[15];

    getPoint(test_param_points_2[POINT_PARAM(0, 0)],
             test_param_points_2[POINT_PARAM(0, 1)],
             test_3d_points_2);

    getPoint(test_param_points_2[POINT_PARAM(1, 0)],
             test_param_points_2[POINT_PARAM(1, 1)],
             test_3d_points_2 + 3);

    getPoint(test_param_points_2[POINT_PARAM(2, 0)],
             test_param_points_2[POINT_PARAM(2, 1)],
             test_3d_points_2 + 6);

    getPoint(test_param_points_2[POINT_PARAM(3, 0)],
             test_param_points_2[POINT_PARAM(3, 1)],
             test_3d_points_2 + 9);

    getPoint(test_param_points_2[POINT_PARAM(4, 0)],
             test_param_points_2[POINT_PARAM(4, 1)],
             test_3d_points_2 + 12);

    double diff_vec[12] {
        test_3d_points_2[0] - (test_3d_points_2 + 3 )[0],
        test_3d_points_2[1] - (test_3d_points_2 + 3 )[1],
        test_3d_points_2[2] - (test_3d_points_2 + 3 )[2],

        test_3d_points_2[0] - (test_3d_points_2 + 6 )[0],
        test_3d_points_2[1] - (test_3d_points_2 + 6 )[1],
        test_3d_points_2[2] - (test_3d_points_2 + 6 )[2],

        test_3d_points_2[0] - (test_3d_points_2 + 9 )[0],
        test_3d_points_2[1] - (test_3d_points_2 + 9 )[1],
        test_3d_points_2[2] - (test_3d_points_2 + 9 )[2],

        test_3d_points_2[0] - (test_3d_points_2 + 12)[0],
        test_3d_points_2[1] - (test_3d_points_2 + 12)[1],
        test_3d_points_2[2] - (test_3d_points_2 + 12)[2]
    };

    double coord_sum = diff_vec[0] + diff_vec[1]  + diff_vec[2] +
                       diff_vec[3] + diff_vec[4]  + diff_vec[5] +
                       diff_vec[6] + diff_vec[7]  + diff_vec[8] +
                       diff_vec[9] + diff_vec[10] + diff_vec[11];

#   undef POINT_PARAM

    return coord_sum < POINT_COMPARISON_PRECISION;
}


bool AbstractSurface::checkDegenerateBoundary2Begin() const
{
#   define POINT_PARAM(number, coord) ((number) * 2 + (coord))

    double period_1 = m_end_1 - m_begin_1;
    double step = period_1 / 4.;

    double test_param_points_1[10] {
        m_begin_1           , m_begin_2,
        m_begin_1 + step    , m_begin_2,
        m_begin_1 + step * 2, m_begin_2,
        m_begin_1 + step * 3, m_begin_2,
        m_begin_1 + step * 4, m_begin_2
    };

    double test_3d_points_1[15];

    getPoint(test_param_points_1[POINT_PARAM(0, 0)],
             test_param_points_1[POINT_PARAM(0, 1)],
             test_3d_points_1);

    getPoint(test_param_points_1[POINT_PARAM(1, 0)],
             test_param_points_1[POINT_PARAM(1, 1)],
             test_3d_points_1 + 3);

    getPoint(test_param_points_1[POINT_PARAM(2, 0)],
             test_param_points_1[POINT_PARAM(2, 1)],
             test_3d_points_1 + 6);

    getPoint(test_param_points_1[POINT_PARAM(3, 0)],
             test_param_points_1[POINT_PARAM(3, 1)],
             test_3d_points_1 + 9);

    getPoint(test_param_points_1[POINT_PARAM(4, 0)],
             test_param_points_1[POINT_PARAM(4, 1)],
             test_3d_points_1 + 12);

    double diff_vec[12] {
        test_3d_points_1[0] - (test_3d_points_1 + 3 )[0],
        test_3d_points_1[1] - (test_3d_points_1 + 3 )[1],
        test_3d_points_1[2] - (test_3d_points_1 + 3 )[2],

        test_3d_points_1[0] - (test_3d_points_1 + 6 )[0],
        test_3d_points_1[1] - (test_3d_points_1 + 6 )[1],
        test_3d_points_1[2] - (test_3d_points_1 + 6 )[2],

        test_3d_points_1[0] - (test_3d_points_1 + 9 )[0],
        test_3d_points_1[1] - (test_3d_points_1 + 9 )[1],
        test_3d_points_1[2] - (test_3d_points_1 + 9 )[2],

        test_3d_points_1[0] - (test_3d_points_1 + 12)[0],
        test_3d_points_1[1] - (test_3d_points_1 + 12)[1],
        test_3d_points_1[2] - (test_3d_points_1 + 12)[2]
    };

    double coord_sum = diff_vec[0] + diff_vec[1]  + diff_vec[2] +
                       diff_vec[3] + diff_vec[4]  + diff_vec[5] +
                       diff_vec[6] + diff_vec[7]  + diff_vec[8] +
                       diff_vec[9] + diff_vec[10] + diff_vec[11];

#   undef POINT_PARAM

    return coord_sum < POINT_COMPARISON_PRECISION;
}


bool AbstractSurface::checkDegenerateBoundary2End() const
{
#   define POINT_PARAM(number, coord) ((number) * 2 + (coord))

    double period_1 = m_end_1 - m_begin_1;
    double step = period_1 / 4.;

    double test_param_points_1[10] {
        m_begin_1           , m_end_2,
        m_begin_1 + step    , m_end_2,
        m_begin_1 + step * 2, m_end_2,
        m_begin_1 + step * 3, m_end_2,
        m_begin_1 + step * 4, m_end_2
    };

    double test_3d_points_1[15];

    getPoint(test_param_points_1[POINT_PARAM(0, 0)],
             test_param_points_1[POINT_PARAM(0, 1)],
             test_3d_points_1);

    getPoint(test_param_points_1[POINT_PARAM(1, 0)],
             test_param_points_1[POINT_PARAM(1, 1)],
             test_3d_points_1 + 3);

    getPoint(test_param_points_1[POINT_PARAM(2, 0)],
             test_param_points_1[POINT_PARAM(2, 1)],
             test_3d_points_1 + 6);

    getPoint(test_param_points_1[POINT_PARAM(3, 0)],
             test_param_points_1[POINT_PARAM(3, 1)],
             test_3d_points_1 + 9);

    getPoint(test_param_points_1[POINT_PARAM(4, 0)],
             test_param_points_1[POINT_PARAM(4, 1)],
             test_3d_points_1 + 12);

    double diff_vec[12] {
        test_3d_points_1[0] - (test_3d_points_1 + 3 )[0],
        test_3d_points_1[1] - (test_3d_points_1 + 3 )[1],
        test_3d_points_1[2] - (test_3d_points_1 + 3 )[2],

        test_3d_points_1[0] - (test_3d_points_1 + 6 )[0],
        test_3d_points_1[1] - (test_3d_points_1 + 6 )[1],
        test_3d_points_1[2] - (test_3d_points_1 + 6 )[2],

        test_3d_points_1[0] - (test_3d_points_1 + 9 )[0],
        test_3d_points_1[1] - (test_3d_points_1 + 9 )[1],
        test_3d_points_1[2] - (test_3d_points_1 + 9 )[2],

        test_3d_points_1[0] - (test_3d_points_1 + 12)[0],
        test_3d_points_1[1] - (test_3d_points_1 + 12)[1],
        test_3d_points_1[2] - (test_3d_points_1 + 12)[2]
    };

    double coord_sum = diff_vec[0] + diff_vec[1]  + diff_vec[2] +
                       diff_vec[3] + diff_vec[4]  + diff_vec[5] +
                       diff_vec[6] + diff_vec[7]  + diff_vec[8] +
                       diff_vec[9] + diff_vec[10] + diff_vec[11];

#   undef POINT_PARAM

    return coord_sum < POINT_COMPARISON_PRECISION;
}


void AbstractSurface::findInitProjectPoint(
    const double *p_point,
    double       *o_point_parameters) const
{
    unsigned int initial_points_u, initial_points_v;
    initial_points_u = 10;
    initial_points_v = 10;
    double initial_xyz[3];
    double min_length_squared = DBL_MAX;
    for (unsigned i = 0; i <= initial_points_u; ++i)
    {
        double u = static_cast<double>(i) / initial_points_u;
        u = (1 - u) * m_begin_1 + u * m_end_1;
        for (unsigned int j = 0; j <= initial_points_v; ++j)
        {
            double v = double(j)/initial_points_v;
            v = (1 - v) * m_begin_2 + v * m_end_2;
            getPoint(u, v, initial_xyz);
            double length_squared = 0;
            for (int coord = 0; coord < 3; ++coord)
                length_squared += (initial_xyz[coord] - p_point[coord]) *
                                  (initial_xyz[coord] - p_point[coord]);

            if (length_squared < min_length_squared)
            {
                min_length_squared = length_squared;
                o_point_parameters[0] = u;
                o_point_parameters[1] = v;
            }
        }
    }
}


int AbstractSurface::projectPoint(
    const double *p_point,
    double       *o_projected_point,
    double       *o_point_parameters) const
{
    double projected_tangent_u[3];
    double projected_tangent_v[3];
    double projected[3*3];
    unsigned counts[2];

    counts[0] = 0;
    counts[1] = 2;

#   define _(k, l, coord) ((coord) + 3 * (counts[k] + l))

//    getPoint(o_point_parameters[0], o_point_parameters[1], o_projected_point);
    getPointAndDerivs(1, o_point_parameters[0], o_point_parameters[1], projected);
    memcpy(o_projected_point, &projected[_(0, 0, 0)], 3 * sizeof(double));
    memcpy(projected_tangent_u, &projected[_(1, 0, 0)], 3 * sizeof(double));
    memcpy(projected_tangent_v, &projected[_(0, 1, 0)], 3 * sizeof(double));

    double tolerance = 10e-10;
    int iters_threshold = 200;
    int iter = 0;
    double dampfer = 0.5; // 0 < dampfer <= 1

    double proj_u, proj_v, len_u, len_v;
    double vector_to_p[3];

    do {
//        getDer1Point(o_point_parameters[0], o_point_parameters[1], projected_tangent_u);
//        getDer2Point(o_point_parameters[0], o_point_parameters[1], projected_tangent_v);
        len_u = sqrt(projected_tangent_u[0] * projected_tangent_u[0] +
                     projected_tangent_u[1] * projected_tangent_u[1] +
                     projected_tangent_u[2] * projected_tangent_u[2]);
        len_v = sqrt(projected_tangent_v[0] * projected_tangent_v[0] +
                     projected_tangent_v[1] * projected_tangent_v[1] +
                     projected_tangent_v[2] * projected_tangent_v[2]);

        for (int i = 0; i < 3; ++i)
        {
            projected_tangent_u[i] /= len_u;
            projected_tangent_v[i] /= len_v;
        }

        for (int i = 0; i < 3; ++i)
            vector_to_p[i] = p_point[i] - o_projected_point[i];

        proj_u = 0;
        proj_v = 0;
        for (int i = 0; i < 3; ++i) {
            proj_u += vector_to_p[i] * projected_tangent_u[i];
            proj_v += vector_to_p[i] * projected_tangent_v[i];
        }

        if (len_u > tolerance)
        {
            o_point_parameters[0] += proj_u / len_u * dampfer;
            if (m_is_closed_1) {
                if (o_point_parameters[0] < m_begin_1) {
                    o_point_parameters[0] = m_end_1;
                } else if (o_point_parameters[0] > m_end_1) {
                    o_point_parameters[0] = m_begin_1;
                }
            }
        }
        if (len_v > tolerance)
        {
            o_point_parameters[1] += proj_v / len_v * dampfer;
            if (m_is_closed_2) {
                if (o_point_parameters[1] < m_begin_2) {
                    o_point_parameters[1] = m_end_2;
                } else if (o_point_parameters[1] > m_end_2) {
                    o_point_parameters[1] = m_begin_2;
                }
            }
        }
//        getPoint(o_point_parameters[0], o_point_parameters[1], o_projected_point);
        getPointAndDerivs(1, o_point_parameters[0], o_point_parameters[1], projected);
        memcpy( o_projected_point, &projected[_(0, 0, 0)], 3 * sizeof(double) );
        memcpy( projected_tangent_u, &projected[_(1, 0, 0)], 3 * sizeof(double) );
        memcpy( projected_tangent_v, &projected[_(0, 1, 0)], 3 * sizeof(double) );

        iter += 1;
    } while ((fabs(proj_u) > tolerance || fabs(proj_v) > tolerance) && iter < iters_threshold);

    if (iter >= iters_threshold)
        return 1;

    return 0;

#   undef _
}

std::shared_ptr<AbstractSurface> AbstractSurface::reflect2Surface() const
{
    return nullptr;
}


unsigned AbstractSurface::getDim() const
{
    return m_dim;
}


std::shared_ptr<AbstractCurve> AbstractSurface::getBoundary1Begin() const
{
    return m_boundary_begin_1;
}


std::shared_ptr<AbstractCurve> AbstractSurface::getBoundary1End() const
{
    return m_boundary_end_1;
}


std::shared_ptr<AbstractCurve> AbstractSurface::getBoundary2Begin() const
{
    return m_boundary_begin_2;
}


std::shared_ptr<AbstractCurve> AbstractSurface::getBoundary2End() const
{
    return m_boundary_end_2;
}

const std::map<AbstractSurface::PoleKeyType,
               AbstractSurface::PoleValType> AbstractSurface::getPoles() const
{
    return m_poles;
}

void AbstractSurface::setPole(unsigned p_nondegenerated_dir,
                              double p_nondegenerated_param,
                              double p_degenerated_param_1,
                              double p_degenerated_param_2)
{
    m_poles.insert({std::make_pair(p_nondegenerated_dir,p_nondegenerated_param),
                    std::make_pair(p_degenerated_param_1,p_degenerated_param_2)});
}
