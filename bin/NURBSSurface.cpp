#include "NURBSSurface.h"
//#include <BasicMathFunctions>
#include <omp.h>


///////////////////////////////////////////////////////////////////////////////
// NURBSSurface
///////////////////////////////////////////////////////////////////////////////


NURBSSurface::NURBSSurface()
    : m_sum_index_1(0)
    , m_sum_index_2(0)
    , m_degree_1(0)
    , m_degree_2(0)
    , m_knots_1_count(0)
    , m_knots_2_count(0)
    , m_control_points_vector_size(0)
    , m_knots_array_1(nullptr)
    , m_knots_array_2(nullptr)
    , m_weights_array(nullptr)
    , m_control_points_array(nullptr)
{ }


NURBSSurface::NURBSSurface(
    unsigned                   p_id,
    unsigned                   p_dim,
    unsigned                   p_sum_index_1,
    unsigned                   p_sum_index_2,
    unsigned                   p_degree_1,
    unsigned                   p_degree_2,
    const std::vector<double> &p_knot_vector_1,
    const std::vector<double> &p_knot_vector_2,
    const std::vector<double> &p_weights,
    const std::vector<double> &p_control_points,
    double                     p_begin_1,
    double                     p_end_1,
    double                     p_begin_2,
    double                     p_end_2,
    bool                       p_is_closed_1,
    bool                       p_is_closed_2,
    std::map<PoleKeyType, PoleValType> p_poles)
    : AbstractSurface(
          p_id,
          p_dim,
          p_begin_1,
          p_end_1,
          p_begin_2,
          p_end_2,
          p_is_closed_1,
          p_is_closed_2,
          p_poles)
    , m_sum_index_1(p_sum_index_1)
    , m_sum_index_2(p_sum_index_2)
    , m_degree_1(p_degree_1)
    , m_degree_2(p_degree_2)
    , m_knot_vector_1(p_knot_vector_1)
    , m_knot_vector_2(p_knot_vector_2)
    , m_weights(p_weights)
    , m_control_points(p_control_points)
{
    initializeDependentData();
}


NURBSSurface::NURBSSurface(
    unsigned                   p_dim,
    unsigned                   p_sum_index_1,
    unsigned                   p_sum_index_2,
    unsigned                   p_degree_1,
    unsigned                   p_degree_2,
    const std::vector<double> &p_knot_vector_1,
    const std::vector<double> &p_knot_vector_2,
    const std::vector<double> &p_weights,
    const std::vector<double> &p_control_points,
    double                     p_begin_1,
    double                     p_end_1,
    double                     p_begin_2,
    double                     p_end_2,
    bool                       p_is_closed_1,
    bool                       p_is_closed_2,
    std::map<PoleKeyType, PoleValType> p_poles)
    : NURBSSurface(
          0,
          p_dim,
          p_sum_index_1,
          p_sum_index_2,
          p_degree_1,
          p_degree_2,
          p_knot_vector_1,
          p_knot_vector_2,
          p_weights,
          p_control_points,
          p_begin_1,
          p_end_1,
          p_begin_2,
          p_end_2,
          p_is_closed_1,
          p_is_closed_2,
          p_poles)
{ }


NURBSSurface::NURBSSurface(const NURBSSurface &p_nurbs_surface)
    : AbstractSurface(p_nurbs_surface)
    , m_sum_index_1(p_nurbs_surface.m_sum_index_1)
    , m_sum_index_2(p_nurbs_surface.m_sum_index_2)
    , m_degree_1(p_nurbs_surface.m_degree_1)
    , m_degree_2(p_nurbs_surface.m_degree_2)
    , m_knot_vector_1(p_nurbs_surface.m_knot_vector_1)
    , m_knot_vector_2(p_nurbs_surface.m_knot_vector_2)
    , m_weights(p_nurbs_surface.m_weights)
    , m_control_points(p_nurbs_surface.m_control_points)
{
    initializeDependentData();
}


/*NURBSSurface::NURBSSurface(const geo::PrimitiveParameters &p_parameters, unsigned &p_n)
    :AbstractSurface(p_parameters, p_n)
{
    m_degree_1       = std::any_cast<unsigned           >(p_parameters[p_n++]);
    m_degree_2       = std::any_cast<unsigned           >(p_parameters[p_n++]);
    m_weights        = std::any_cast<std::vector<double>>(p_parameters[p_n++]);
    m_control_points = std::any_cast<std::vector<double>>(p_parameters[p_n++]);
    m_knot_vector_1  = std::any_cast<std::vector<double>>(p_parameters[p_n++]);
    m_knot_vector_2  = std::any_cast<std::vector<double>>(p_parameters[p_n++]);

    m_sum_index_1 = m_knot_vector_1.size() - m_degree_1 - 2;
    m_sum_index_2 = m_knot_vector_2.size() - m_degree_2 - 2;

    initializeDependentData();
}*/


const NURBSSurface& NURBSSurface::operator=(const NURBSSurface &p_nurbs_surface)
{
    if (this == &p_nurbs_surface)
        return *this;

    AbstractSurface::operator=(p_nurbs_surface);

    m_sum_index_1 = p_nurbs_surface.m_sum_index_1;
    m_sum_index_2 = p_nurbs_surface.m_sum_index_2;
    m_degree_1 = p_nurbs_surface.m_degree_1;
    m_degree_2 = p_nurbs_surface.m_degree_2;
    m_knot_vector_1 = p_nurbs_surface.m_knot_vector_1;
    m_knot_vector_2 = p_nurbs_surface.m_knot_vector_2;
    m_weights = p_nurbs_surface.m_weights;
    m_control_points = p_nurbs_surface.m_control_points;

    initializeDependentData();

    return *this;
}


void NURBSSurface::initializeDependentData()
{
    // Кэширование данных
    m_order_1 = m_degree_1 + 1;
    m_order_2 = m_degree_2 + 1;
    m_knots_1_count = static_cast<unsigned>(m_knot_vector_1.size());
    m_knots_2_count = static_cast<unsigned>(m_knot_vector_2.size());
    m_control_points_1_count = m_sum_index_1 + 1;
    m_control_points_2_count = m_sum_index_2 + 1;
    m_control_points_vector_size = static_cast<unsigned>(m_control_points.size());

    m_weights_transpose.resize(
        m_control_points_1_count * m_control_points_2_count);
    m_weighted_control_points.resize(
        m_control_points_1_count * m_control_points_2_count * (m_dim + 1));
    m_weighted_control_points_transpose.resize(
        m_control_points_1_count*m_control_points_2_count * (m_dim + 1));

    int threads_count = omp_get_max_threads();
    m_cox_de_boor_algorithm.resize(threads_count);
    m_nurbs_numerator_point.resize(threads_count * (m_dim + 1));
    m_nurbs_numerator_der_point.resize(threads_count * (m_dim + 1));
    m_v.resize(threads_count * (m_dim + 1));
    m_v2.resize(threads_count * m_dim);
    m_control_points_1.resize(threads_count * (m_dim + 1) * m_control_points_1_count);
    m_control_points_2.resize(threads_count * (m_dim + 1) * m_control_points_2_count);
    m_weights_1.resize(threads_count * m_control_points_1_count);
    m_weights_2.resize(threads_count * m_control_points_2_count);
    m_Sw.resize(threads_count * (m_dim + 1));

    cacheArrays();

    std::vector<unsigned> t(4);
    t[0] = 1;
    m_bin.push_back(t);
    t[1] = 1;
    m_bin.push_back(t);
    t[1] = 2;
    t[2] = 1;
    m_bin.push_back(t);
    t[1] = 3;
    t[2] = 3;
    t[3] = 1;
    m_bin.push_back(t);
    t.~vector();
    //basicMathFunctions::bci(m_bin, m_degree_1 + m_degree_2);

    transpose(m_weights_array,
              1,
              m_control_points_1_count,
              m_control_points_2_count,
              m_weights_array_transpose);
    getWeightedControlPoints(m_weighted_control_points_array);
    getWeightedControlPointsTranspose(m_weighted_control_points_transpose_array);

    m_numerator_cox_de_boor_algorithm1.resize(m_control_points_2_count);
    m_denominator_cox_de_boor_algorithm1.resize(m_control_points_2_count);

    for (unsigned j = 0; j < m_control_points_2_count; ++j)
    {
        m_numerator_cox_de_boor_algorithm1[j] = std::make_unique<CoxDeBoorAlgorithm> (
            m_dim + 1,
            m_order_1,
            m_knots_1_count,
            m_knots_array_1,
            m_control_points_1_count,
            m_weighted_control_points_array + (m_dim + 1) * m_control_points_1_count * j);

        m_denominator_cox_de_boor_algorithm1[j] = std::make_unique<CoxDeBoorAlgorithm> (
            1,
            m_order_1,
            m_knots_1_count,
            m_knots_array_1,
            m_control_points_1_count,
            m_weights_array + m_control_points_1_count * j);
    }

    m_numerator_cox_de_boor_algorithm2.resize(m_control_points_1_count);
    m_denominator_cox_de_boor_algorithm2.resize(m_control_points_1_count);

    for (unsigned i = 0; i < m_control_points_1_count; ++i)
    {
        m_numerator_cox_de_boor_algorithm2[i] = std::make_unique<CoxDeBoorAlgorithm> (
            m_dim + 1,
            m_order_2,
            m_knots_2_count,
            m_knots_array_2,
            m_control_points_2_count,
            m_weighted_control_points_transpose_array + (m_dim + 1) * m_control_points_2_count * i);

        m_denominator_cox_de_boor_algorithm2[i] = std::make_unique<CoxDeBoorAlgorithm> (
            1,
            m_order_2,
            m_knots_2_count,
            m_knots_array_2,
            m_control_points_2_count,
            m_weights_array_transpose + m_control_points_2_count * i);
    }
}


void NURBSSurface::cacheArrays()
{
    m_knots_array_1 = m_knot_vector_1.data();
    m_knots_array_2 = m_knot_vector_2.data();
    m_weights_array = m_weights.data();
    m_control_points_array = m_control_points.data();

    m_cox_de_boor_algorithm_array = m_cox_de_boor_algorithm.data();
    m_nurbs_numerator_point_array = m_nurbs_numerator_point.data();
    m_nurbs_numerator_der_point_array = m_nurbs_numerator_der_point.data();
    m_v_array = m_v.data();
    m_v2_array = m_v2.data();
    m_control_points_1_array = m_control_points_1.data();
    m_control_points_2_array = m_control_points_2.data();
    m_weights_1_array = m_weights_1.data();
    m_weights_2_array = m_weights_2.data();
    m_weights_array_transpose = m_weights_transpose.data();
    m_weighted_control_points_array = m_weighted_control_points.data();
    m_weighted_control_points_transpose_array = m_weighted_control_points_transpose.data();
    m_Sw_array = m_Sw.data();
}


void NURBSSurface::getPoint(
    double  p_param_1,
    double  p_param_2,
    double *o_point) const
{
//    double denominator_value;
    double period_1 = std::fabs(m_end_1 - m_begin_1);
    double period_2 = std::fabs(m_end_2 - m_begin_2);
    double *Sw = &m_Sw_array[omp_get_thread_num() * (m_dim + 1)];

    // Циклический пересчет параметров для закнутых поверхностей
    double param_1 =
        m_is_closed_1 ?
        (
            (p_param_1 < m_begin_1) ?
                p_param_1 +
                static_cast<unsigned>(std::fabs(m_end_1 - p_param_1) /
                period_1) *
                period_1
            : (
                (p_param_1 > m_end_1) ?
                    p_param_1 -
                    static_cast<unsigned>(std::fabs(p_param_1 - m_begin_1) /
                    period_1) *
                    period_1
                : p_param_1
            )
        ) :
        p_param_1;

    double param_2 =
        m_is_closed_2 ?
        (
            (p_param_2 < m_begin_2) ?
                p_param_2 +
                static_cast<unsigned>(std::fabs(m_end_2 - p_param_2) /
                period_2) *
                period_2
            : (
                (p_param_2 > m_end_2) ?
                    p_param_2 -
                    static_cast<unsigned>(std::fabs(p_param_2 - m_begin_2) /
                    period_2) *
                    period_2
                : p_param_2
            )
        ) :
        p_param_2;

    getNumerator(param_1, param_2, Sw);
//    getDenominator(param_1, param_2, denominator_value);

    for (unsigned d = 0; d < m_dim; ++d)
        o_point[d] = Sw[d] / Sw[m_dim];
}


void NURBSSurface::getDer1Point(
    double  p_param_1,
    double  p_param_2,
    double *o_der_1_point) const
{
    double *nurbs_numerator_point = &m_nurbs_numerator_point_array[
        omp_get_thread_num() * (m_dim + 1)];
    double *nurbs_numerator_der1_point = &m_nurbs_numerator_der_point_array[
        omp_get_thread_num() * (m_dim + 1)];
//    double  nurbs_denominator_value;
//    double  nurbs_denominator_der1_value;

    double period_1 = std::fabs(m_end_1 - m_begin_1);
    double period_2 = std::fabs(m_end_2 - m_begin_2);

    // Циклический пересчет параметров для закнутых поверхностей
    double param_1 =
        m_is_closed_1 ?
        (
            (p_param_1 < m_begin_1) ?
                p_param_1 +
                static_cast<unsigned>(std::fabs(m_end_1 - p_param_1) /
                period_1) *
                period_1
            : (
                (p_param_1 > m_end_1) ?
                    p_param_1 -
                    static_cast<unsigned>(std::fabs(p_param_1 - m_begin_1) /
                    period_1) *
                    period_1
                : p_param_1
            )
        ) :
        p_param_1;

    double param_2 =
        m_is_closed_2 ?
        (
            (p_param_2 < m_begin_2) ?
                p_param_2 +
                static_cast<unsigned>(std::fabs(m_end_2 - p_param_2) /
                period_2) *
                period_2
            : (
                (p_param_2 > m_end_2) ?
                    p_param_2 -
                    static_cast<unsigned>(std::fabs(p_param_2 - m_begin_2) /
                    period_2) *
                    period_2
                : p_param_2
            )
        ) :
        p_param_2;

    getNumerator(param_1, param_2, nurbs_numerator_point);
    getDer1Numerator(param_1, param_2, nurbs_numerator_der1_point);
//    getDenominator(param_1, param_2, nurbs_denominator_value);
//    getDer1Denominator(param_1, param_2, nurbs_denominator_der1_value);

//    for (unsigned d = 0; d < m_dim; ++d)
//        o_der_1_point[d] = (nurbs_numerator_der1_point[d] * nurbs_denominator_value -
//                            nurbs_numerator_point[d] * nurbs_denominator_der1_value) /
//                           (nurbs_denominator_value * nurbs_denominator_value);
    for (unsigned d = 0; d < m_dim; ++d)
        o_der_1_point[d] = (nurbs_numerator_der1_point[d] * nurbs_numerator_point[m_dim] -
                            nurbs_numerator_point[d] * nurbs_numerator_der1_point[m_dim]) /
                           (nurbs_numerator_point[m_dim] * nurbs_numerator_point[m_dim]);
}


void NURBSSurface::getDer2Point(
    double  p_param_1,
    double  p_param_2,
    double *o_der_2_point) const
{
    double *nurbs_numerator_point = &m_nurbs_numerator_point_array[
        omp_get_thread_num() * (m_dim + 1)];
    double *nurbs_numerator_der2_point = &m_nurbs_numerator_der_point_array[
        omp_get_thread_num() * (m_dim + 1)];
//    double  nurbs_denominator_value;
//    double  nurbs_denominator_der2_value;

    double period_1 = std::fabs(m_end_1 - m_begin_1);
    double period_2 = std::fabs(m_end_2 - m_begin_2);

    // Циклический пересчет параметров для замкнутых поверхностей
    double param_1 =
        m_is_closed_1 ?
        (
            (p_param_1 < m_begin_1) ?
                p_param_1 +
                static_cast<unsigned>(std::fabs(m_end_1 - p_param_1) /
                period_1) *
                period_1
            : (
                (p_param_1 > m_end_1) ?
                    p_param_1 -
                    static_cast<unsigned>(std::fabs(p_param_1 - m_begin_1) /
                    period_1) *
                    period_1
                : p_param_1
            )
        ) :
        p_param_1;

    double param_2 =
        m_is_closed_2 ?
        (
            (p_param_2 < m_begin_2) ?
                p_param_2 +
                static_cast<unsigned>(std::fabs(m_end_2 - p_param_2) /
                period_2) *
                period_2
            : (
                (p_param_2 > m_end_2) ?
                    p_param_2 -
                    static_cast<unsigned>(std::fabs(p_param_2 - m_begin_2) /
                    period_2) *
                    period_2
                : p_param_2
            )
        ) :
        p_param_2;

    getNumerator(param_1, param_2, nurbs_numerator_point);
    getDer2Numerator(param_1, param_2, nurbs_numerator_der2_point);
//    getDenominator(param_1, param_2, nurbs_denominator_value);
//    getDer2Denominator(param_1, param_2, nurbs_denominator_der2_value);

//    for (unsigned d = 0; d < m_dim; ++d)
//        o_der_2_point[d] = (nurbs_numerator_der2_point[d] * nurbs_denominator_value -
//                            nurbs_numerator_point[d] * nurbs_denominator_der2_value) /
//                           (nurbs_denominator_value * nurbs_denominator_value);
    for (unsigned d = 0; d < m_dim; ++d)
        o_der_2_point[d] = (nurbs_numerator_der2_point[d] * nurbs_numerator_point[m_dim] -
                            nurbs_numerator_point[d] * nurbs_numerator_der2_point[m_dim]) /
                           (nurbs_numerator_point[m_dim] * nurbs_numerator_point[m_dim]);
}


void NURBSSurface::getPointAndDerivs(unsigned p_der_order,
                                     double p_param_1,
                                     double p_param_2,
                                     double *o_der_point) const
{
    if ((p_der_order > m_max_der_order)) // Требуется увеличить значение m_max_der_order
        abort();

        double period_1 = std::fabs(m_end_1 - m_begin_1);
    double period_2 = std::fabs(m_end_2 - m_begin_2);

    // Циклический пересчет параметров для замкнутых поверхностей
    p_param_1 =
        m_is_closed_1 ?
        (
            (p_param_1 < m_begin_1) ?
                p_param_1 +
                static_cast<unsigned>(std::fabs(m_end_1 - p_param_1) /
                period_1) *
                period_1
            : (
                (p_param_1 > m_end_1) ?
                    p_param_1 -
                    static_cast<unsigned>(std::fabs(p_param_1 - m_begin_1) /
                    period_1) *
                    period_1
                : p_param_1
            )
        ) :
        p_param_1;

    p_param_2 =
        m_is_closed_2 ?
        (
            (p_param_2 < m_begin_2) ?
                p_param_2 +
                static_cast<unsigned>(std::fabs(m_end_2 - p_param_2) /
                period_2) *
                period_2
            : (
                (p_param_2 > m_end_2) ?
                    p_param_2 -
                    static_cast<unsigned>(std::fabs(p_param_2 - m_begin_2) /
                    period_2) *
                    period_2
                : p_param_2
            )
        ) :
        p_param_2;

    unsigned m;

//    double *control_points_2 = new double[m_dim * m_control_points_2_count];
//    double *weighed_control_points = new double[m_control_points_vector_size];

    double *v = &m_v_array[omp_get_thread_num() * (m_dim + 1)];
    double *v2 = &m_v2_array[omp_get_thread_num() * m_dim];
    unsigned *counts = &m_counts_array[omp_get_thread_num() * (m_max_der_order + 1)];
    double wders0;
    double wders;

    counts[0] = 0;
    for (unsigned k = 1; k <= p_der_order; k++)
        counts[k] = counts[k - 1] + p_der_order + 2 - k;

#   define _(k, l, coord) ((coord) + m_dim * (counts[k] + l))
    for (unsigned k = m_degree_1 + m_degree_2 + 1; k <= p_der_order; k++)
        for (unsigned l = 0; l <= p_der_order - k; l++)
            for (m = 0; m < m_dim; ++m)
                o_der_point[_(k, l, m)] = 0;

    for (unsigned l = m_degree_1 + m_degree_2 + 1; l <= p_der_order; l++)
        for (unsigned k = 0; k <= p_der_order - l; k++)
            for (m = 0; m < m_dim; ++m)
                o_der_point[_(k, l, m)] = 0;

    p_der_order = std::min(p_der_order, m_degree_1 + m_degree_2);

    /// Расчет проводится по формулам из книги:
    /// Les Piegl, Wayne Tiller. The NURBS Book. 2nd Edition. Paragraph 4.5
    /// Русский перевод:
    /// Ю.Д. Шевелев. Математические основы задач проектирования:
    /// Учебное пособие по курсу "Математическое и программное обеспечение САПР"
    /// М., 2005. Параграф 2.11

    // ALGORITHM A4.4
    getDenominator(p_param_1, p_param_2, wders0);
    for (unsigned k = 0; k <= p_der_order; k++)
        for (unsigned l = 0; l <= p_der_order - k; l++)
        {
            getDerNumerator(k, l, p_param_1, p_param_2, v);
            for (unsigned j = 1; j <= l; j++)
            {
                getDerDenominator(0, j, p_param_1, p_param_2, wders);
                for (m = 0; m < m_dim; ++m)
                    v[m] = v[m] - m_bin[l][j] * wders * o_der_point[_(k, l - j, m)];//на втором проходе цикла
            }

            for (unsigned i = 1; i <= k; i++)
            {
                getDerDenominator(i, 0, p_param_1, p_param_2, wders);
                for (m = 0; m < m_dim; ++m)
                    v[m] = v[m] - m_bin[k][i] * wders * o_der_point[_(k - i, l, m)];

                for (m = 0; m < m_dim; ++m)
                    v2[m] = 0;

                for (unsigned j = 1; j <= l; j++)
                {
                    getDerDenominator(i, j, p_param_1, p_param_2, wders);
                    for (m = 0; m < m_dim; ++m)
                        v2[m] = v2[m] + m_bin[l][j] * wders * o_der_point[_(k - i, l - j, m)];
                }

                for (m = 0; m < m_dim; ++m)
                    v[m] = v[m] - m_bin[k][i] * v2[m];
            }

            for (m = 0; m < m_dim; ++m)
                o_der_point[_(k, l, m)] = v[m] / wders0;
    }
#   undef _
}


void NURBSSurface::findInitProjectPoint(
    const double *p_point,
    double       *o_point_parameters) const
{
    unsigned int m_initial_points_u, m_initial_points_v;
    m_initial_points_u = m_control_points_1_count;
    m_initial_points_v = m_control_points_2_count;
    double initial_xyz[3];
    double min_length_squared = DBL_MAX;

    for (unsigned i = 0; i <= m_initial_points_u; ++i)
    {
        double u = double(i)/m_initial_points_u;
        u = (1 - u) * m_begin_1 + u * m_end_1;
        for (unsigned j = 0; j <= m_initial_points_v; ++j)
        {
            double v = static_cast<double>(j) / m_initial_points_v;
            v = (1 - v) * m_begin_2 + v * m_end_2;
            getPoint(u, v, initial_xyz);
            double length_squared = 0;

            for (unsigned coord = 0; coord < 3; ++coord)
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


std::shared_ptr<AbstractSurface> NURBSSurface::reflect2Surface() const
{
    std::vector<double> control_points(m_control_points_vector_size);
    std::vector<double> weights(m_control_points_vector_size);

#   define _WEIGHTS_INDEX(i, j) ((i) + m_control_points_1_count * (j))
#   define _CONTROL_POINTS_INDEX(i, j, coord) ((coord) + m_dim * ((i) + m_control_points_1_count * (j)))

    for (unsigned i = 0; i < m_control_points_1_count; ++i)
        for (unsigned j = 0; j < m_control_points_2_count; ++j)
            for (unsigned d = 0; d < m_dim; ++d)
            {
                control_points[_CONTROL_POINTS_INDEX(m_control_points_1_count - 1 - i, j, d)] =
                    m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
                weights[_WEIGHTS_INDEX(m_control_points_1_count - 1 - i, j)] =
                        m_weights_array[_WEIGHTS_INDEX(i, j)];
            }

#   undef _WEIGHTS_INDEX
#   undef _CONTROL_POINTS_INDEX

    return std::make_shared<NURBSSurface>(
        m_id,
        m_dim,
        m_sum_index_1,
        m_sum_index_2,
        m_degree_1,
        m_degree_2,
        m_knot_vector_1,
        m_knot_vector_2,
        weights,
        control_points,
        m_begin_1,
        m_end_1,
        m_begin_2,
        m_end_2,
        m_is_closed_1,
        m_is_closed_2);
}


void NURBSSurface::setVariableData(
        unsigned                   p_degree_1,
        unsigned                   p_degree_2,
        const std::vector<double> &p_knot_vector_1,
        const std::vector<double> &p_knot_vector_2,
        const std::vector<double> &p_control_points)
{
    m_degree_1 = p_degree_1;
    m_degree_2 = p_degree_2;
    m_knot_vector_1 = p_knot_vector_1;
    m_knot_vector_2 = p_knot_vector_2;
    m_control_points = p_control_points;

    initializeDependentData();
}

void NURBSSurface::setControlPoints(const std::vector<double> &p_control_points)
{
    m_control_points = p_control_points;
    m_control_points_array = m_control_points.data();

    initializeDependentData();
}


unsigned NURBSSurface::degree1()
{
    return m_degree_1;
}


unsigned NURBSSurface::degree2()
{
    return m_degree_2;
}


const std::vector<double> &NURBSSurface::weights()
{
    return m_weights;
}


const std::vector<double> &NURBSSurface::controlPoints()
{
    return m_control_points;
}


const std::vector<double> &NURBSSurface::knotVector1()
{
    return m_knot_vector_1;
}


const std::vector<double> &NURBSSurface::knotVector2()
{
    return m_knot_vector_2;
}


void NURBSSurface::calculateBoundaryCurves()
{
    std::vector<double> control_points_begin_1(m_dim * m_control_points_2_count);
    std::vector<double> control_points_end_1(m_dim * m_control_points_2_count);
    std::vector<double> control_points_begin_2(m_dim * m_control_points_1_count);
    std::vector<double> control_points_end_2(m_dim * m_control_points_1_count);
    std::vector<double> weights_begin_1(m_control_points_2_count);
    std::vector<double> weights_end_1(m_control_points_2_count);
    std::vector<double> weights_begin_2(m_control_points_1_count);
    std::vector<double> weights_end_2(m_control_points_1_count);

#   define _WEIGHTS_INDEX(i, j) ((i) + m_control_points_1_count * (j))
#   define _CONTROL_POINTS_INDEX(i, j, coord) ((coord) + m_dim*((i) + m_control_points_1_count*(j)))
#   define _(point, coord) (m_dim * (point) + (coord))


    for (unsigned i = 0; i < m_control_points_1_count; ++i)
        for (unsigned j = 0; j < m_control_points_2_count; ++j)
            for (unsigned d = 0; d < m_dim; ++d)
            {
                if (i == 0) // left boundary
                {
                    control_points_begin_1[_(j, d)] = m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
                    weights_begin_1[j] = m_weights_array[_WEIGHTS_INDEX(i, j)];
                }
                if (i == m_control_points_1_count - 1) // right boundary
                {
                    control_points_end_1[_(j, d)] = m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
                    weights_end_1[j] = m_weights_array[_WEIGHTS_INDEX(i, j)];
                }
                if (j == 0) // bottom boundary
                {
                    control_points_begin_2[_(i, d)] = m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
                    weights_begin_2[i] = m_weights_array[_WEIGHTS_INDEX(i, j)];
                }
                if (j == m_control_points_2_count - 1) // bottom boundary
                {
                    control_points_end_2[_(i, d)] = m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
                    weights_end_2[i] = m_weights_array[_WEIGHTS_INDEX(i, j)];
                }
            }

#   undef _WEIGHTS_INDEX
#   undef _CONTROL_POINTS_INDEX
#   undef _

    m_boundary_begin_1 = std::make_shared<NURBSCurve>(
        m_dim,
        m_degree_2,
        m_knot_vector_2,
        weights_begin_1,
        control_points_begin_1,
        std::vector<double> (),
        m_begin_2,
        m_end_2,
        m_is_closed_2);

    m_boundary_end_1 = std::make_shared<NURBSCurve>(
        m_dim,
        m_degree_2,
        m_knot_vector_2,
        weights_end_1,
        control_points_end_1,
        std::vector<double> (),
        m_begin_2,
        m_end_2,
        m_is_closed_2);

    m_boundary_begin_2 = std::make_shared<NURBSCurve>(
        m_dim,
        m_degree_1,
        m_knot_vector_1,
        weights_begin_2,
        control_points_begin_2,
        std::vector<double> (),
        m_begin_1,
        m_end_1,
        m_is_closed_1);

    m_boundary_end_2 = std::make_shared<NURBSCurve>(
        m_dim,
        m_degree_1,
        m_knot_vector_1,
        weights_end_2,
        control_points_end_2,
        std::vector<double> (),
        m_begin_1,
        m_end_1,
        m_is_closed_1);
}


int NURBSSurface::getNumerator(double p_param_1, double p_param_2, double *o_point) const
{
    double *control_points_2 = &m_control_points_2_array[
        omp_get_thread_num() * (m_dim + 1) * m_control_points_2_count];
//    double *weighed_control_points = &m_weighed_control_points_array[omp_get_thread_num() * m_control_points_vector_size];

//    getWeighedControlPoints(weighed_control_points);

//    CoxDeBoorAlgorithm cox_de_boor_algorithm;

    for (unsigned j = 0; j < m_control_points_2_count; ++j)
    {
//        cox_de_boor_algorithm.init(
//            m_dim,
//            m_order_1,
//            m_knots_1_count,
//            m_knots_array_1,
//            m_control_points_1_count,
//            weighed_control_points + m_dim * m_control_points_1_count * j);
//        cox_de_boor_algorithm.getPoint(p_param_1, control_points_2 + m_dim * j);
//        cox_de_boor_algorithm.clear();

        m_numerator_cox_de_boor_algorithm1[j]->getPoint(p_param_1, control_points_2 + (m_dim + 1) * j);
    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
        m_dim + 1,
        m_order_2,
        m_knots_2_count,
        m_knots_array_2,
        m_control_points_2_count,
        control_points_2);
    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_2, o_point);

    return 0;
}


int NURBSSurface::getDenominator(double p_param_1, double p_param_2, double &o_value) const
{
    double *control_points_2 = &m_weights_2_array[omp_get_thread_num() * m_control_points_2_count];

//    CoxDeBoorAlgorithm cox_de_boor_algorithm;
    

    for (unsigned j = 0; j < m_control_points_2_count; ++j)
    {
//        cox_de_boor_algorithm.init(
//            1,
//            m_order_1,
//            m_knots_1_count,
//            m_knots_array_1,
//            m_control_points_1_count,
//            m_weights_array + m_control_points_1_count * j);
//        cox_de_boor_algorithm.getPoint(p_param_1, control_points_2 + j);
//        cox_de_boor_algorithm.clear();

        m_denominator_cox_de_boor_algorithm1[j]->getPoint(p_param_1, control_points_2 + j);

    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
        1,
        m_order_2,
        m_knots_2_count,
        m_knots_array_2,
        m_control_points_2_count,
        control_points_2);
    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_2, &o_value);

    return 0;
}


int NURBSSurface::getDer1Numerator(double p_param_1, double p_param_2, double *o_point) const
{
    double *control_points_2 = &m_control_points_2_array[
        omp_get_thread_num() * (m_dim + 1) * m_control_points_2_count];
//    double *weighed_control_points = &m_weighed_control_points_array[omp_get_thread_num() * m_control_points_vector_size];

//    getWeighedControlPoints(weighed_control_points);

//    CoxDeBoorAlgorithm cox_de_boor_algorithm;

    for (unsigned j = 0; j < m_control_points_2_count; ++j)
    {
//        cox_de_boor_algorithm.init(
//            m_dim,
//            m_order_1,
//            m_knots_1_count,
//            m_knots_array_1,
//            m_control_points_1_count,
//            weighed_control_points + m_dim * m_control_points_1_count * j);
//        cox_de_boor_algorithm.getDerPoint(1, p_param_1, control_points_2 + m_dim * j);
//        cox_de_boor_algorithm.clear();

        m_numerator_cox_de_boor_algorithm1[j]->getDerPoint(1, p_param_1, control_points_2 + (m_dim + 1) * j);
    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
        m_dim + 1,
        m_order_2,
        m_knots_2_count,
        m_knots_array_2,
        m_control_points_2_count,
        control_points_2);
    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_2, o_point);

    return 0;
}


int NURBSSurface::getDer2Numerator(double p_param_1, double p_param_2, double *o_point) const
{
    double *control_points_1 = &m_control_points_1_array[
        omp_get_thread_num() * (m_dim + 1) * m_control_points_1_count];
//    double *weighed_control_points_transpose = &m_weighed_control_points_transpose_array[omp_get_thread_num() * m_control_points_vector_size];

//    getWeighedControlPointsTranspose(weighed_control_points_transpose);

//    CoxDeBoorAlgorithm cox_de_boor_algorithm;

    for (unsigned i = 0; i < m_control_points_1_count; ++i)
    {
//        cox_de_boor_algorithm.init(
//            m_dim,
//            m_order_2,
//            m_knots_2_count,
//            m_knots_array_2,
//            m_control_points_2_count,
//            weighed_control_points_transpose + m_dim * m_control_points_2_count * i);
//        cox_de_boor_algorithm.getDerPoint(1, p_param_2, control_points_1 + m_dim * i);
//        cox_de_boor_algorithm.clear();

        m_numerator_cox_de_boor_algorithm2[i]->getDerPoint(1, p_param_2, control_points_1 + (m_dim + 1) * i);
    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
        m_dim + 1,
        m_order_1,
        m_knots_1_count,
        m_knots_array_1,
        m_control_points_1_count,
        control_points_1);
    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_1, o_point);

    return 0;
}


int NURBSSurface::getDerNumerator(unsigned p_der1_order,
                                  unsigned p_der2_order,
                                  double p_param_1,
                                  double p_param_2,
                                  double *o_point) const
{
    double *control_points_2 = &m_control_points_2_array[
        omp_get_thread_num() * (m_dim + 1) * m_control_points_2_count];
//    double *weighed_control_points = &m_weighed_control_points_array[omp_get_thread_num() * m_control_points_vector_size];

//    getWeighedControlPoints(weighed_control_points);

//    CoxDeBoorAlgorithm cox_de_boor_algorithm;

    for (unsigned j = 0; j < m_control_points_2_count; ++j)
    {
//        cox_de_boor_algorithm.init(
//            m_dim,
//            m_order_1,
//            m_knots_1_count,
//            m_knots_array_1,
//            m_control_points_1_count,
//            weighed_control_points + m_dim * m_control_points_1_count * j);
        if (p_der1_order == 0)
//            cox_de_boor_algorithm.getPoint(p_param_1, control_points_2 + m_dim * j);
            m_numerator_cox_de_boor_algorithm1[j]->getPoint(p_param_1, control_points_2 + (m_dim + 1) * j);
        else
//            cox_de_boor_algorithm.getDerPoint(p_der1_order, p_param_1, control_points_2 + m_dim * j);
            m_numerator_cox_de_boor_algorithm1[j]->getDerPoint(p_der1_order, p_param_1, control_points_2 + (m_dim + 1) * j);
//        cox_de_boor_algorithm.clear();
    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
        m_dim + 1,
        m_order_2,
        m_knots_2_count,
        m_knots_array_2,
        m_control_points_2_count,
        control_points_2);
    if (p_der2_order == 0)
        m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_2, o_point);
    else
        m_cox_de_boor_algorithm_array[omp_get_thread_num()].getDerPoint(p_der2_order, p_param_2, o_point);

    return 0;
}


//int NURBSSurface::getDer1Denominator(double p_param_1, double p_param_2, double &o_value) const
//{
//    double *control_points_2 = &m_control_points_2_array[omp_get_thread_num() * m_dim * m_control_points_2_count];

////    CoxDeBoorAlgorithm cox_de_boor_algorithm;

//    for (unsigned j = 0; j < m_control_points_2_count; ++j)
//    {
////        cox_de_boor_algorithm.init(
////            1,
////            m_order_1,
////            m_knots_1_count,
////            m_knots_array_1,
////            m_control_points_1_count,
////            m_weights_array + m_control_points_1_count * j);
////        cox_de_boor_algorithm.getDerPoint(1, p_param_1, control_points_2 + j);
////        cox_de_boor_algorithm.clear();

//        m_denominator_cox_de_boor_algorithm1[j]->getDerPoint(1, p_param_1, control_points_2 + j);

//    }

//    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
//        1,
//        m_order_2,
//        m_knots_2_count,
//        m_knots_array_2,
//        m_control_points_2_count,
//        control_points_2);
//    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_2, &o_value);

//    return 0;
//}


//int NURBSSurface::getDer2Denominator(double p_param_1, double p_param_2, double &o_value) const
//{
//    double *control_points_1 = &m_control_points_1_array[omp_get_thread_num() * m_dim * m_control_points_1_count];
////    double *weights_array_transpose = &m_weighed_control_points_transpose_array[omp_get_thread_num() * m_control_points_vector_size];

////    transpose(m_weights_array,
////              1,
////              m_control_points_1_count,
////              m_control_points_2_count,
////              weights_array_transpose);

////    CoxDeBoorAlgorithm cox_de_boor_algorithm;

//    for (unsigned i = 0; i < m_control_points_1_count; ++i)
//    {
////        cox_de_boor_algorithm.init(
////            1,
////            m_order_2,
////            m_knots_2_count,
////            m_knots_array_2,
////            m_control_points_2_count,
////            weights_array_transpose + m_control_points_2_count * i);
////        cox_de_boor_algorithm.getDerPoint(1, p_param_2, control_points_1 + i);
////        cox_de_boor_algorithm.clear();

//        m_denominator_cox_de_boor_algorithm2[i]->getDerPoint(1, p_param_2, control_points_1 + i);
//    }

//    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
//        1,
//        m_order_1,
//        m_knots_1_count,
//        m_knots_array_1,
//        m_control_points_1_count,
//        control_points_1);
//    m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_1, &o_value);

//    return 0;
//}


int NURBSSurface::getDerDenominator(unsigned p_der1_order,
                                    unsigned p_der2_order,
                                    double p_param_1,
                                    double p_param_2,
                                    double &o_value) const
{
    double *control_points_2 = &m_weights_2_array[
        omp_get_thread_num() * m_control_points_2_count];

//    CoxDeBoorAlgorithm cox_de_boor_algorithm;

    for (unsigned j = 0; j < m_control_points_2_count; ++j)
    {
//        cox_de_boor_algorithm.init(
//            1,
//            m_order_1,
//            m_knots_1_count,
//            m_knots_array_1,
//            m_control_points_1_count,
//            m_weights_array + m_control_points_1_count * j);
        if (p_der1_order == 0)
//            cox_de_boor_algorithm.getPoint(p_param_1, control_points_2 + j);
            m_denominator_cox_de_boor_algorithm1[j]->getPoint(p_param_1, control_points_2 + j);
        else
//            cox_de_boor_algorithm.getDerPoint(p_der1_order, p_param_1, control_points_2 + j);
            m_denominator_cox_de_boor_algorithm1[j]->getDerPoint(p_der1_order, p_param_1, control_points_2 + j);
//        cox_de_boor_algorithm.clear();
    }

    m_cox_de_boor_algorithm_array[omp_get_thread_num()].init(
        1,
        m_order_2,
        m_knots_2_count,
        m_knots_array_2,
        m_control_points_2_count,
        control_points_2);
    if (p_der2_order == 0)
        m_cox_de_boor_algorithm_array[omp_get_thread_num()].getPoint(p_param_2, &o_value);
    else
        m_cox_de_boor_algorithm_array[omp_get_thread_num()].getDerPoint(p_der2_order, p_param_2, &o_value);

    return 0;
}


int NURBSSurface::getWeightedControlPoints(
    double *o_weighted_control_points) const
{
#   define _WEIGHTS_INDEX(i, j) ((i) + m_control_points_1_count * (j))
#   define _CONTROL_POINTS_INDEX(i, j, coord) ((coord) + m_dim * ((i) + m_control_points_1_count * (j)))
#   define _WEIGHTED_CONTROL_POINTS_INDEX(i, j, coord) ((coord) + (m_dim + 1) * ((i) + m_control_points_1_count * (j)))

    for (unsigned i = 0; i < m_control_points_1_count; ++i)
        for (unsigned j = 0; j < m_control_points_2_count; ++j) {
            for (unsigned d = 0; d < m_dim; ++d)
                o_weighted_control_points[_WEIGHTED_CONTROL_POINTS_INDEX(i, j, d)] =
                    m_weights_array[_WEIGHTS_INDEX(i, j)] *
                    m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
            o_weighted_control_points[_WEIGHTED_CONTROL_POINTS_INDEX(i, j, m_dim)] =
                m_weights_array[_WEIGHTS_INDEX(i, j)];
        }

#   undef _WEIGHTS_INDEX
#   undef _CONTROL_POINTS_INDEX
#   undef _WEIGHTED_CONTROL_POINTS_INDEX
    return 0;
}


int NURBSSurface::getWeightedControlPointsTranspose(
    double *o_weighed_control_points_transpose) const
{
#   define _WEIGHTS_INDEX(i, j) ((i) + m_control_points_1_count * (j))
#   define _CONTROL_POINTS_INDEX(i, j, coord) ((coord) + m_dim * ((i) + m_control_points_1_count * (j)))
#   define _DEST(i, j, coord) ((coord) + (m_dim + 1) * ((i) + m_control_points_2_count * (j)))

    for (unsigned i = 0; i < m_control_points_1_count; ++i)
        for (unsigned j = 0; j < m_control_points_2_count; ++j) {
            for (unsigned d = 0; d < m_dim; ++d)
                o_weighed_control_points_transpose[_DEST(j, i, d)] =
                    m_weights_array[_WEIGHTS_INDEX(i, j)] *
                    m_control_points_array[_CONTROL_POINTS_INDEX(i, j, d)];
            o_weighed_control_points_transpose[_DEST(j, i, m_dim)] =
                m_weights_array[_WEIGHTS_INDEX(i, j)];
        }

#   undef _WEIGHTS_INDEX
#   undef _CONTROL_POINTS_INDEX
#   undef _DEST
    return 0;
}


void NURBSSurface::transpose(
    const double *p_points,
    unsigned      p_dim,
    unsigned      p_count_1,
    unsigned      p_count_2,
    double       *o_res) const
{
#   define _SOURCE(i, j, coord) ((coord) + p_dim * ((i) + p_count_1 * (j)))
#   define _DEST(i, j, coord) ((coord) + p_dim * ((i) + p_count_2 * (j)))

    for (unsigned i = 0; i < p_count_1; ++i)
        for (unsigned j = 0; j < p_count_2; ++j)
            for (unsigned d = 0; d < p_dim; ++d)
                o_res[_DEST(j, i, d)] = p_points[_SOURCE(i, j, d)];

#   undef _SOURCE
#   undef _DEST
}
