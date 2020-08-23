#include "NURBSCurve.h"
#include <omp.h>


///////////////////////////////////////////////////////////////////////////////
// NURBSCurve
///////////////////////////////////////////////////////////////////////////////


NURBSCurve::NURBSCurve()
    : m_degree(0)
    , m_order(0)
    , m_knots_count(0)
    , m_control_points_vector_size(0)
    , m_control_points_count(0)
    , m_knots_array(nullptr)
    , m_weights_array(nullptr)
    , m_control_points_array(nullptr)
{}


NURBSCurve::NURBSCurve(
    unsigned                   p_id,
    unsigned                   p_dim,
    unsigned                   p_degree,
    const std::vector<double> &p_knot_vector,
    const std::vector<double> &p_weights,
    const std::vector<double> &p_control_points,
    const std::vector<double> &p_normal,
    double                     p_begin,
    double                     p_end,
    bool                       p_is_closed)
    : AbstractCurve(
          p_id,
          p_dim,
          p_begin,
          p_end,
          p_is_closed)
    , m_degree(p_degree)
    , m_knot_vector(p_knot_vector)
    , m_weights(p_weights)
    , m_control_points(p_control_points)
    , m_normal(p_normal)
{
    initializeDependentData();
}


NURBSCurve::NURBSCurve(
    unsigned                   p_dim,
    unsigned                   p_degree,
    const std::vector<double> &p_knot_vector,
    const std::vector<double> &p_weights,
    const std::vector<double> &p_control_points,
    const std::vector<double> &p_normal,
    double                     p_begin,
    double                     p_end,
    bool                       p_is_closed)
    : NURBSCurve(
          0,
          p_dim,
          p_degree,
          p_knot_vector,
          p_weights,
          p_control_points,
          p_normal,
          p_begin,
          p_end,
          p_is_closed)
{ }

/*NURBSCurve::NURBSCurve(const geo::PrimitiveParameters &p_parameters, unsigned &p_n)
    : AbstractCurve(p_parameters, p_n)
{
    m_degree         = std::any_cast<unsigned           >(p_parameters[p_n++]);
    m_weights        = std::any_cast<std::vector<double>>(p_parameters[p_n++]);
    m_control_points = std::any_cast<std::vector<double>>(p_parameters[p_n++]);
    m_knot_vector    = std::any_cast<std::vector<double>>(p_parameters[p_n++]);
    m_normal         = std::any_cast<std::vector<double>>(p_parameters[p_n++]);

    initializeDependentData();
}
*/

NURBSCurve::NURBSCurve(const NURBSCurve &p_nurbs_curve)
    : AbstractCurve(p_nurbs_curve)
    , m_degree(p_nurbs_curve.m_degree)
    , m_knot_vector(p_nurbs_curve.m_knot_vector)
    , m_weights(p_nurbs_curve.m_weights)
    , m_control_points(p_nurbs_curve.m_control_points)
    , m_normal(p_nurbs_curve.m_normal)
    , m_weighed_control_points(p_nurbs_curve.m_weighed_control_points)
    , m_nurbs_numerator_point(p_nurbs_curve.m_nurbs_numerator_point)
    , m_nurbs_numerator_der_point(p_nurbs_curve.m_nurbs_numerator_der_point)
{
    initializeDependentData();
}


const NURBSCurve& NURBSCurve::operator=(const NURBSCurve &p_nurbs_curve)
{
    if (this == &p_nurbs_curve)
        return *this;

    AbstractCurve::operator=(p_nurbs_curve);

    m_degree = p_nurbs_curve.m_degree;
    m_knot_vector = p_nurbs_curve.m_knot_vector;
    m_weights = p_nurbs_curve.m_weights;
    m_control_points = p_nurbs_curve.m_control_points;
    m_normal = p_nurbs_curve.m_normal;
    m_weighed_control_points = p_nurbs_curve.m_weighed_control_points;
    m_nurbs_numerator_point = p_nurbs_curve.m_nurbs_numerator_point;
    m_nurbs_numerator_der_point = p_nurbs_curve.m_nurbs_numerator_der_point;

    initializeDependentData();

    return *this;
}


void NURBSCurve::initializeDependentData()
{
    // Кэширование данных

    m_order = m_degree + 1;
    m_knots_count = m_knot_vector.size();
    m_control_points_vector_size = m_control_points.size();
    m_control_points_count = static_cast<unsigned>(m_control_points_vector_size / m_dim);

    m_weighed_control_points.resize(m_control_points_vector_size);

    int threads_count = omp_get_max_threads();
    m_nurbs_numerator_point.resize(threads_count * m_dim);
    m_nurbs_numerator_der_point.resize(threads_count * m_dim);

    cacheArrays();
    //basicMathFunctions::bci(m_bin, m_degree);
    // Вычисление взвешенных контрольных точек
    getWeighedControlPoints(m_weighed_control_points_array);

    // Создание объектов класса алгоритма Кокса - де Бура
    m_numerator_cox_de_boor_algorithm = std::make_unique<CoxDeBoorAlgorithm> (
            m_dim,
            m_order,
            m_knots_count,
            m_knots_array,
            m_control_points_count,
            m_weighed_control_points_array);

    // Создание объектов класса алгоритма Кокса - де Бура
    m_denominator_cox_de_boor_algorithm = std::make_unique<CoxDeBoorAlgorithm> (
        1,
        m_order,
        m_knots_count,
        m_knots_array,
        m_control_points_count,
        m_weights_array);
}


void NURBSCurve::cacheArrays()
{
    m_knots_array = m_knot_vector.data();
    m_weights_array = m_weights.data();
    m_control_points_array = m_control_points.data();

    m_weighed_control_points_array = m_weighed_control_points.data();
    m_nurbs_numerator_point_array = m_nurbs_numerator_point.data();
    m_nurbs_numerator_der_point_array = m_nurbs_numerator_der_point.data();
}


void NURBSCurve::getPoint(double p_param, double *o_point) const
{
//    double *weighed_control_points = &m_weighed_control_points_array[omp_get_thread_num() * m_control_points_vector_size];
    double  denominator_value;

//#   pragma omp parallel sections shared(weighed_control_points, denominator_value)
    {
//#       pragma omp section // Секция расчета числителя
        {
//            // Вычисление взвешенных контрольных точек
//            getWeighedControlPoints(weighed_control_points);

//            // Создание объектов класса алгоритма Кокса - де Бура
//            const CoxDeBoorAlgorithm numerator_cox_de_boor_algorithm(
//                m_dim,
//                m_order,
//                m_knots_count,
//                m_knots_array,
//                m_control_points_count,
//                weighed_control_points);

            // Вычисление значения числителя в заданной точке
            m_numerator_cox_de_boor_algorithm->getPoint(p_param, o_point);
        }
//#       pragma omp section // Секция расчета знаменателя
        {
//            // Создание объектов класса алгоритма Кокса - де Бура
//            const CoxDeBoorAlgorithm denominator_cox_de_boor_algorithm(
//                1,
//                m_order,
//                m_knots_count,
//                m_knots_array,
//                m_control_points_count,
//                m_weights_array);

            // Вычисление значения знаменателя в заданной точке
            m_denominator_cox_de_boor_algorithm->getPoint(p_param, &denominator_value);
        }
    }

    // Вычисление точки
    for (unsigned d = 0; d < m_dim; ++d)
        o_point[d] /= denominator_value;
}


void NURBSCurve::getDerPoint(double p_Param, double *o_der_point) const
{
    double *nurbs_numerator_point = &m_nurbs_numerator_point_array[omp_get_thread_num() * m_dim];
    double *nurbs_numerator_der_point = &m_nurbs_numerator_der_point_array[omp_get_thread_num() * m_dim];
    double  nurbs_denominator_value;
    double  nurbs_denominator_der_value;

//    // Вычисление взвешенных контрольных точек
//    getWeighedControlPoints(weighed_control_points);

//    // Создание объектов класса алгоритма Кокса - де Бура
//    const CoxDeBoorAlgorithm nurbs_numerator_cox_de_boor_algorithm(
//        m_dim,
//        m_order,
//        m_knots_count,
//        m_knots_array,
//        m_control_points_count,
//        weighed_control_points);
//    const CoxDeBoorAlgorithm nurbs_denominator_cox_de_boor_algorithm(
//        1,
//        m_order,
//        m_knots_count,
//        m_knots_array,
//        m_control_points_count,
//        m_weights_array);

//#   pragma omp parallel sections \
//           shared(nurbs_numerator_point, nurbs_numerator_der_point, nurbs_denominator_value, nurbs_denominator_der_value)
    // nurbs_numerator_cox_de_boor_algorithm, nurbs_denominator_cox_de_boor_algorithm)
    {
        // Секция расчета числителя
//#       pragma omp section
        m_numerator_cox_de_boor_algorithm->getPoint(p_Param, nurbs_numerator_point);

        // Секция расчета первой производной числителя
//#       pragma omp section
        m_numerator_cox_de_boor_algorithm->getDerPoint(1, p_Param, nurbs_numerator_der_point);

        // Секция расчета знаменателя
//#       pragma omp section
        m_denominator_cox_de_boor_algorithm->getPoint(p_Param, &nurbs_denominator_value);

        // Секция расчета первой производной знаменателя
//#       pragma omp section
        m_denominator_cox_de_boor_algorithm->getDerPoint(1, p_Param, &nurbs_denominator_der_value);
    }

    // Вычисление первой производной NURBS-кривой в данной точке
    for (unsigned d = 0; d < m_dim; ++d)
        o_der_point[d] = (nurbs_numerator_der_point[d] * nurbs_denominator_value - nurbs_numerator_point[d] * nurbs_denominator_der_value) /
                         (nurbs_denominator_value * nurbs_denominator_value);
}


void NURBSCurve::getPointAndDerivs(unsigned p_der_order, double p_param, double *o_der_point) const
{
    unsigned j;

    double  wders0;
    double  wdersi;

#   define _(point, coord) (m_dim * (point) + (coord))
    for (unsigned k = m_degree + 1; k <= p_der_order; ++k)
        for (j = 0; j < m_dim; ++j)
            o_der_point[_(k,j)] = 0;

    p_der_order = std::min(p_der_order, m_degree);

//    // Вычисление взвешенных контрольных точек
//    getWeighedControlPoints(weighed_control_points);

//    // Создание объектов класса алгоритма Кокса - де Бура
//    const CoxDeBoorAlgorithm Aders(
//        m_dim,
//        m_order,
//        m_knots_count,
//        m_knots_array,
//        m_control_points_count,
//        weighed_control_points);
//    const CoxDeBoorAlgorithm wders(
//        1,
//        m_order,
//        m_knots_count,
//        m_knots_array,
//        m_control_points_count,
//        m_weights_array);

    /// Расчет проводится по формулам из книги:
    /// Les Piegl, Wayne Tiller. The NURBS Book. 2nd Edition. Paragraph 4.3
    /// Русский перевод:
    /// Ю.Д. Шевелев. Математические основы задач проектирования:
    /// Учебное пособие по курсу "Математическое и программное обеспечение САПР"
    /// М., 2005. Параграф 2.7

    // ALGORITHM A4.2
    m_denominator_cox_de_boor_algorithm->getPoint(p_param, &wders0);
    for (unsigned k = 0; k <= p_der_order; ++k) {
        if (k == 0)
            m_numerator_cox_de_boor_algorithm->getPoint(p_param, &o_der_point[_(k,0)]);
        else {
            m_numerator_cox_de_boor_algorithm->getDerPoint(k, p_param, &o_der_point[_(k,0)]);
        }
        for (unsigned i=1; i<=k; i++) {
            m_denominator_cox_de_boor_algorithm->getDerPoint(i, p_param, &wdersi);
            for (j = 0; j < m_dim; ++j)
                o_der_point[_(k,j)] = o_der_point[_(k,j)] -
                        m_bin[k][i] * wdersi * o_der_point[_(k-1,j)];
        }
        for (j = 0; j < m_dim; ++j)
            o_der_point[_(k,j)] = o_der_point[_(k,j)] / wders0;
    }
#   undef _
}


unsigned NURBSCurve::degree()
{
    return m_degree;
}


const std::vector<double> &NURBSCurve::weights()
{
    return m_weights;
}


const std::vector<double> &NURBSCurve::controlPoints()
{
    return m_control_points;
}


const std::vector<double> &NURBSCurve::knotVector()
{
    return m_knot_vector;
}


std::vector<double> NURBSCurve::getOpenKnotVector(unsigned p_sum_index, unsigned p_degree)
{
    /// Расчет проводится по формуле (8.63) из книги:
    /// Д. Херн, М.П. Бейкер.
    /// Компьютерная графика и стандарт OpenGL, 3-е издание. : Пер. с англ.
    /// М. : Издательский дом "Вильямс", 2005. 1168 с.

    std::vector<double> knot_vector(p_sum_index + p_degree + 2);

    for (unsigned j = 0; j < p_sum_index + p_degree + 2; ++j) {
        if (j<p_degree+1) {
            knot_vector[j] = 0;
        } else if (j>p_sum_index) {
            knot_vector[j] = p_sum_index-p_degree+1;
        } else {
            knot_vector[j] = j-p_degree;
        }
    }

    return knot_vector;
}


int NURBSCurve::getWeighedControlPoints(double *o_WeighedControlPoints) const
{
#   define _(point, coord) (m_dim * (point) + (coord))
    for (unsigned i = 0; i < m_control_points_count; ++i)
        for (unsigned d = 0; d < m_dim; ++d)
            o_WeighedControlPoints[_(i, d)] = m_weights_array[i] * m_control_points_array[_(i, d)];
#   undef _
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
// Операторы сравнения кривых
///////////////////////////////////////////////////////////////////////////////


/// Оператор сравнения двух NURBS-кривых
uint8_t operator==(const NURBSCurve &p_curve_1, const NURBSCurve &p_curve_2)
{
#   define _(point, coord) (point) * p_curve_1.m_dim + (coord)

    if (&p_curve_1 == &p_curve_2)
        return static_cast<uint8_t>(NURBSCurve::EqualityCodes::forward_equal);

    if (p_curve_1.m_control_points_count != p_curve_2.m_control_points_count)
        return static_cast<uint8_t>(NURBSCurve::EqualityCodes::unequal);

    const double precision = 1.E-6;
    double delta_square{ 0. };
    bool is_forward{ true };

    std::vector<double> delta_vec(p_curve_1.m_dim);

    // Проверка совпадения контрольных точек
    for (unsigned i = 0; i < p_curve_1.m_control_points_count; ++i)
    {
        for (unsigned d = 0; d < p_curve_1.m_dim; ++d)
            delta_vec[d] = p_curve_2.m_control_points_array[_(i, d)] -
                           p_curve_1.m_control_points_array[_(i, d)];

        delta_square = 0.;
        for (unsigned d = 0; d < p_curve_1.m_dim; ++d)
            delta_square += delta_vec[d] * delta_vec[d];

        if (delta_square > precision * precision)
        {
            is_forward = false;
            break;
        }
    }

    if (!is_forward)
    {
        for (unsigned i = 0; i < p_curve_1.m_control_points_count; ++i)
        {
            for (unsigned d = 0; d < p_curve_1.m_dim; ++d)
                delta_vec[d] =
                    p_curve_2.m_control_points_array[_(p_curve_2.m_control_points_count - i - 1, d)] -
                    p_curve_1.m_control_points_array[_(i, d)];

            delta_square = 0.;
            for (unsigned d = 0; d < p_curve_1.m_dim; ++d)
                delta_square += delta_vec[d] * delta_vec[d];

            if (delta_square > precision * precision)
                return static_cast<uint8_t>(NURBSCurve::EqualityCodes::unequal);
        }
    }

    // Проверка совпадения весов.
    // Если контрольные точки совпадают в прямом направлении,
    // то и веса должны совпадать в том же порядке

    if (is_forward)
    {
        for (unsigned i = 0; i < p_curve_1.m_control_points_count; ++i)
        {
            if (std::fabs(p_curve_2.m_weights_array[i] -
                          p_curve_1.m_weights_array[i]) > precision)
                return static_cast<uint8_t>(NURBSCurve::EqualityCodes::unequal);
        }

        return static_cast<uint8_t>(NURBSCurve::EqualityCodes::forward_equal);
    }
    else
    {
        for (unsigned i = 0; i < p_curve_1.m_control_points_count; ++i)
        {
            if (std::fabs(p_curve_2.m_weights_array[p_curve_2.m_control_points_count - i - 1] -
                          p_curve_1.m_weights_array[i]) > precision)
                return static_cast<uint8_t>(NURBSCurve::EqualityCodes::unequal);
        }

        return static_cast<uint8_t>(NURBSCurve::EqualityCodes::backward_equal);
    }

#   undef _
}
