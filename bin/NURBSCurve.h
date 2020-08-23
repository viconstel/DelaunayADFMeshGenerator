#ifndef NURBS_CURVE_H
#define NURBS_CURVE_H

#include "AbstractCurve.h"
#include <algorithm>
#include <cmath>

#include "CoxDeBoorAlgorithm.h"
//#include <BasicMathFunctions>
//#include <GeometryData>


/// Класс NURBS-кривых
class NURBSCurve : public AbstractCurve
{
public:
    /// \brief Коды равенства двух NURBS-кривых
    enum EqualityCodes : uint8_t
    {
        unequal,       ///< Не равны
        forward_equal, ///< Равны в точности
        backward_equal ///< Равны при инвертации порядка контрольных точек
    };

public:
    /// \brief Конструктор по умолчанию
    NURBSCurve();

    /// \brief Конструктор кривой c ID
    /// \param [in] p_id Идентификатор кривой
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_degree Степень базисных функций.
    /// \param [in] p_knot_vector Вектор узлов.
    /// \param [in] p_weights Веса.
    /// \param [in] p_control_points Контрольные точки.
    /// \param [in] p_normal Вектор единичной нормали (для плоской кривой).
    /// \param [in] p_begin Начальное значение параметра кривой.
    /// \param [in] p_end Конечное значение параметра кривой.
    /// \param [in] p_is_closed true, если кривая замкнутая, иначе false.
    NURBSCurve(
        unsigned                   p_id,
        unsigned                   p_dim,
        unsigned                   p_degree,
        const std::vector<double> &p_knot_vector,
        const std::vector<double> &p_weights,
        const std::vector<double> &p_control_points,
        const std::vector<double> &p_normal,
        double                     p_begin = 0.,
        double                     p_end = 1.,
        bool                       p_is_closed = false);

    /// \brief Конструктор кривой без ID
    /// \param [in] p_id Идентификатор кривой
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_degree Степень базисных функций.
    /// \param [in] p_knot_vector Вектор узлов.
    /// \param [in] p_weights Веса.
    /// \param [in] p_control_points Контрольные точки.
    /// \param [in] p_normal Вектор единичной нормали (для плоской кривой).
    /// \param [in] p_begin Начальное значение параметра кривой.
    /// \param [in] p_end Конечное значение параметра кривой.
    /// \param [in] p_is_closed true, если кривая замкнутая, иначе false.
    NURBSCurve(
        unsigned                   p_dim,
        unsigned                   p_degree,
        const std::vector<double> &p_knot_vector,
        const std::vector<double> &p_weights,
        const std::vector<double> &p_control_points,
        const std::vector<double> &p_normal,
        double                     p_begin = 0.,
        double                     p_end = 1.,
        bool                       p_is_closed = false);

    /// \brief Конструктор кривой на основе данных контейнера
    /// \param [in] p_parameters Контейнер откуда будут загружаться данные.
    /// \param [in,out] p_n Текущая позиция в контейнере откуда будут загружаться данные.
   
    //NURBSCurve(const geo::PrimitiveParameters &p_parameters, unsigned &p_n);

    /// \brief Конструктор копирования
    NURBSCurve(const NURBSCurve &p_nurbs_curve);

    /// \brief Оператор присваивания
    const NURBSCurve& operator=(const NURBSCurve &p_nurbs_curve);

    /// \brief Метод для вычисления положения точки на NURBS-кривой по параметру кривой
    /// \param [in] p_param Значение параметра кривой
    /// \param [out] o_point Искомая точка
    void getPoint(double p_param, double *o_point) const override;

    /// \brief Метод для вычисления производной данного порядка при заданном значении параметра кривой.
    /// \param [in] p_param Значение параметра кривой.
    /// \param [out] o_point Значение производной по параметру в точке на NURBS-поверхности
    void getDerPoint(double p_param, double *o_der_point) const override;

    /// \brief Метод для вычисления производных от 0-го до d-ого порядка включительно.
    /// \param [in] p_der_order Порядок производной d.
    /// \param [in] p_param Значение параметра кривой.
    /// \param [out] o_point Значения производных по параметру в точке на NURBS-кривой
    void getPointAndDerivs(unsigned p_der_order, double p_param, double *o_der_point) const override;

    /// \brief Вернуть степень NURBS кривой
    unsigned degree();

    /// \brief Вернуть веса NURBS кривой
    const std::vector<double> &weights();

    /// \brief Вернуть контрольные точки NURBS кривой
    const std::vector<double> &controlPoints();

    /// \brief Вернуть вектор узлов NURBS кривой
    const std::vector<double> &knotVector();

    /// \brief Компиляторонезависимая реализация type_info.name()
    
    /*
    const std::string& type() override
    {
        return m_type;
    }
    static inline const std::string m_type = "NURBSCurve"; // TODO: данные-члены должны быть закрытыми
    */


    /// \brief Метод для расчета открытого вектора узлов.
    /// \param [in] p_sum_index Верхний индекс суммы.
    /// \param [in] p_degree Степень базисных функций.
    static std::vector<double> getOpenKnotVector(unsigned p_sum_index,
                                                 unsigned p_degree);
private:

    /// \brief Вычисление произведения контрольных точек на вес
    /// \param [Out] o_WeighedControlPoints Контрольные точки, умноженные на вес
    /// \return 0 или код ошибки
    int getWeighedControlPoints(double *o_WeighedControlPoints) const;

    /// \brief Метод кэширует внутреннее представление векторов
    void cacheArrays();

    void initializeDependentData();

private:
    unsigned            m_degree;         ///< Степень базисных функций
    std::vector<double> m_weights;        ///< Весовые коэффициенты
    std::vector<double> m_control_points; ///< Контрольные точки
    std::vector<double> m_knot_vector;    ///< Вектор узлов
    std::vector<double> m_normal;         ///< Единичная нормаль (для плоской кривой)

    std::vector<std::vector<unsigned>>  m_bin;   ///< Binomial coefficients C^k_n

    std::vector<double> m_weighed_control_points;
    std::vector<double> m_nurbs_numerator_point;
    std::vector<double> m_nurbs_numerator_der_point;

    std::unique_ptr<CoxDeBoorAlgorithm> m_numerator_cox_de_boor_algorithm;
    std::unique_ptr<CoxDeBoorAlgorithm> m_denominator_cox_de_boor_algorithm;

    // Кэшированные данные

    unsigned m_order;
    unsigned m_knots_count;
    unsigned m_control_points_vector_size;
    unsigned m_control_points_count;

    const double *m_knots_array;
    const double *m_weights_array;
    const double *m_control_points_array;

    double *m_weighed_control_points_array;
    double *m_nurbs_numerator_point_array;
    double *m_nurbs_numerator_der_point_array;

private:
    friend uint8_t operator==(
        const NURBSCurve &p_curve_1,
        const NURBSCurve &p_curve_2);
};


///////////////////////////////////////////////////////////////////////////////
// Операторы сравнения кривых
///////////////////////////////////////////////////////////////////////////////


/// \brief Оператор сравнения двух NURBS-кривых
/// \param p_curve_1 Первая кривая
/// \param p_curve_2 Вторая кривая
/// \return Код результата сравнения: 0 - неравны, равны - в противном случае
uint8_t operator==(const NURBSCurve &p_curve_1, const NURBSCurve &p_curve_2);


#endif // NURBS_CURVE_H
