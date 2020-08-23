#ifndef NURBS_SURFACE_H
#define NURBS_SURFACE_H

#include "NURBSCurve.h"
#include "AbstractSurface.h"


/// Класс NURBS-поверхности
class NURBSSurface : public AbstractSurface
{
public:

    /// \brief Конструктор по умолчанию
    NURBSSurface();

    /// \brief Конструктор поверхности с идентификатором
    /// \param [in] p_id Идентификатор поверхности.
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_sum_index_1 Верхний индекс суммы, такой что (p_SumIndex1 + 1) / Dim - кол-во контрольных точек вдоль первого направления.
    /// \param [in] p_sum_index_2 Верхний индекс суммы, такой что (p_SumIndex2 + 1) / Dim - кол-во контрольных точек вдоль второго направления.
    /// \param [in] p_degree_1 Степень первого множества базисных функций.
    /// \param [in] p_degree_2 Степень второго множества базисных функций.
    /// \param [in] p_knot_vector_1 Первый вектор узловых точек.
    /// \param [in] p_knot_vector_2 Второй вектор узловых точек.
    /// \param [in] p_weights Веса.
    /// \param [in] p_control_points Контрольные точки.
    /// \param [in] p_begin_1 Начальное значение параметра в первом направлении.
    /// \param [in] p_end_1 Конечное значение параметра в первом направлении.
    /// \param [in] p_begin_2 Начальное значение параметра во втором направлении.
    /// \param [in] p_end_2 Конечное значение параметра во втором направлении.
    /// \param [in] p_is_closed_1 TRUE, если поверхность замкнута по первому параметру, иначе - FALSE.
    /// \param [in] p_is_closed_2 TRUE, если поверхность замкнута по вnорому параметру, иначе - FALSE.
    /// \param [in] p_poles Полюсы поверхности.
    NURBSSurface(
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
        double                     p_begin_1 = 0.,
        double                     p_end_1 = 1.,
        double                     p_begin_2 = 0.,
        double                     p_end_2 = 1.,
        bool                       p_is_closed_1 = false,
        bool                       p_is_closed_2 = false,
        std::map<PoleKeyType, PoleValType> p_poles = {});

    /// \brief Конструктор поверхности без идентификатора
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_sum_index_1 Верхний индекс суммы, такой что (p_SumIndex1 + 1) / Dim - кол-во контрольных точек вдоль первого направления.
    /// \param [in] p_sum_index_2 Верхний индекс суммы, такой что (p_SumIndex2 + 1) / Dim - кол-во контрольных точек вдоль второго направления.
    /// \param [in] p_degree_1 Степень первого множества базисных функций.
    /// \param [in] p_degree_2 Степень второго множества базисных функций.
    /// \param [in] p_knot_vector_1 Первый вектор узловых точек.
    /// \param [in] p_knot_vector_2 Второй вектор узловых точек.
    /// \param [in] p_weights Веса.
    /// \param [in] p_control_points Контрольные точки.
    /// \param [in] p_begin_1 Начальное значение параметра в первом направлении.
    /// \param [in] p_end_1 Конечное значение параметра в первом направлении.
    /// \param [in] p_begin_2 Начальное значение параметра во втором направлении.
    /// \param [in] p_end_2 Конечное значение параметра во втором направлении.
    /// \param [in] p_is_closed_1 TRUE, если поверхность замкнута по первому параметру, иначе - FALSE.
    /// \param [in] p_is_closed_2 TRUE, если поверхность замкнута по вnорому параметру, иначе - FALSE.
    /// \param [in] p_poles Полюсы поверхности.
    NURBSSurface(
        unsigned                   p_dim,
        unsigned                   p_sum_index_1,
        unsigned                   p_sum_index_2,
        unsigned                   p_degree_1,
        unsigned                   p_degree_2,
        const std::vector<double> &p_knot_vector_1,
        const std::vector<double> &p_knot_vector_2,
        const std::vector<double> &p_weights,
        const std::vector<double> &p_control_points,
        double                     p_begin_1 = 0.,
        double                     p_end_1 = 1.,
        double                     p_begin_2 = 0.,
        double                     p_end_2 = 1.,
        bool                       p_is_closed_1 = false,
        bool                       p_is_closed_2 = false,
            std::map<PoleKeyType, PoleValType> p_poles = {});

    /// \brief Конструктор поверхности на основе данных контейнера
    /// \param [in] p_parameters Контейнер откуда будут загружаться данные.
    /// \param [in,out] p_n Текущая позиция в контейнере откуда будут загружаться данные.
   // NURBSSurface(const geo::PrimitiveParameters &p_parameters, unsigned &p_n);

    /// \brief Конструктор копирования
    NURBSSurface(const NURBSSurface &p_nurbs_surface);

    /// \brief Оператор присваивания
    const NURBSSurface& operator=(const NURBSSurface &p_nurbs_surface);

    /// \brief Вычисление положения точки на NURBS
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Точка на поверхности
    void getPoint(double p_param_1, double p_param_2, double *o_point) const override;

    /// \brief Метод для вычисления производной NURBS-поверхности в данной точке по первому параметру
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_DerPoint Значение производной в точке на NURBS-поверхности
    void getDer1Point(double p_param_1, double p_param_2, double *o_der_1_point) const override;

    /// \brief Метод для вычисления производной NURBS-поверхности в данной точке по второму параметру
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_DerPoint Значение производной в точке на NURBS-поверхности
    void getDer2Point(double p_param_1, double p_param_2, double *o_der_2_point) const override;

    /// \brief Метод для вычисления (k,l)-производных от 0-го до d-ого порядка (0 <= k + l <= d).
    /// \param [in] p_der_order Значение d.
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Значения производных по параметру в точке на NURBS-поверхности
    void getPointAndDerivs(unsigned p_der_order,
                           double p_param_1,
                           double p_param_2,
                           double *o_der_point) const override;

    /// \brief Нахождение начального приближения проекции точки
    /// \param [in] p_point Координаты проецируемой точки
    /// \param [out] o_point_parameters Параметры начального приближения проекции
    void findInitProjectPoint(
        const double *p_point,
        double       *o_point_parameters) const override;

    /// \brief Метод возвращает поверхность, отражённую относительно оси 2 (ось 1 инвертируется)
    std::shared_ptr<AbstractSurface> reflect2Surface() const override;

    /// \brief Установить новые значения для построения NURBS-поверхности
    /// \param [in] p_degree_1 Степень первого множества базисных функций.
    /// \param [in] p_degree_2 Степень второго множества базисных функций.
    /// \param [in] p_knot_vector_1 Первый вектор узловых точек.
    /// \param [in] p_knot_vector_2 Второй вектор узловых точек.
    /// \param [in] p_control_points Задающие точки, через которые должен пройти сплайн.
    void setVariableData(
            unsigned                   p_degree_1,
            unsigned                   p_degree_2,
            const std::vector<double> &p_knot_vector_1,
            const std::vector<double> &p_knot_vector_2,
            const std::vector<double> &p_control_points);

    /// \brief Установить новые значения задающих точек, через которые должен пройти сплайн
    /// \param [in] p_control_points Задающие точки, через которые должен пройти сплайн.
    void setControlPoints(
            const std::vector<double> &p_control_points);

    /// \brief Вернуть степень NURBS поверхности по направлению 1
    unsigned degree1();

    /// \brief Вернуть степень NURBS поверхности по направлению 2
    unsigned degree2();

    /// \brief Вернуть веса NURBS поверхности
    const std::vector<double> &weights();

    /// \brief Вернуть контрольные точки NURBS поверхности
    const std::vector<double> &controlPoints();

    /// \brief Вернуть вектор узлов NURBS поверхности по направлению 1
    const std::vector<double> &knotVector1();

    /// \brief Вернуть вектор узлов NURBS поверхности по направлению 2
    const std::vector<double> &knotVector2();

    /// \brief Компиляторонезависимая реализация type_info.name()
    /*
    const std::string& type() override
    {
        return m_type;
    }
    static inline const std::string m_type = "NURBSSurface"; // TODO: Поля класса д.б. закрытыми
    */

    void calculateBoundaryCurves() override;

private:

    /// \brief Вычисление числителя
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Значение числителя в данной точке
    /// \return 0 или код ошибки
    int getNumerator(double p_param_1, double p_param_2, double *o_point) const;

    /// \brief Вычисление знаменателя
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_value Значение знаменателя в данной точке
    /// \return 0 или код ошибки
    int getDenominator(double p_param_1, double p_param_2, double &o_value) const;

    /// \brief Вычисление производной числителя по первому параметру
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Значение числителя в данной точке
    /// \return 0 или код ошибки
    int getDer1Numerator(double p_param_1, double p_param_2, double *o_point) const;

    /// \brief Вычисление производной числителя по второму параметру
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_Param2 Значение второго параметра
    /// \param [out] o_point Значение числителя в данной точке
    /// \return 0 или код ошибки
    int getDer2Numerator(double p_param_1, double p_param_2, double *o_point) const;

    /// \brief Вычисление (k,l)-производной числителя
    /// \param [in] p_der_order Порядок производной по первому параметру k.
    /// \param [in] p_der_order Порядок производной по второму параметру l.
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Значение числителя в данной точке
    /// \return 0 или код ошибки
    int getDerNumerator(unsigned p_der1_order,
                        unsigned p_der2_order,
                        double p_param_1,
                        double p_param_2,
                        double *o_point) const;

//    /// \brief Вычисление производной знаменателя по первому параметру
//    /// \param [in] p_param_1 Значение первого параметра
//    /// \param [in] p_param_2 Значение второго параметра
//    /// \param [out] o_value Значение знаменателя в данной точке
//    /// \return 0 или код ошибки
//    int getDer1Denominator(double p_param_1, double p_param_2, double &o_value) const;

//    /// \brief Вычисление производной знаменателя по второму параметру
//    /// \param [in] p_param_1 Значение первого параметра
//    /// \param [in] p_param_2 Значение второго параметра
//    /// \param [out] o_value Значение знаменателя в данной точке
//    /// \return 0 или код ошибки
//    int getDer2Denominator(double p_param_1, double p_param_2, double &o_value) const;

    /// \brief Вычисление (k,l)-производной знаменателя
    /// \param [in] p_der_order Порядок производной по первому параметру k.
    /// \param [in] p_der_order Порядок производной по второму параметру l.
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_value Значение знаменателя в данной точке
    /// \return 0 или код ошибки
    int getDerDenominator(unsigned p_der1_order,
                          unsigned p_der2_order,
                          double p_param_1,
                          double p_param_2,
                          double &o_value) const;

    /// \brief Вычисление произведения контрольных точек на вес
    /// \param [Out] o_weighed_control_points Контрольные точки, умноженные на вес
    /// \return 0 или код ошибки
    int getWeightedControlPoints(double *o_weighted_control_points) const;

    /// \brief Вычисление произведения контрольных точек на вес с перестановкой.
    /// \param [Out] o_weighed_control_points_transpose Контрольные точки, умноженные на вес
    /// \return 0 или код ошибки
    int getWeightedControlPointsTranspose(double *o_weighted_control_points_transpose) const;

    /// \brief Метод изменяет порядок хранения контрольных точек и весов
    /// \param [in] p_points Массив точек
    /// \param [in] p_dim Размерность пространства
    /// \param [in] p_count_1 Количество точек вдоль первого направления
    /// \param [in] p_count_2 Количество точек вдоль второго направления
    /// \param [out] o_res Моссив точек в новом порядке
    void transpose(const double *p_points, unsigned p_dim, unsigned p_count_1, unsigned p_count_2, double *o_res) const;

    /// \brief Метод кэширует внутреннее представление векторов
    void cacheArrays();

    void initializeDependentData();

private:
    std::vector<double> m_knot_vector_1;
    std::vector<double> m_knot_vector_2;
    std::vector<double> m_weights;
    std::vector<double> m_control_points;
    unsigned            m_sum_index_1;
    unsigned            m_sum_index_2;
    unsigned            m_degree_1;
    unsigned            m_degree_2;

    std::vector<double> m_nurbs_numerator_point;
    std::vector<double> m_nurbs_numerator_der_point;
    std::vector<double> m_v;
    std::vector<double> m_v2;
    std::vector<double> m_control_points_1;
    std::vector<double> m_control_points_2;
    std::vector<double> m_weights_1;
    std::vector<double> m_weights_2;
    std::vector<double> m_weights_transpose;
    std::vector<double> m_weighted_control_points;
    std::vector<double> m_weighted_control_points_transpose;
    std::vector<double> m_Sw;

    std::vector<std::vector<unsigned>>  m_bin;   ///< Binomial coefficients C^k_n

    std::vector<CoxDeBoorAlgorithm> m_cox_de_boor_algorithm;
    std::vector<std::unique_ptr<CoxDeBoorAlgorithm>> m_numerator_cox_de_boor_algorithm1;
    std::vector<std::unique_ptr<CoxDeBoorAlgorithm>> m_numerator_cox_de_boor_algorithm2;
    std::vector<std::unique_ptr<CoxDeBoorAlgorithm>> m_denominator_cox_de_boor_algorithm1;
    std::vector<std::unique_ptr<CoxDeBoorAlgorithm>> m_denominator_cox_de_boor_algorithm2;

    // Кэшированные данные

    unsigned m_order_1;
    unsigned m_order_2;
    unsigned m_knots_1_count;
    unsigned m_knots_2_count;
    unsigned m_control_points_vector_size;
    unsigned m_control_points_1_count;
    unsigned m_control_points_2_count;

    const double *m_knots_array_1;
    const double *m_knots_array_2;
    const double *m_weights_array;
    const double *m_control_points_array;

    double *m_nurbs_numerator_point_array;
    double *m_nurbs_numerator_der_point_array;
    double *m_v_array;
    double *m_v2_array;
    double *m_control_points_1_array;
    double *m_control_points_2_array;
    double *m_weights_1_array;
    double *m_weights_2_array;
    double *m_weights_array_transpose;
    double *m_weighted_control_points_array;
    double *m_weighted_control_points_transpose_array;
    CoxDeBoorAlgorithm *m_cox_de_boor_algorithm_array;
    double *m_Sw_array;
};


#endif // NURBS_SURFACE_H

