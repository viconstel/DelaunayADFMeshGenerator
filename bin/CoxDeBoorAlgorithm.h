#ifndef COX_DE_BOOR_ALGORITHM_H
#define COX_DE_BOOR_ALGORITHM_H

#include <vector>
#include <memory>


/// \brief Реализация алгоритма Кокса - де Бура для вычисления B-сплайновой кривой и ее производных
class CoxDeBoorAlgorithm
{
public:

    CoxDeBoorAlgorithm();

    /// \brief Конструктор.
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_spline_order Порядок сплайна (степень базисных функций + 1).
    /// \param [in] p_knots_count Количество узловых точек
    /// \param [in] p_knot_vector Узловые точки.
    /// \param [in] p_control_points_count Количество контрольных точек.
    /// \param [in] p_control_points Контрольные точки.
    CoxDeBoorAlgorithm(unsigned      p_dim,
                       unsigned      p_spline_order,
                       unsigned      p_knots_count,
                       const double *p_knot_vector,
                       unsigned      p_control_points_count,
                       const double *p_control_points);

    /// \brief Метод инициализации
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_spline_order Порядок сплайна (степень базисных функций + 1).
    /// \param [in] p_knots_count Количество узловых точек
    /// \param [in] p_knot_vector Узловые точки.
    /// \param [in] p_control_points_count Количество контрольных точек.
    /// \param [in] p_control_points Контрольные точки.
    void init(unsigned      p_dim,
              unsigned      p_spline_order,
              unsigned      p_knots_count,
              const double *p_knot_vector,
              unsigned      p_control_points_count,
              const double *p_control_points);

//    /// \brief Очистить
//    void clear();

    /// \brief Метод для вычисления положения точки на B-сплайновой кривой по параметру кривой
    /// \param [in] p_param Значение параметра кривой
    /// \param [out] o_point Точка на B-сплайновой кривой
    /// \return 0 или код ошибки
    int getPoint(double p_param, double *o_point) const;
    
    /// \brief Метод для вычисления производной данного порядка при заданном значении параметра B-сплайновой кривой.
    /// \param [in] p_der_order Порядок производной.
    /// \param [in] p_param Значение параметра кривой.
    /// \param [out] o_der_point Значение производной в точке на B-сплайновой кривой
    /// \return 0 или код ошибки
    int getDerPoint(unsigned p_der_order, double p_param, double *o_der_point) const;

private:

    ///< Метод определяет такой индекс i, что выполняется условие m_KnotVector[i] <= p_Param < m_KnotVector[i + 1]
    unsigned getParamInd(double p_param) const;

    void initializeDependentData();

private:

    unsigned      m_dim;                  ///< Размерность пространсва
    unsigned      m_spline_order;         ///< Порядок сплайнов
    unsigned      m_knots_count;          ///< Количество узловых точек
    const double *m_knot_vector;          ///< Вектор узлов
    unsigned      m_control_points_count; ///< Количество контрольных точек
    const double *m_control_points;       ///< Контрольные точки

    unsigned m_max_dim;
    unsigned m_max_der_order;
    int m_threads_count;
    std::vector<double> m_compressed_control_points;
    std::vector<double> m_recursion_control_points;
    std::vector<CoxDeBoorAlgorithm> m_cox_de_boor_algorithm;

    double *m_compressed_control_points_array;
    double *m_recursion_control_points_array;
    CoxDeBoorAlgorithm *m_cox_de_boor_algorithm_array;
};


/// \brief Реализация алгоритма расчета базисных функций N
/// Реализован по книге:
/// Les Piegl, Wayne Tiller. The NURBS Book. 2nd Edition. Paragraph 2.5
class NFunction
{
public:

    NFunction();

    /// \brief Конструктор.
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_spline_order Порядок сплайна (степень базисных функций + 1).
    /// \param [in] p_knots_count Количество узловых точек
    /// \param [in] p_knot_vector Узловые точки.
    NFunction(unsigned      p_dim,
              unsigned      p_spline_order,
              unsigned      p_knots_count,
              const double *p_knot_vector);

    ~NFunction();

    /// \brief Метод инициализации
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_order Порядок сплайна (степень базисных функций + 1).
    /// \param [in] p_knots_count Количество узловых точек
    /// \param [in] p_knot_vector Узловые точки.
    /// \param [in] p_control_points_count Количество контрольных точек.
    /// \param [in] p_control_points Контрольные точки.
    void init(unsigned      p_dim,
              unsigned      p_order,
              unsigned      p_knots_count,
              const double *p_knot_vector);

    /// \brief Очистить
    void clear();

//    /// \brief Функтор для вычисления значения функции N при заданном параметре
//    /// \param [in] p_k Номер функции N_k
//    /// \param [in] p_param Значение параметра
//    /// \return Значение функции N
//    double operator()(unsigned p_k, double p_param) const;

    /// \brief Метод для вычисления отличных от нуля базисных функций при заданном параметре
    /// \param [in] p_param Значение параметра кривой
    /// \param [out] o_N Значения базисных функций, отличных от нуля
    /// \param [out] o_param_ind Индекс первой базисной функции, отличной от нуля
    /// \return 0 или код ошибки
    int getValues(const double p_param, const double* &o_N, unsigned &o_param_ind) const;


    /// \brief Метод для вычисления отличных от нуля базисных функций и их производных при заданном параметре
    /// \param [in] p_param Значение параметра кривой
    /// \param [out] o_N Значения базисных функций, отличных от нуля
    /// \param [out] o_param_ind Индекс первой базисной функции, отличной от нуля
    /// \return 0 или код ошибки
    int getDerValues(const double p_param, const double* &o_ders, unsigned &o_param_ind) const;

private:

    ///< Метод определяет такой индекс i, что выполняется условие
    ///< m_KnotVector[i] <= p_Param < m_KnotVector[i + 1]
    ///< Предполагается, что вектор узлов является открытым
    unsigned getParamInd(double p_param) const;

    /// \brief Метод кэширует внутреннее представление векторов
    void cacheArrays();

    void initializeDependentData();

private:

    unsigned      m_dim;                  ///< Размерность пространсва
    unsigned      m_order;         ///< Порядок сплайнов
    unsigned      m_knots_count;          ///< Количество узловых точек
    const double *m_knot_vector;          ///< Вектор узлов

    std::vector<double> m_N;
    std::vector<double> m_left;
    std::vector<double> m_right;
    std::vector<double> m_ndu;
    std::vector<double> m_a;
    std::vector<double> m_ders;

    double       *m_N_array;                    ///< Функции N^{q-1}
    double       *m_left_array;
    double       *m_right_array;
    double       *m_ndu_array;
    double       *m_a_array;
    double       *m_ders_array;
};


#endif // COX_DE_BOOR_ALGORITHM_H
