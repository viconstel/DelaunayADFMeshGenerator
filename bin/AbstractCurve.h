#ifndef ABSTRACT_CURVE_H
#define ABSTRACT_CURVE_H

//#include <GeometryData>
#include <string>
#include <vector>

/// \brief Класс абстрактной параметрической кривой
class AbstractCurve
{
public:
    /// \brief Конструктор абстрактной кривой без ID
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_begin Начальное значение параметра кривой.
    /// \param [in] p_end Конечное значение параметра кривой.
    /// \param [in] p_is_closed true, если кривая замкнутая, иначе false.
    AbstractCurve(
        unsigned p_dim = 0,
        double   p_begin = 0.,
        double   p_end = 1.,
        bool     p_is_closed = false);

    /// \brief Конструктор абстрактной кривой с ID
    /// \param [in] p_id Идентификатор кривой
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_begin Начальное значение параметра кривой.
    /// \param [in] p_end Конечное значение параметра кривой.
    /// \param [in] p_is_closed true, если кривая замкнутая, иначе false.
    explicit AbstractCurve(
        unsigned p_id,
        unsigned p_dim,
        double   p_begin = 0.,
        double   p_end = 1.,
        bool     p_is_closed = false);

    /// \brief Конструктор кривой на основе данных контейнера
    /// \param [in] p_parameters Контейнер откуда будут загружаться данные.
    /// \param [in,out] p_n Текущая позиция в контейнере откуда будут загружаться данные.
   // explicit AbstractCurve(const geo::PrimitiveParameters &p_parameters, unsigned &p_n);

    /// \brief Конструктор копирования
    AbstractCurve(const AbstractCurve &p_curve) = default;

    /// \brief Деструктор
    virtual ~AbstractCurve();

    /// \brief Оператор присваивания
    AbstractCurve& operator=(const AbstractCurve &p_curve) = default;

    /// \brief Метод для вычисления положения точки на кривой по параметру кривой
    /// \param [in] p_param Значение параметра кривой
    /// \param [out] o_point Искомая точка
    virtual void getPoint(double p_param, double *o_point) const = 0;

    /// \brief Метод для вычисления производной данного порядка при заданном значении параметра кривой.
    /// \param [in] p_param Значение параметра кривой.
    /// \param [out] o_point Значение производной по параметру в точке на NURBS-поверхности
    virtual void getDerPoint(double p_param, double *o_der_point) const = 0;

    /// \brief Метод для вычисления производных от 0-го до d-ого порядка включительно.
    /// \param [in] p_der_order Порядок производной d.
    /// \param [in] p_param Значение параметра кривой.
    /// \param [out] o_point Значения производных по параметру в точке на кривой
    virtual void getPointAndDerivs(unsigned p_der_order,
                                   double p_param,
                                   double *o_der_point) const = 0;
    
    /// \brief Возвращение значения параметра кривой, отвечающего ее началу
    double getBegin() const;

    /// \brief Возвращение значения параметра кривой, отвечающего ее концу
    double getEnd() const;

    /// \brief TRUE, если кривая замкнутая и FALSE - в противном случае
    bool isClosed() const;

    /// \brief Возвращение значения размерности пространства
    unsigned getDim() const;

    /// \brief Возвращение идентификатора кривой
    unsigned getId() const;

    /// \brief Компиляторонезависимая реализация type_info.name()
    /*
    virtual const std::string& type()
    {
        return m_type;
    }
    static inline const std::string m_type = "AbstractCurve"; // TODO: данные-члены должны быть закрытыми
    */
protected:
    unsigned m_id{ 0 };            ///< Идентификатор кривой
    unsigned m_dim;                ///< Размерность пространсва
    double   m_begin{ 0. };        ///< Начальное значение параметра кривой
    double   m_end{ 1. };          ///< Конечное значение параметра кривой
    bool     m_is_closed{ false }; ///< TRUE, если кривая замкнутая и FALSE - в противном случае

    unsigned *m_counts_array;

    unsigned m_max_der_order;
private:
    void initializeDependentData();
    std::vector<unsigned> m_counts;
};


#endif // CURVE_H
