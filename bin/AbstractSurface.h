#ifndef ABSTRACT_SURFACE_H
#define ABSTRACT_SURFACE_H

#include "CoxDeBoorAlgorithm.h"
#include <vector>
#include <algorithm>
#include "AbstractCurve.h"
//#include <LineCurve>
//#include <GeometryData>
#include <unordered_map>
#include <map>

/// \brief Класс абстрактной параметричекой поверхности
class AbstractSurface
{
public:
    using PoleKeyType = std::pair<unsigned, double>;
    using PoleValType = std::pair<double, double>;
public:
    /// \brief Конструктор поверхности без идентификатора
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_begin_1 Начальное значение параметра в первом направлении.
    /// \param [in] p_end_1 Конечное значение параметра в первом направлении.
    /// \param [in] p_begin_2 Начальное значение параметра во втором направлении.
    /// \param [in] p_end_2 Конечное значение параметра во втором направлении.
    /// \param [in] p_is_closed_1 TRUE, если поверхность замкнута по первому параметру, иначе - FALSE.
    /// \param [in] p_is_closed_2 TRUE, если поверхность замкнута по второму параметру, иначе - FALSE.
    AbstractSurface(
        unsigned p_dim         = 0,
        double   p_begin_1     = 0.,
        double   p_end_1       = 1.,
        double   p_begin_2     = 0.,
        double   p_end_2       = 1.,
        bool     p_is_closed_1 = false,
        bool     p_is_closed_2 = false,
            std::map<PoleKeyType, PoleValType> p_poles = {});

    /// \brief Конструктор поверхности с идентификатором
    /// \param [in] p_id Идентификатор поверхности
    /// \param [in] p_dim Размерность пространства.
    /// \param [in] p_begin_1 Начальное значение параметра в первом направлении.
    /// \param [in] p_end_1 Конечное значение параметра в первом направлении.
    /// \param [in] p_begin_2 Начальное значение параметра во втором направлении.
    /// \param [in] p_end_2 Конечное значение параметра во втором направлении.
    /// \param [in] p_is_closed_1 TRUE, если поверхность замкнута по первому параметру, иначе - FALSE.
    /// \param [in] p_is_closed_2 TRUE, если поверхность замкнута по второму параметру, иначе - FALSE.
    explicit AbstractSurface(
        unsigned p_id,
        unsigned p_dim         = 0,
        double   p_begin_1     = 0.,
        double   p_end_1       = 1.,
        double   p_begin_2     = 0.,
        double   p_end_2       = 1.,
        bool     p_is_closed_1 = false,
        bool     p_is_closed_2 = false,
        std::map<PoleKeyType, PoleValType> p_poles = {});

    /// \brief Конструктор поверхности на основе данных контейнера
    /// \param [in] p_parameters Контейнер откуда будут загружаться данные.
    /// \param [in,out] p_n Текущая позиция в контейнере откуда будут загружаться данные.
    
    //explicit AbstractSurface(const geo::PrimitiveParameters &p_parameters, unsigned &p_n);

    /// \brief Конструктор копирования
    AbstractSurface(const AbstractSurface &p_surface) = default;

    /// \brief Деструктор
    virtual ~AbstractSurface();

    /// \brief Оператор присваивания
    AbstractSurface& operator=(const AbstractSurface &p_surface) = default;

    /// \brief Вычисление положения точки на поверхности
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Точка на поверхности
    virtual void getPoint(double p_param_1, double p_param_2, double *o_point) const = 0;

    /// \brief Метод для вычисления производной поверхности в данной точке по первому параметру
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_DerPoint Значение производной в точке на поверхности
    virtual void getDer1Point(double p_param_1, double p_param_2, double *o_der_1_point) const = 0;

    /// \brief Метод для вычисления производной поверхности в данной точке по второму параметру
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_DerPoint Значение производной в точке на поверхности
    virtual void getDer2Point(double p_param_1, double p_param_2, double *o_der_2_point) const = 0;

    /// \brief Метод для вычисления нормали поверхности в данной точке
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_norm_point Значение нормали в точке на поверхности
   // virtual void getNormalPoint(double p_param_1, double p_param_2, double *o_norm_point) const;

    /// \brief Метод для вычисления (k,l)-производных от 0-го до d-ого порядка (0 <= k + l <= d).
    /// \param [in] p_der_order Значение d.
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_point Значения производных по параметру в точке на поверхности
    virtual void getPointAndDerivs(unsigned p_der_order,
                                   double p_param_1,
                                   double p_param_2,
                                   double *o_der_point) const = 0;

    /// \brief Метод для вычисления производной нормали поверхности в данной точке
    /// \param [in] p_der_order Значение d.
    /// \param [in] p_param_1 Значение первого параметра
    /// \param [in] p_param_2 Значение второго параметра
    /// \param [out] o_norm_point Значение производной нормали в точке на поверхности
   /* virtual void getDerNormalPoint(unsigned p_der_order,
                                   double   p_param_1,
                                   double   p_param_2,
                                   double   *o_der_norm_point) const;

                                   */
    /// \brief Возвращение начальных значений параметров поверхности
    /// \return Начальные значения параметров поверхности
    std::pair<double, double> getBegin() const;

    /// \brief Возвращение конечных значений параметров поверхности
    /// \return Конечные значения параметров поверхности
    std::pair<double, double> getEnd() const;

    /// \brief Возвращение идентификатора поверхности
    unsigned getId() const;

    /// \brief TRUE, если поверхность замкнута по первому параметру, иначе - FALSE.
    bool isClosed1() const;

    /// \brief TRUE, если поверхность замкнута по второму параметру, иначе - FALSE.
    bool isClosed2() const;

    /// \brief Метод проверяет замкнутость поверхности по первому параметру.
   // bool checkClosed1() const;

    /// \brief Метод проверяет замкнутость поверхности по второму параметру.
   // bool checkClosed2() const;

    /// \brief Преверка вырожденности по первому параметру при u = u0
    bool checkDegenerateBoundary1Begin() const;

    /// \brief Преверка вырожденности по первому параметру при u = u1
    bool checkDegenerateBoundary1End() const;

    /// \brief Преверка вырожденности по второму параметру при v = v0
    bool checkDegenerateBoundary2Begin() const;

    /// \brief Преверка вырожденности по второму параметру при v = v1
    bool checkDegenerateBoundary2End() const;

    /// \brief Нахождение начального приближения проекции точки
    /// \param [in] p_point Координаты проецируемой точки
    /// \param [out] o_point_parameters Параметры начального приближения проекции
    virtual void findInitProjectPoint(
        const double *p_point,
        double       *o_point_parameters) const;

    /// \brief Проецирование точки на поверхность
    /// \param [in] p_point Координаты проецируемой точки
    /// \param [out] o_projected_point Координаты проекции на поверхность
    /// \param [in,out] o_point_parameters На входе начальное приближение, на выходе параметрические координаты проекции.
    /// \return 0 или код ошибки
    virtual int projectPoint(
        const double *p_point,
        double       *o_projected_point,
        double       *o_point_parameters) const;

    /// \brief Метод возвращает поверхность, отражённую относительно оси 2 (ось 1 инвертируется)
    virtual std::shared_ptr<AbstractSurface> reflect2Surface() const;

    /// \brief Возвращение значения размерности пространства
    unsigned getDim() const;

    /// \brief Компиляторонезависимая реализация type_info.name()
    /*virtual const std::string& type()
    {
        return m_type;
    }
    static inline const std::string m_type = "AbstractSurface";
    */

    /// \brief Рассчитывает границы поверхности
    virtual void calculateBoundaryCurves() = 0;

    /// \brief Возвращает левую границу поверхности
    std::shared_ptr<AbstractCurve> getBoundary1Begin() const;

    /// \brief Возвращает правую границу поверхности
    std::shared_ptr<AbstractCurve> getBoundary1End() const;

    /// \brief Возвращает верхнюю границу поверхности
    std::shared_ptr<AbstractCurve> getBoundary2Begin() const;

    /// \brief Возвращает нижнюю границу поверхности
    std::shared_ptr<AbstractCurve> getBoundary2End() const;

    /// \brief Возвращает полюса поверхности
    const std::map<PoleKeyType, PoleValType> getPoles() const;

    /// \brief Устанавливает полюс поверхности
    /// \param [in] p_nondegenerated_dir Направление по которому нет вырождения (1 или 2)
    /// \param [in] p_nondegenerated_param Значение невырожденного параметра
    /// \param [in] p_degenerated_param_1 Значение первого вырожденного параметра для расчета нормали
    /// \param [in] p_degenerated_param_2 Значение второго вырожденного параметра для расчета нормали
    void setPole(unsigned p_nondegenerated_dir,
                 double p_nondegenerated_param,
                 double p_degenerated_param_1,
                 double p_degenerated_param_2);

protected:
    unsigned m_id;
    unsigned m_dim;
    double   m_begin_1;
    double   m_end_1;
    double   m_begin_2;
    double   m_end_2;
    bool     m_is_closed_1;
    bool     m_is_closed_2;

    std::shared_ptr<AbstractCurve> m_boundary_begin_1;
    std::shared_ptr<AbstractCurve> m_boundary_end_1;
    std::shared_ptr<AbstractCurve> m_boundary_begin_2;
    std::shared_ptr<AbstractCurve> m_boundary_end_2;

    double *m_r1_array;
    double *m_r2_array;
    unsigned *m_counts_array;

    unsigned m_max_der_order;

private:
    void initializeDependentData();
    std::vector<double> m_r1;
    std::vector<double> m_r2;
    std::vector<unsigned> m_counts;
    std::vector<unsigned> m_counts_base;
    std::vector<double> m_der_point1;
    std::vector<double> m_der_point2;

    unsigned *m_counts_base_array;
    double *m_der_point1_array;
    double *m_der_point2_array;

    // полюса поверхности
    std::map<PoleKeyType, PoleValType> m_poles;
};


#endif // SURFACE_H
