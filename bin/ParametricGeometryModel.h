#ifndef GEOMETRY_MODEL_H
#define GEOMETRY_MODEL_H

//#include <AbstractGeometryModel>
#include "AbstractCurve.h"
#include "AbstractSurface.h"

#include <vector>
#include <unordered_map>
#include <string>
#include <memory>
#include <algorithm>


/// \brief Класс исключения
class EntityNotFoundException : public std::exception
{
public:
    EntityNotFoundException(const std::string &p_msg)
        : m_msg(p_msg)
    {}

    char const* what() const noexcept override
    {
        return m_msg.c_str();
    }

private:
    std::string m_msg;
};


/// \brief Класс границы поверхности
class SurfaceBoundary
{
public:
    /// \brief Формат элементов контейнера: <ID кривой, <Указатель на кривую, Флаг обхода>>
    /// Для флага: TRUE - направление обхода кривой измениять не следует,
    /// FALSE - напрвление обхода следует изменить.
    using OrientedCurve = std::pair<std::shared_ptr<AbstractCurve>, bool>;
    using SurfaceBoundaryContType = std::vector<std::pair<unsigned, OrientedCurve>>;

public:
    /// \brief Конструктор выделяет память под контейнер кривых
    SurfaceBoundary(unsigned p_curve_count);

    /// \brief Добавление кривой в конец контейнера кривых
    void addCurve(
        unsigned                       p_curve_id,
        bool                           p_is_direction_change,
        std::shared_ptr<AbstractCurve> p_curve);

    /// \brief Получение кривой по ID
    OrientedCurve getCurve(
        unsigned p_curve_id);

    /// \brief Получение списка всех кривых данной границы
    const SurfaceBoundaryContType& getCurves() const;

    /// \brief Получение списка всех кривых данной границы
    SurfaceBoundaryContType& getCurves();

private:
    SurfaceBoundaryContType m_surface_boundary; ///< Граница поверхности
};


/// \brief Класс ограниченной поверхности
class BoundedSurface
{
public:
    /// \brief Формат элементов контейнера: <ID границы, Граница>
    using BoundedSurfaceContType = std::vector<std::pair<unsigned, SurfaceBoundary>>;

public:
    BoundedSurface() = default;

    BoundedSurface(
        std::shared_ptr<AbstractSurface> p_surface,
        unsigned                         p_boundary_count);

    /// \brief Задание NURBS-поверхности
    void setSurface(
        std::shared_ptr<AbstractSurface> p_surface);

    /// \brief Выделение памяти под контейнер границ
    void reserve(unsigned p_boundary_count);

    /// \brief Получение NURBS-поверхности
    const std::shared_ptr<AbstractSurface> getSurface() const;

    /// \brief Получение NURBS-поверхности
    std::shared_ptr<AbstractSurface> getSurface();

    /// \brief Добавление границы в конец контейнера
    void addBoundary(
        unsigned               p_boundary_id,
        const SurfaceBoundary &p_boundary);

    /// \brief Получение границы по ID
    const SurfaceBoundary& getBoundary(
        unsigned p_boundary_id);

    /// \brief Получение списка всех границ для данной поверхности
    const BoundedSurfaceContType& getBoundaries() const;

    /// \brief Получение списка всех границ для данной поверхности
    BoundedSurfaceContType& getBoundaries();

private:
    std::shared_ptr<AbstractSurface> m_surface;            ///< Параметрическая поверхность
    BoundedSurfaceContType           m_surface_boundaries; ///< Границы поверхности
};


/// \brief Класс трехмерной области
class Body
{
public:
    /// \brief Формат элементов контейнера: <ID тела, <Указатель на ограниченную поверхность, Флаг ориентации>>.
    /// Для флага: TRUE - направление обхода границ задает внутреннюю нормаль к гране тела, FALSE - внешню.
    using OrientedBoundedSurface =
        std::pair<std::shared_ptr<BoundedSurface>, bool>;
    using BodyContType =
        std::unordered_map<unsigned, OrientedBoundedSurface>;

public:
    /// \brief Добавление ограниченной поверхности
    bool addBoundedSurface(
        unsigned                        p_bounded_surface_id,
        std::shared_ptr<BoundedSurface> p_bounded_surface,
        bool                            p_orientation_flag);

    /// \brief Получение ограниченной поверхности по ее ID
    const OrientedBoundedSurface& getBoundedSurface(unsigned p_bounded_surface_id) const;

    /// \brief Получение ограниченных поверхностей для данного тела
    const BodyContType& getBoundedSurfaces() const;

    /// \brief Получение ограниченных поверхностей для данного тела
    BodyContType& getBoundedSurfaces();

private:
    BodyContType m_body; ///< Тело
};


/// \brief Класс парметрической геометрической модели
/*class ParametricGeometryModel : public AbstractGeometryModel
{
public:
    /// \brief Формат элементов контейнера: <ID тела, Тело>
    using ParametricGeometryModelContType = std::unordered_map<unsigned, Body>;

    /// \brief Формат элементов контейнера: <ID кривой, Кривая>
    using ModelCurvesContainer =
        std::unordered_map<unsigned,
                           std::shared_ptr<AbstractCurve>>;

    /// \brief Формат элементов контейнера: <ID грани, Грань>
    using ModelSurfacesContainer =
        std::unordered_map<unsigned,
                           std::shared_ptr<BoundedSurface>>;
    
public:
    /// \brief Добавление тела
    bool addBody(unsigned p_body_id, const Body &p_body);

    /// \brief Получение тела по его ID
    const Body& getBody(unsigned p_body_id);

    /// \brief Получение всех тел
    const ParametricGeometryModelContType& getBodies() const;

    /// \brief Получение всех тел
    ParametricGeometryModelContType& getBodies();

    /// \brief Получение контейнера всех кривых
    ModelCurvesContainer getAllModelCurves();
    
    /// \brief Получение контейнера всех ограниченных поверхностей
    ModelSurfacesContainer getAllModelBoundedSurfaces();

private:
    ParametricGeometryModelContType m_bodies; ///< Модель
};
*/
#endif // GEOMETRY_MODEL_H
