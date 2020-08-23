#ifndef IGES_MODEL
#define IGES_MODEL

// std: c++ 17

//#include "IGESData.h"
//#include <unt_StdDefines>
#include "IGESDataTypes.h"
#include "ParametricGeometryModel.h"
#include "NURBSCurve.h"
#include "NURBSSurface.h"
#include <vector>


/// \brief Класс данных IGES-файла
class IGESData
{
public:
    /// \brief Котейнер произвольных сущностей IGES-файла с ID
    template<class _EntityType>
    using IGESMap = std::map<unsigned, _EntityType>;

    /// \brief Тип списка NURBS-кривых модели.
    /// Ключ - Sequence Number
    using NURBSCurves = IGESMap<NURBSCurve>;

    /// \brief Тип списка NURBS-поверхностей модели.
    /// Ключ - Sequence Number
    using NURBSSurfaces = IGESMap<NURBSSurface>;

    /// \brief Список кривых границы с указанием направления обхода
    /// First - ID кривой, Second - флаг изменения направления обхода
    using SurfaceBoundaryCurves = std::vector<std::pair<unsigned, bool>>;

    /// \brief Список границ данной поверхности
    /// Ключ - ID границы, значение - список кривых, формирующих границу.
    using SurfaceBoundaries = IGESMap<SurfaceBoundaryCurves>;

    /// \brief Тип данных для хранения информации о границах поверхностей и последовательности их обхода
    /// Ключ - ID поверхности, значение - список границ.
    using SurfacesBoundaryData = IGESMap<SurfaceBoundaries>;

    /// \brief Тип данных для задания множества кривых, совпадающей с данной
    /// \detail Формат: <ID кривой, код равенства>
    using MatchingCurvesGroup = std::vector<std::pair<unsigned, uint8_t>>;

    /// \brief Контейнер, задающий множество кривых, совпадающих с данной
    /// \detail Формат: <ID кривой, множество кривых, совпадающих с данной>
    using MatchingCurvesGroups = std::unordered_map<unsigned, MatchingCurvesGroup>;

 public:
    /// \brief Метод инициализации модели
    /// \param [in] p_DESections Список секций Data Entry (DE).
    /// \param [in] p_PDSections Список секций Parameter Data (PD).
    void init(DataEntrySections &&p_de_sections, ParameterDataSections &&p_pd_sections);

    /// \brief Метод освобождает внутренние ресурсы
    void reset();

    /// \brief Метод создает список NURBS-кривых (тип 126) модели по данным секции PD
    /// \return Список NURBS-кривых модели
    const NURBSCurves& getNURBSCurves();

    /// \brief Метод создает список NURBS-поверхностей (тип 128) модели по данным секции PD
    /// \return Список NURBS-кривых модели
    const NURBSSurfaces& getNURBSSurfaces();

    /// \brief Метод возввращает данные о граничных кривых поверхностей.
    const SurfacesBoundaryData& getSurfacesBoundaryData();

    /// \brief Метод создает контейнер с информацией о совпадающих кривых
    /// \param [in] p_curves Кривые, для которых требуется провести анализ
    /// \param [out] o_groups Контейнер групп
    template<class _CurveType>
    void getMatchingCurves(
        const IGESMap<_CurveType> &p_curves,
        MatchingCurvesGroups      &o_groups) const;

private:
    DataEntrySections     m_de_sections; ///< Текстовые данные секции Data Entry (DE)
    ParameterDataSections m_pd_sections; ///< Текстовые данные секции Parameter Data (PD)

    NURBSCurves          m_nurbs_curves_cache;           ///< Кэш NURBS-кривых
    NURBSSurfaces        m_nurbs_surfaces_cache;         ///< Кэш NURBS-поверхностей
    SurfacesBoundaryData m_surfaces_boundary_data_cache; ///< Кэш граничных данных
};


template<class _CurveType>
void IGESData::getMatchingCurves(
     const IGESMap<_CurveType> &p_curves,
     MatchingCurvesGroups      &o_groups) const
{   
    const unsigned curves_count = p_curves.size();
    
    if (!o_groups.empty())
        o_groups.clear();
    o_groups.reserve(curves_count);
        
    std::vector<unsigned> curve_ids;
    curve_ids.reserve(curves_count);

    for (const auto &[id, nurbs_curve] : p_curves)
        curve_ids.push_back(id);
    
#   pragma omp paralel for
    for (unsigned i = 0; i < curves_count; ++i)
    {
        const _CurveType &curve_1 = p_curves.find(curve_ids[i])->second;

        MatchingCurvesGroup new_group;

        for (unsigned j = 0; j < curves_count; ++j)
        {
            if (i == j)
                continue;

            const _CurveType &curve_2 = p_curves.find(curve_ids[j])->second;

            if (uint8_t code = (curve_1 == curve_2); code != 0)
                new_group.push_back({ curve_ids[j], code });
        }

#       pragma omp critical
        o_groups.insert({ curve_ids[i], std::move(new_group) });
    }
}

#endif // IGES_MODEL
