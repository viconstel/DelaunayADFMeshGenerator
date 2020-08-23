#include "IGESData.h"
#include <algorithm>


void
IGESData::
init(
    DataEntrySections     &&p_de_sections,
    ParameterDataSections &&p_pd_sections)
{
    reset();

    m_de_sections = std::move(p_de_sections);
    m_pd_sections = std::move(p_pd_sections);
}


void
IGESData::
reset()
{
    m_de_sections.clear();
    m_pd_sections.clear();

    m_nurbs_curves_cache.clear();
    m_nurbs_surfaces_cache.clear();
    m_surfaces_boundary_data_cache.clear();
}


const IGESData::NURBSCurves&
IGESData::
getNURBSCurves()
{
    if (!m_nurbs_curves_cache.empty())
        return m_nurbs_curves_cache;

    auto pd_iter = m_pd_sections.cbegin();
    while (true)
    {
        pd_iter = std::find_if(pd_iter, m_pd_sections.cend(),
            [](const ParameterDataSections::value_type &p_pd) {
                return p_pd.second.m_entity_type == 126;
            });
        
        if (pd_iter != m_pd_sections.end())
        {
            auto     &pd_params          = pd_iter->second;
            unsigned  degree             = pd_params["M"].toInteger();
            unsigned  upper_index_of_sum = pd_params["K"].toInteger();
            unsigned  knots_count        = (degree + 1) + (upper_index_of_sum + 1);

            std::vector<double> knots_vector(knots_count);
            std::vector<double> weights(upper_index_of_sum + 1);
            std::vector<double> control_points(3 * (upper_index_of_sum + 1));

            for (int i = 0; i < knots_count; ++i)
            {
                auto res = pd_params["T(" + std::to_string(-static_cast<int>(degree) + i) + ")"];
                double val = res.toReal();
                knots_vector[i] = val;
            }

            for (int i = 0; i < upper_index_of_sum + 1; ++i)
                weights[i] = pd_params["W(" + std::to_string(i) + ")"].toReal();

#           define _(point, coord) (3 * (point) + (coord))
            for (int i = 0; i < upper_index_of_sum + 1; ++i)
            {
                control_points[_(i, 0)] = pd_params["X(" + std::to_string(i) + ")"].toReal();
                control_points[_(i, 1)] = pd_params["Y(" + std::to_string(i) + ")"].toReal();
                control_points[_(i, 2)] = pd_params["Z(" + std::to_string(i) + ")"].toReal();
            }
#           undef _
            
            std::vector<double> normal(3);
            normal[0] = pd_params["XNORM"].toReal();
            normal[1] = pd_params["YNORM"].toReal();
            normal[2] = pd_params["ZNORM"].toReal();
            
            m_nurbs_curves_cache.insert({
                static_cast<unsigned>(pd_iter->first),
                {
                    3,
                    degree,
                    knots_vector,
                    weights,
                    control_points,
                    normal,
                    pd_params["V(0)"].toReal(),
                    pd_params["V(1)"].toReal(),
                    static_cast<bool>(pd_params["PROP2"].toInteger())
                }
            });

            ++pd_iter;
        }
        else break;
    }

    return m_nurbs_curves_cache;
}


const IGESData::NURBSSurfaces&
IGESData::
getNURBSSurfaces()
{
    if (!m_nurbs_surfaces_cache.empty())
        return m_nurbs_surfaces_cache;

    auto pd_iter = m_pd_sections.cbegin();
    while (true)
    {
        pd_iter = std::find_if(pd_iter, m_pd_sections.cend(),
            [](const ParameterDataSections::value_type &p_pd) {
            return p_pd.second.m_entity_type == 128;
        });

        if (pd_iter != m_pd_sections.end())
        {
            auto     &pd_params = pd_iter->second;
            unsigned  degree_1 = pd_params["M1"].toInteger();
            unsigned  degree_2 = pd_params["M2"].toInteger();
            unsigned  upper_index_of_sum_1 = pd_params["K1"].toInteger();
            unsigned  upper_index_of_sum_2 = pd_params["K2"].toInteger();
            unsigned  knots_count_1 = (degree_1 + 1) + (upper_index_of_sum_1 + 1);
            unsigned  knots_count_2 = (degree_2 + 1) + (upper_index_of_sum_2 + 1);
            unsigned  weights_count = (upper_index_of_sum_1 + 1) * (upper_index_of_sum_2 + 1);

            std::vector<double> knots_vector_1(knots_count_1);
            std::vector<double> knots_vector_2(knots_count_2);
            std::vector<double> weights(weights_count);
            std::vector<double> control_points(3 * weights_count);

            for (int i = 0; i < knots_count_1; ++i)
            {
                auto res = pd_params["S(" + std::to_string(-static_cast<int>(degree_1) + i) + ")"];
                knots_vector_1[i] = res.toReal();
            }

            for (int i = 0; i < knots_count_2; ++i)
            {
                auto res = pd_params["T(" + std::to_string(-static_cast<int>(degree_2) + i) + ")"];
                knots_vector_2[i] = res.toReal();
            }

#           define _(i, j) ((i) + (upper_index_of_sum_1 + 1) * (j))
            for (int i = 0; i < upper_index_of_sum_1 + 1; ++i)
                for (int j = 0; j < upper_index_of_sum_2 + 1; ++j)
                    weights[_(i, j)] = pd_params["W(" + std::to_string(i) + "," + std::to_string(j) + ")"].toReal();
#           undef _
#           define _(i, j, coord) ((coord) + 3 * ((i) + (upper_index_of_sum_1 + 1) * (j)))
            for (int i = 0; i < upper_index_of_sum_1 + 1; ++i)
                for (int j = 0; j < upper_index_of_sum_2 + 1; ++j)
                {
                    control_points[_(i, j, 0)] = pd_params["X(" + std::to_string(i) + "," + std::to_string(j) + ")"].toReal();
                    control_points[_(i, j, 1)] = pd_params["Y(" + std::to_string(i) + "," + std::to_string(j) + ")"].toReal();
                    control_points[_(i, j, 2)] = pd_params["Z(" + std::to_string(i) + "," + std::to_string(j) + ")"].toReal();
                }
#           undef _

            m_nurbs_surfaces_cache.insert({
                static_cast<unsigned>(pd_iter->first),
                {
                    3,
                    upper_index_of_sum_1,
                    upper_index_of_sum_2,
                    degree_1,
                    degree_2,
                    knots_vector_1,
                    knots_vector_2,
                    weights,
                    control_points,
                    pd_params["U(0)"].toReal(),
                    pd_params["U(1)"].toReal(),
                    pd_params["V(0)"].toReal(),
                    pd_params["V(1)"].toReal(),
                    static_cast<bool>(pd_params["PROP1"].toInteger()),
                    static_cast<bool>(pd_params["PROP2"].toInteger())
                } 
            });

            ++pd_iter;
        }
        else break;
    }

    return m_nurbs_surfaces_cache;
}


const IGESData::SurfacesBoundaryData&
IGESData::
getSurfacesBoundaryData()
{
    if (!m_surfaces_boundary_data_cache.empty())
        return m_surfaces_boundary_data_cache;
    
    // Получение всех NURBS-кривых
    const NURBSCurves &nurbs_curves = getNURBSCurves();

    // Определение совпадающих кривых
    MatchingCurvesGroups curve_groups;
    getMatchingCurves(nurbs_curves, curve_groups);

    // Получение данных о всех границах (тип 143)
    const auto &de_boundaries = m_de_sections.find(143)->second;

    // Обход всех граничных данных (тип 143)
    auto de_boundaries_iter = de_boundaries.begin();
    auto de_boundaries_iter_end = de_boundaries.end();
    for (; de_boundaries_iter != de_boundaries_iter_end; ++de_boundaries_iter)
    {
        // Получение всех граничных данных (тип 143) для данной поверхности
        const auto &pd_surface_boundaries_data = m_pd_sections.find(de_boundaries_iter->second.m_pd_pointer)->second;

        // Указатель на DE-секцию для данной поверхности
        unsigned surface_de_pointer = pd_surface_boundaries_data["SPTR"].toInteger();

        // Определение PD идентификатора поверхности
        const auto &de_nurbs_surfaces = m_de_sections.find(128)->second;
        unsigned surface_pd_pointer = std::find_if(
            de_nurbs_surfaces.begin(),
            de_nurbs_surfaces.end(),
            [&surface_de_pointer](const auto &p_data_entry_section) {
                return p_data_entry_section.first == surface_de_pointer;
            })->second.m_pd_pointer;

        // Определение количества границ у данной поверхности
        unsigned boundaries_count = pd_surface_boundaries_data["N"].toInteger();
        SurfaceBoundaries surface_boundaries;

        // Для каждой границы (тип 141) данной поверхности
        for (unsigned i = 0; i < boundaries_count; ++i)
        {
            // Указатель на DE-секцию (тип 141) текущей границы
            unsigned de_boundary_pointer = pd_surface_boundaries_data["BDPT(" + std::to_string(i + 1) + ")"].toInteger();

            // Определение PD идентификатора текущей границы (тип 141)
            const auto &de_boundaries = m_de_sections.find(141)->second;
            unsigned pd_boundary_pointer = std::find_if(
                de_boundaries.begin(),
                de_boundaries.end(),
                [&de_boundary_pointer](const auto &p_data_entry_section) {
                    return p_data_entry_section.first == de_boundary_pointer;
                })->second.m_pd_pointer;

            // Получение граничных данных (тип 141) данной границы
            const auto &pd_boundary_data = m_pd_sections.find(pd_boundary_pointer)->second;
            
            // Определение количества кривых, составляющих данную границу
            unsigned boundary_curves_count = pd_boundary_data["N"].toInteger();
            std::vector<std::pair<unsigned, bool>> surface_boundary_curves(boundary_curves_count);

            for (unsigned j = 0; j < boundary_curves_count; ++j)
            {
                // Указатель на DE-секцию для данной кривой
                unsigned de_curve_pointer = pd_boundary_data["CRVPT(" + std::to_string(j + 1) + ")"].toInteger();

                // Поиск PD идентификатора NURBS-кривой
                const auto &de_curves = m_de_sections.find(126)->second;
                unsigned pd_curve_pointer = std::find_if(
                    de_curves.begin(),
                    de_curves.end(),
                    [&de_curve_pointer](const auto &data_entry_section) {
                        return data_entry_section.first == de_curve_pointer;
                    })->second.m_pd_pointer;

                unsigned is_change_bypass_direction = 
                    pd_boundary_data["SENSE(" + std::to_string(j + 1) + ")"].toInteger();

                // Учет полностью совпадающих кривых
                const auto &group = *curve_groups.find(pd_curve_pointer);
                if (!group.second.empty())
                {
                    const MatchingCurvesGroup::value_type &min_pair = *std::min_element(
                        group.second.begin(),
                        group.second.end(),
                        [](const auto &p_val_1, const auto &p_val_2) {
                            return p_val_1.first < p_val_2.first;
                        });

                    unsigned id = std::min(pd_curve_pointer, min_pair.first);

                    if (id != pd_curve_pointer)
                    {
                        surface_boundary_curves[j].first = id;
                        
                        // Исходная кривая - сохранить; новая кривая - совпадает с исходной
                        if (is_change_bypass_direction == 1 &&
                            min_pair.second == static_cast<uint8_t>(NURBSCurve::forward_equal))
                            surface_boundary_curves[j].second = false;
                        // Исходная кривая - сохранить; новая кривая - совпадает с исходной при инвертации
                        else if (is_change_bypass_direction == 1 &&
                                 min_pair.second == static_cast<uint8_t>(NURBSCurve::backward_equal))
                            surface_boundary_curves[j].second = true;
                        // Исходная кривая - изменить; новая кривая - совпадает с исходной
                        else if (is_change_bypass_direction == 2 &&
                                 min_pair.second == static_cast<uint8_t>(NURBSCurve::forward_equal))
                            surface_boundary_curves[j].second = true;
                        // Исходная кривая - изменить; новая кривая - совпадает с исходной при инвертации
                        else if (is_change_bypass_direction == 2 &&
                                 min_pair.second == static_cast<uint8_t>(NURBSCurve::backward_equal))
                            surface_boundary_curves[j].second = false;

                        continue;
                    }
                }

                surface_boundary_curves[j].first = pd_curve_pointer;

                if (is_change_bypass_direction == 1)
                    surface_boundary_curves[j].second = false;
                else if (is_change_bypass_direction == 2)
                    surface_boundary_curves[j].second = true;
            }

            // Добавление данных о кривой в список кривых, формирующих границу
            surface_boundaries.emplace(pd_boundary_pointer, std::move(surface_boundary_curves));
        }

        // Добавление данных о границах поверхности в общий контейнер
        m_surfaces_boundary_data_cache.emplace(surface_pd_pointer, surface_boundaries);
    }

    return m_surfaces_boundary_data_cache;
}
