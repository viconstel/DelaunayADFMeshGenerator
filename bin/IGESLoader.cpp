#include "IGESLoader.h"
#include "AbstarctIGESParamsParsers.h"
#include <algorithm>
#include <fstream>
#include <sstream>
//#include <Log>
//#include <spdlog/spdlog.h>


#define IGES_LINE_LENGTH 81


int IGESLoader::loadGeo()
{
    //INFO("Start of reading geometry IGES-file");

    std::vector<std::string> file_content;
    char line_buf[IGES_LINE_LENGTH];
    std::string field;
    unsigned lines_count{ 0 };
    unsigned s_lines_count{ 0 };
    unsigned g_lines_count{ 0 };
    unsigned de_lines_count{ 0 };
    unsigned pd_lines_count{ 0 };

    // Чтение файла
    const std::string &geo_file_name = getGeoFileName();
    std::ifstream in(geo_file_name);
    if (!in)
    {
       // ERR("Geometry file {} isn't found", geo_file_name.c_str());
        return -1;
    }

    while (in)
    {
        in.getline(line_buf, std::numeric_limits<std::streamsize>::max());
        if (0 != std::strlen(line_buf))
            file_content.push_back(std::string(line_buf, line_buf + IGES_LINE_LENGTH));
    }

    // Определение количества строк в секциях
    lines_count = file_content.size();
    field = file_content[lines_count - 1].substr(1, 7);
    s_lines_count = std::stoi(field);
    field = file_content[lines_count - 1].substr(9, 7);
    g_lines_count = std::stoi(field);
    field = file_content[lines_count - 1].substr(17, 7);
    de_lines_count = std::stoi(field);
    field = file_content[lines_count - 1].substr(25, 7);
    pd_lines_count = std::stoi(field);

    // Чтение данных секции Data Entry (DE)
    DataEntrySections de_sections;
    unsigned beg = s_lines_count + g_lines_count;
    unsigned end = s_lines_count + g_lines_count + de_lines_count;
    unsigned entity_type;
    unsigned sequence_number{ 1 };

    auto checkField = [](const std::string &p_field_str)
    {
        return p_field_str.find_first_not_of(" ") != p_field_str.npos;
    };

    for (int i = beg; i < end; i += 2)
    {
        field = file_content[i].substr(0, 8);
        entity_type = checkField(field) ? std::stoi(field) : 0;

        field = file_content[i].substr(8, 8);
        de_sections[entity_type][sequence_number].m_pd_pointer = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(16, 8);
        de_sections[entity_type][sequence_number].m_structure = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(24, 8);
        de_sections[entity_type][sequence_number].m_line_font_pattern = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(32, 8);
        de_sections[entity_type][sequence_number].m_level = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(40, 8);
        de_sections[entity_type][sequence_number].m_view = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(48, 8);
        de_sections[entity_type][sequence_number].m_transformation_matrix_pointer = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(56, 8);
        de_sections[entity_type][sequence_number].m_label_display_associativity = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i].substr(64, 8);
        de_sections[entity_type][sequence_number].m_status_number = field;
        field = file_content[i].substr(72, 8);

        field = file_content[i + 1].substr(0, 8);
        field = file_content[i + 1].substr(8, 8);
        de_sections[entity_type][sequence_number].m_line_weight_number = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i + 1].substr(16, 8);
        de_sections[entity_type][sequence_number].m_color_number = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i + 1].substr(24, 8);
        de_sections[entity_type][sequence_number].m_parameter_line_count = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i + 1].substr(32, 8);
        de_sections[entity_type][sequence_number].m_form_number = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i + 1].substr(40, 8);
        field = file_content[i + 1].substr(48, 8);
        field = file_content[i + 1].substr(56, 8);
        de_sections[entity_type][sequence_number].m_entity_label = field;
        field = file_content[i + 1].substr(64, 8);
        de_sections[entity_type][sequence_number].m_entity_subscript_number = checkField(field) ? std::stoi(field) : 0;
        field = file_content[i + 1].substr(72, 8);

        sequence_number += 2;
    }
    
    // Чтение данных секции Parameter Data (PD)
    ParameterDataSections pd_sections;
    beg = s_lines_count + g_lines_count + de_lines_count;
    end = s_lines_count + g_lines_count + de_lines_count + pd_lines_count;
    std::string params_str;
    std::shared_ptr<AbstarctIGESParamsParsers> parser;

    int i = beg;
    while (i < end)
    {
        field = file_content[i].substr(73, 8);
        sequence_number = std::stoi(field);
        
        field = file_content[i].substr(64, 8);
        pd_sections[sequence_number].m_de_pointer = std::stoi(field);

        field = file_content[i].substr(0, 3);
        pd_sections[sequence_number].m_entity_type = std::stoi(field);

        params_str = file_content[i].substr(4, 60);
        while (params_str.find(";") == params_str.npos)
            params_str += file_content[++i].substr(0, 64);

        params_str.erase(std::remove(params_str.begin(), params_str.end(), ' '), params_str.end());
        params_str.erase(std::remove(params_str.begin(), params_str.end(), ';'), params_str.end());

        parser = AbstarctIGESParamsParsers::Make(pd_sections[sequence_number].m_entity_type);
        if (nullptr != parser)
            pd_sections[sequence_number].m_params = parser->Parse(params_str);
        //else
           // WARN("The unsupported entity type in the IGES-file is detected: {}",
                 //pd_sections[sequence_number].m_entity_type);

        ++i;
    }

    m_iges_data.init(std::move(de_sections), std::move(pd_sections));

  //  INFO("Reading the geometry IGES-file is complete successfully");
    return 0;
}


/*std::shared_ptr<AbstractGeometryModel> IGESLoader::getGeoModel(GeoModels p_geo_model)
{
    // TODO: Добавить поддержку твердотельных объектов (тип 186)

    std::shared_ptr<ParametricGeometryModel> new_model;

    if (p_geo_model != GeoModels::parametric_geometry_model)
    {
        WARN("Geometry file format does not support this geometry model");
        return nullptr;
    }

    const IGESData::SurfacesBoundaryData &surfaces_boundary_data = m_iges_data.getSurfacesBoundaryData();
    const IGESData::NURBSSurfaces &nurbs_surfaces = m_iges_data.getNURBSSurfaces();
    const IGESData::NURBSCurves &nurbs_curves = m_iges_data.getNURBSCurves();

    Body new_body;
    for (const auto &surface_boundary_data : surfaces_boundary_data) // Обход поверхностей
    {
        auto surface_iter = nurbs_surfaces.find(surface_boundary_data.first);
        std::shared_ptr<AbstractSurface> nurbs_surface(new NURBSSurface(surface_iter->second));
        std::shared_ptr<BoundedSurface> new_bounded_surface = std::make_shared<BoundedSurface>();
        new_bounded_surface->setSurface(nurbs_surface);
        new_bounded_surface->reserve(surface_boundary_data.second.size());

        for (const auto &boundary : surface_boundary_data.second) // Обход границ
        {
            SurfaceBoundary new_surface_boundary(boundary.second.size());

            for (const auto &curve : boundary.second) // Обход кривых границы
            {
                auto curve_iter = nurbs_curves.find(curve.first);
                std::shared_ptr<AbstractCurve> nurbs_curve(new NURBSCurve(curve_iter->second));
                new_surface_boundary.addCurve(curve.first, curve.second, nurbs_curve);
            }

            new_bounded_surface->addBoundary(boundary.first, new_surface_boundary);
        }

        new_body.addBoundedSurface(surface_boundary_data.first,
                                   new_bounded_surface,
                                   false);
    }

    new_model = std::make_shared<ParametricGeometryModel>();
    new_model->addBody(1, new_body);

    return new_model;
}*/


IGESData& IGESLoader::GetIGESModel()
{
    return m_iges_data;
}
