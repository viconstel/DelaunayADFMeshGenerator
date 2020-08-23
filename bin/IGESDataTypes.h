#ifndef UNT_IGES_DATA_TYPES
#define UNT_IGES_DATA_TYPES

#include <string>
#include <map>
#include <stdexcept>


/// Типы данных IGES
enum class IGESDataTypes
{
    INTEGER,
    REAL,
    STRING,
    POINTER,
    LANGUAGE_STATEMENT,
    LOGICAL
};


/// Тип свойств параметра Parameter Data (PD)
struct PDParameter
{
    IGESDataTypes m_type;
    std::string   m_str_value;

    int toInteger() const
    {
        return std::stoi(m_str_value);
    }

    double toReal() const
    {
        return std::stod(m_str_value);
    }

    int toPointer() const
    {
        return std::stoi(m_str_value);
    }

    bool toLogical() const
    {
        return std::stoi(m_str_value);
    }
};


/// Тип секции Data Entry (DE)
struct DataEntrySection
{
    int         m_pd_pointer;
    int         m_structure;
    int         m_line_font_pattern;
    int         m_level;
    int         m_view;
    int         m_transformation_matrix_pointer;
    int         m_label_display_associativity;
    std::string m_status_number;
    char        m_line_weight_number;
    char        m_color_number;
    int         m_parameter_line_count;
    int         m_form_number;
    std::string m_entity_label;
    int         m_entity_subscript_number;
};


/// Тип перечня секций Data Entry (DE).
/// Ключи: Entity Type, Sequence Number
using DataEntrySections = std::map<int, std::map<int, DataEntrySection>>;


/// Тип секции Parameter Data (PD)
struct stc_ParameterDataSection
{
    int                                m_entity_type;
    std::map<std::string, PDParameter> m_params;
    int                                m_de_pointer;

    PDParameter& operator[](const std::string &p_param_sid)
    {
        return m_params.find(p_param_sid)->second;
    }

    const PDParameter& operator[](const std::string &p_param_sid) const
    {
        auto res = m_params.find(p_param_sid);
        if (res != m_params.cend())
            return res->second;
        else throw std::out_of_range("Parameter with key '" + p_param_sid + "' isn't found");
    }
};


/// Тип перечня секций Parameter Data (PD).
/// Ключ - PD Sequence Number
using ParameterDataSections = std::map<int, stc_ParameterDataSection>;


#endif // UNT_IGES_DATA_TYPES