#ifndef ABSTRACT_IGES_PARAMS_PARSERS
#define ABSTRACT_IGES_PARAMS_PARSERS

#include <vector>
#include <memory>
#include "IGESDataTypes.h"

/// Интерфейс классов-парсеров параметров секции Parameter Data (PD) IGES-файла
class AbstarctIGESParamsParsers
{
public:

    static std::shared_ptr<AbstarctIGESParamsParsers> Make(int p_Type);

    virtual ~AbstarctIGESParamsParsers() {}
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr) = 0;

protected:
    std::vector<std::string> Split(const std::string &p_Str, char p_Delim);
};

#endif // ABSTRACT_IGES_PARAMS_PARSERS
