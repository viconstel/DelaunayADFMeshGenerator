#include "AbstarctIGESParamsParsers.h"
#include <sstream>



// Конкретные классы парсеров


/// Circular Arc Entity (Type 100)
class cls_Type100
    : public AbstarctIGESParamsParsers
{
public:
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr)
    {
        std::map<std::string, PDParameter> l_Params;
        std::vector<std::string> l_StrParamsValues = Split(p_ParamsStr, ',');

        l_Params["ZT"] = { IGESDataTypes::REAL, l_StrParamsValues[0] };
        l_Params["X1"] = { IGESDataTypes::REAL, l_StrParamsValues[1] };
        l_Params["Y1"] = { IGESDataTypes::REAL, l_StrParamsValues[2] };
        l_Params["X2"] = { IGESDataTypes::REAL, l_StrParamsValues[3] };
        l_Params["Y2"] = { IGESDataTypes::REAL, l_StrParamsValues[4] };
        l_Params["X3"] = { IGESDataTypes::REAL, l_StrParamsValues[5] };
        l_Params["Y3"] = { IGESDataTypes::REAL, l_StrParamsValues[6] };

        return l_Params;
    }
};


/// Rational B-Spline Curve Entity (Type 126)
class cls_Type126
    : public AbstarctIGESParamsParsers
{
public:
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr)
    {
        std::map<std::string, PDParameter> l_Params;
        std::vector<std::string> l_StrParamsValues = Split(p_ParamsStr, ',');

        l_Params["K"]     = { IGESDataTypes::INTEGER, l_StrParamsValues[0] };
        l_Params["M"]     = { IGESDataTypes::INTEGER, l_StrParamsValues[1] };
        l_Params["PROP1"] = { IGESDataTypes::INTEGER, l_StrParamsValues[2] };
        l_Params["PROP2"] = { IGESDataTypes::INTEGER, l_StrParamsValues[3] };
        l_Params["PROP3"] = { IGESDataTypes::INTEGER, l_StrParamsValues[4] };
        l_Params["PROP4"] = { IGESDataTypes::INTEGER, l_StrParamsValues[5] };

        int l_K = l_Params["K"].toInteger();
        int l_M = l_Params["M"].toInteger();
        int l_N = 1 + l_K - l_M;
        int l_A = l_N + 2 * l_M;

        for (int i = 0; i < l_A + 1; ++i)
            l_Params["T(" + std::to_string(-l_M + i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[6 + i] };
        
        for (int i = 0; i < l_K + 1; ++i)
            l_Params["W(" + std::to_string(i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[7 + l_A + i] };

        for (int p = 0, i = 0; p < 3 * l_K + 1; p += 3, ++i)
        {
            l_Params["X(" + std::to_string(i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[8  + l_A + l_K + p] };
            l_Params["Y(" + std::to_string(i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[9  + l_A + l_K + p] };
            l_Params["Z(" + std::to_string(i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[10 + l_A + l_K + p] };
        }

        l_Params["V(0)"]  = { IGESDataTypes::REAL, l_StrParamsValues[11 + l_A + 4 * l_K] };
        l_Params["V(1)"]  = { IGESDataTypes::REAL, l_StrParamsValues[12 + l_A + 4 * l_K] };
        l_Params["XNORM"] = { IGESDataTypes::REAL, l_StrParamsValues[13 + l_A + 4 * l_K] };
        l_Params["YNORM"] = { IGESDataTypes::REAL, l_StrParamsValues[14 + l_A + 4 * l_K] };
        l_Params["ZNORM"] = { IGESDataTypes::REAL, l_StrParamsValues[15 + l_A + 4 * l_K] };
        
        return l_Params;
    }
};


/// Rational B-Spline Surface Entity (Type 128)
class cls_Type128
    : public AbstarctIGESParamsParsers
{
public:
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr)
    {
        std::map<std::string, PDParameter> l_Params;
        std::vector<std::string> l_StrParamsValues = Split(p_ParamsStr, ',');

        l_Params["K1"] = { IGESDataTypes::INTEGER, l_StrParamsValues[0] };
        l_Params["K2"] = { IGESDataTypes::INTEGER, l_StrParamsValues[1] };
        l_Params["M1"] = { IGESDataTypes::INTEGER, l_StrParamsValues[2] };
        l_Params["M2"] = { IGESDataTypes::INTEGER, l_StrParamsValues[3] };
        l_Params["PROP1"] = { IGESDataTypes::INTEGER, l_StrParamsValues[4] };
        l_Params["PROP2"] = { IGESDataTypes::INTEGER, l_StrParamsValues[5] };
        l_Params["PROP3"] = { IGESDataTypes::INTEGER, l_StrParamsValues[6] };
        l_Params["PROP4"] = { IGESDataTypes::INTEGER, l_StrParamsValues[7] };
        l_Params["PROP5"] = { IGESDataTypes::INTEGER, l_StrParamsValues[8] };

        int l_K1 = l_Params["K1"].toInteger();
        int l_K2 = l_Params["K2"].toInteger();
        int l_M1 = l_Params["M1"].toInteger();
        int l_M2 = l_Params["M2"].toInteger();
        int l_N1 = 1 + l_K1 - l_M1;
        int l_N2 = 1 + l_K2 - l_M2;
        int l_A  = l_N1 + 2 * l_M1;
        int l_B  = l_N2 + 2 * l_M2;
        int l_C  = (1 + l_K1) * (1 + l_K2);

        for (int i = 0; i < l_A + 1; ++i)
            l_Params["S(" + std::to_string(-l_M1 + i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[9 + i] };

        for (int i = 0; i < l_B + 1; ++i)
            l_Params["T(" + std::to_string(-l_M2 + i) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[10 + l_A + i] };
        
        for (int j = 0; j < l_K2 + 1; ++j)
            for (int i = 0; i < l_K1 + 1; ++i)
                l_Params["W(" + std::to_string(i) + "," + std::to_string(j) + ")"] = { IGESDataTypes::REAL, l_StrParamsValues[11 + l_A + l_B + i + j * (l_K1 + 1)] };

        for (int p2 = 0, j = 0; p2 < 3 * l_K2 + 1; p2 += 3, ++j)
            for (int p1 = 0, i = 0; p1 < 3 * l_K1 + 1; p1 += 3, ++i)
            {
                l_Params["X(" + std::to_string(i) + "," + std::to_string(j) + ")"] =
                    { IGESDataTypes::REAL, l_StrParamsValues[11 + l_A + l_B + l_C + p1 + p2 * (l_K1 + 1)] };
                l_Params["Y(" + std::to_string(i) + "," + std::to_string(j) + ")"] =
                    { IGESDataTypes::REAL, l_StrParamsValues[12 + l_A + l_B + l_C + p1 + p2 * (l_K1 + 1)] };
                l_Params["Z(" + std::to_string(i) + "," + std::to_string(j) + ")"] =
                    { IGESDataTypes::REAL, l_StrParamsValues[13 + l_A + l_B + l_C + p1 + p2 * (l_K1 + 1)] };
            }

        l_Params["U(0)"] = { IGESDataTypes::REAL, l_StrParamsValues[11 + l_A + l_B + 4 * l_C] };
        l_Params["U(1)"] = { IGESDataTypes::REAL, l_StrParamsValues[12 + l_A + l_B + 4 * l_C] };
        l_Params["V(0)"] = { IGESDataTypes::REAL, l_StrParamsValues[13 + l_A + l_B + 4 * l_C] };
        l_Params["V(1)"] = { IGESDataTypes::REAL, l_StrParamsValues[14 + l_A + l_B + 4 * l_C] };

        return l_Params;
    }
};


/// Boundary Entity (Type 141)
class cls_Type141
    : public AbstarctIGESParamsParsers
{
public:
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr)
    {
        std::map<std::string, PDParameter> l_Params;
        std::vector<std::string> l_StrParamsValues = Split(p_ParamsStr, ',');
        l_Params["TYPE"]     = { IGESDataTypes::INTEGER, l_StrParamsValues[0] };
        l_Params["PREF"]     = { IGESDataTypes::INTEGER, l_StrParamsValues[1] };
        l_Params["SPTR"]     = { IGESDataTypes::POINTER, l_StrParamsValues[2] };
        l_Params["N"]        = { IGESDataTypes::INTEGER, l_StrParamsValues[3] };

        unsigned l_N = l_Params["N"].toInteger();
        unsigned l_K, l_Sum{ 0 };
        for (int i = 0; i < l_N; ++i)
        {
            l_Params["CRVPT(" + std::to_string(i + 1) + ")"] = { IGESDataTypes::POINTER, l_StrParamsValues[4 + 3 * i + l_Sum] };
            l_Params["SENSE(" + std::to_string(i + 1) + ")"] = { IGESDataTypes::INTEGER, l_StrParamsValues[5 + 3 * i + l_Sum] };
            l_Params["K(" + std::to_string(i + 1) + ")"]     = { IGESDataTypes::INTEGER, l_StrParamsValues[6 + 3 * i + l_Sum] };

            l_K = l_Params["K(" + std::to_string(i + 1) + ")"].toInteger();
            for (int j = 0; j < l_K; ++j)
                l_Params["PSCPT(" + std::to_string(i + 1) + "," + std::to_string(j + 1) + ")"] = { IGESDataTypes::POINTER, l_StrParamsValues[7 + 3 * i + l_Sum + j] };
            l_Sum += l_K;
        }
        
        return l_Params;
    }
};


/// Bounded Surface Entity (Type 143)
class cls_Type143
    : public AbstarctIGESParamsParsers
{
public:
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr)
    {
        std::map<std::string, PDParameter> l_Params;
        std::vector<std::string> l_StrParamsValues = Split(p_ParamsStr, ',');
        l_Params["TYPE"] = { IGESDataTypes::INTEGER, l_StrParamsValues[0] };
        l_Params["SPTR"] = { IGESDataTypes::POINTER, l_StrParamsValues[1] };
        l_Params["N"]    = { IGESDataTypes::INTEGER, l_StrParamsValues[2] };
        unsigned l_N     = l_Params["N"].toInteger();
        for (int i = 0; i < l_N; ++i)
            l_Params["BDPT(" + std::to_string(i + 1) + ")"] = { IGESDataTypes::POINTER, l_StrParamsValues[3 + i] };
        return l_Params;
    }
};


/// Color Definition Entity (Type 314)
class cls_Type314
    : public AbstarctIGESParamsParsers
{
public:
    virtual std::map<std::string, PDParameter> Parse(const std::string &p_ParamsStr)
    {
        std::map<std::string, PDParameter> l_Params;
        std::vector<std::string> l_StrParamsValues = Split(p_ParamsStr, ',');
        l_Params["CC1"] = { IGESDataTypes::REAL,   l_StrParamsValues[0] };
        l_Params["CC2"] = { IGESDataTypes::REAL,   l_StrParamsValues[1] };
        l_Params["CC3"] = { IGESDataTypes::REAL,   l_StrParamsValues[2] };
        l_Params["CNAME"] = { IGESDataTypes::STRING, l_StrParamsValues[3] };
        return l_Params;
    }
};



// Реализация методов класса AbstarctIGESParamsParsers



std::shared_ptr<AbstarctIGESParamsParsers> AbstarctIGESParamsParsers::Make(int p_Type)
{
    switch (p_Type)
    {
    case 100:
        return std::make_shared<cls_Type100>();
    case 126:
        return std::make_shared<cls_Type126>();
    case 128:
        return std::make_shared<cls_Type128>();
    case 141:
        return std::make_shared<cls_Type141>();
    case 143:
        return std::make_shared<cls_Type143>();
    case 314:
        return std::make_shared<cls_Type314>();
    default:
        return nullptr;
    }
}


std::vector<std::string> AbstarctIGESParamsParsers::Split(const std::string &p_Str, char p_Delim)
{
    std::vector<std::string> l_SubStrings;
    std::istringstream l_In(p_Str);
    std::string l_Temp;
    while (l_In)
    {
        std::getline(l_In, l_Temp, p_Delim);
        l_SubStrings.push_back(l_Temp);
    }
    return l_SubStrings;
}
