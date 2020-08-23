#ifndef IFC_GEO_LOADER
#define IFC_GEO_LOADER

//#include <AbstractGeometryModel>
#include <string>
#include <memory>


/// \brief Базовый класс загрузчиков геометрии
class AbstractGeoLoader
{
public:
    enum GeoModels
    {
        parametric_geometry_model
    };

public:
    /// \brief Конструктор
    AbstractGeoLoader();

    /// \brief Деструктор
    virtual ~AbstractGeoLoader()
    { }

    /// \brief Метод загрузки геометрии
    virtual int loadGeo() = 0;

    /// \brief Метод создает модель геометрии по загруженным данным
    /// \param p_geo_model Тип модели геометрии
    /// \return Указатель на модель или пустой указатель, если формат не поддерживает данную модель
    
    //virtual std::shared_ptr<AbstractGeometryModel> getGeoModel(GeoModels p_geo_model) = 0;

    /// \brief Метод устанавливает имя файла геометрии
    void setGeoFileName(const std::string &p_geo_file_name)
    {
        m_geo_file_name = p_geo_file_name;
    }

    /// \brief Метод возвращает имя файла геометрии
    const std::string& getGeoFileName() const
    {
        return m_geo_file_name;
    }

private:
    std::string m_geo_file_name;
};


#endif // IFC_GEO_LOADER
