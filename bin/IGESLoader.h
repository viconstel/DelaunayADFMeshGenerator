#ifndef CLS_IGES_LOADER
#define CLS_IGES_LOADER

#include "AbstractGeoLoader.h"
#include "IGESData.h"
#include "IGESDataTypes.h"
#include <map>


/// \brief Загрузчик файла IGES
class IGESLoader
    : public AbstractGeoLoader
{
public:
    /// \brief Метод загрузки геометрии
    int loadGeo() override;

    /// \brief Метод создает модель геометрии по загруженным данным
    /// \param p_geo_model Тип модели геометрии
    /// \return Указатель на модель или пустой указатель, если формат не поддерживает данную модель
   
    //std::shared_ptr<AbstractGeometryModel> getGeoModel(GeoModels p_geo_model) override;

    /// \brief Метод возвращает все загруженные данные из файла геометрии
    IGESData& GetIGESModel();

private:
    IGESData m_iges_data;
};

#endif // CLS_IGES_LOADER