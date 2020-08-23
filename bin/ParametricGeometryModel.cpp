#include "ParametricGeometryModel.h"


///////////////////////////////////////////////////////////////////////////////
// SurfaceBoundary
///////////////////////////////////////////////////////////////////////////////


SurfaceBoundary::SurfaceBoundary(
    unsigned p_curve_count)
{
    m_surface_boundary.reserve(p_curve_count);
}


void
SurfaceBoundary::addCurve(
    unsigned                       p_curve_id,
    bool                           p_is_direction_change,
    std::shared_ptr<AbstractCurve> p_curve)
{
    return m_surface_boundary.push_back({
        p_curve_id,
        { p_curve, p_is_direction_change }
    });
}


SurfaceBoundary::OrientedCurve
SurfaceBoundary::getCurve(
    unsigned p_curve_id)
{
    auto curve_iter = std::find_if(
        m_surface_boundary.begin(),
        m_surface_boundary.end(),
        [&](const auto &p_data){
            return p_data.first == p_curve_id;
        });

    if (curve_iter != m_surface_boundary.end())
        return curve_iter->second;

    throw EntityNotFoundException("The curve with ID = " +
                                  std::to_string(p_curve_id) +
                                  " was not found");
}


const SurfaceBoundary::SurfaceBoundaryContType&
SurfaceBoundary::getCurves() const
{
    return m_surface_boundary;
}


SurfaceBoundary::SurfaceBoundaryContType&
SurfaceBoundary::getCurves()
{
    return m_surface_boundary;
}


///////////////////////////////////////////////////////////////////////////////
// BoundedSurface
///////////////////////////////////////////////////////////////////////////////


BoundedSurface::BoundedSurface(
    std::shared_ptr<AbstractSurface> p_surface,
    unsigned                         p_boundary_count)
    : m_surface(p_surface)
{
    m_surface_boundaries.reserve(p_boundary_count);
}


void
BoundedSurface::setSurface(
    std::shared_ptr<AbstractSurface> p_surface)
{
    m_surface = p_surface;
}


void BoundedSurface::reserve(
    unsigned p_boundary_count)
{
    m_surface_boundaries.reserve(p_boundary_count);
}

/// \brief Получение NURBS-поверхности
const std::shared_ptr<AbstractSurface>
BoundedSurface::getSurface() const
{
    return m_surface;
}

/// \brief Получение NURBS-поверхности
std::shared_ptr<AbstractSurface>
BoundedSurface::getSurface()
{
    return m_surface;
}

/// \brief Добавление границы в конец контейнера
void
BoundedSurface::addBoundary(
    unsigned               p_boundary_id,
    const SurfaceBoundary &p_boundary)
{
    return m_surface_boundaries.push_back({ p_boundary_id, p_boundary });
}


const SurfaceBoundary&
BoundedSurface::getBoundary(
    unsigned p_boundary_id)
{
    BoundedSurfaceContType::iterator boundary_iter = std::find_if(
        m_surface_boundaries.begin(),
        m_surface_boundaries.end(),
        [&](const std::pair<unsigned, SurfaceBoundary> &p_data){
            return p_data.first == p_boundary_id;
        });

    if (boundary_iter != m_surface_boundaries.end())
        return boundary_iter->second;

    throw EntityNotFoundException("The boundary with the ID = " +
                                  std::to_string(p_boundary_id) +
                                  " was not found");
}


const BoundedSurface::BoundedSurfaceContType&
BoundedSurface::getBoundaries() const
{
    return m_surface_boundaries;
}


BoundedSurface::BoundedSurfaceContType&
BoundedSurface::getBoundaries()
{
    return m_surface_boundaries;
}


///////////////////////////////////////////////////////////////////////////////
// Body
///////////////////////////////////////////////////////////////////////////////


bool
Body::addBoundedSurface(
    unsigned                        p_bounded_surface_id,
    std::shared_ptr<BoundedSurface> p_bounded_surface,
    bool                            p_orientation_flag)
{
    return m_body.insert({
        p_bounded_surface_id,
        { p_bounded_surface, p_orientation_flag }
    }).second;
}


const Body::OrientedBoundedSurface&
Body::getBoundedSurface(
    unsigned p_bounded_surface_id) const
{
    auto bounded_surface_iter = m_body.find(p_bounded_surface_id);
    if (bounded_surface_iter != m_body.end())
        return bounded_surface_iter->second;

    throw EntityNotFoundException("The bounded surface with the ID = " +
                                  std::to_string(p_bounded_surface_id) +
                                  " was not found");
}


const Body::BodyContType&
Body::getBoundedSurfaces() const
{
    return m_body;
}


Body::BodyContType&
Body::getBoundedSurfaces()
{
    return m_body;
}


///////////////////////////////////////////////////////////////////////////////
// ParametricGeometryModel
///////////////////////////////////////////////////////////////////////////////

/*
bool
ParametricGeometryModel::addBody(
    unsigned    p_body_id,
    const Body &p_body)
{
    return m_bodies.emplace(p_body_id, p_body).second;
}


const Body&
ParametricGeometryModel::getBody(
    unsigned p_body_id)
{
    auto body_iter = m_bodies.find(p_body_id);
    if (body_iter != m_bodies.end())
        return body_iter->second;

    throw EntityNotFoundException("The body with the ID = " +
                                  std::to_string(p_body_id) +
                                  " was not found");
}


const ParametricGeometryModel::ParametricGeometryModelContType&
ParametricGeometryModel::getBodies() const
{
    return m_bodies;
}


ParametricGeometryModel::ParametricGeometryModelContType&
ParametricGeometryModel::getBodies()
{
    return m_bodies;
}


ParametricGeometryModel::ModelCurvesContainer
ParametricGeometryModel::getAllModelCurves()
{
    ModelCurvesContainer all_model_curves;
    for (const auto &[id, body] : m_bodies)
    {
        const Body::BodyContType &bounded_surfaces = body.getBoundedSurfaces();
        for (const auto &[bounded_surface_id, bounded_surface] : bounded_surfaces)
        {
            const BoundedSurface::BoundedSurfaceContType &boundaries =
                bounded_surface.first->getBoundaries();
            for (const auto &[boundary_id, boundary] : boundaries)
            {
                const SurfaceBoundary::SurfaceBoundaryContType &curves =
                    boundary.getCurves();
                for (const auto &[curve_id, curve] : curves)
                    all_model_curves.insert({ curve_id, curve.first });
            }
        }
    }
    return all_model_curves;
}


ParametricGeometryModel::ModelSurfacesContainer
ParametricGeometryModel::getAllModelBoundedSurfaces()
{
    ModelSurfacesContainer all_model_surfaces;
    for (const auto &[body_id, body] : m_bodies)
    {
        const Body::BodyContType &bounded_surfaces = body.getBoundedSurfaces();
        for (const auto &[bounded_surface_id, oriented_bounded_surface] : bounded_surfaces)
            all_model_surfaces.insert({
                bounded_surface_id,
                oriented_bounded_surface.first
            });
    }
    return all_model_surfaces;
}
*/

