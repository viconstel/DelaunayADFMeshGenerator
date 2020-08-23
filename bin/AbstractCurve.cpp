#include "AbstractCurve.h"
#include <omp.h>

///////////////////////////////////////////////////////////////////////////////
// AbstractCurve
///////////////////////////////////////////////////////////////////////////////


AbstractCurve::AbstractCurve(
    unsigned p_dim,
    double   p_begin,
    double   p_end,
    bool     p_is_closed)
    : AbstractCurve(
          0,
          p_dim,
          p_begin,
          p_end,
          p_is_closed)
{ }


AbstractCurve::AbstractCurve(
    unsigned p_id,
    unsigned p_dim,
    double   p_begin,
    double   p_end,
    bool     p_is_closed)
    : m_id(p_id)
    , m_dim(p_dim)
    , m_begin(p_begin)
    , m_end(p_end)
    , m_is_closed(p_is_closed)
{
    initializeDependentData();
}

/*AbstractCurve::AbstractCurve(const geo::PrimitiveParameters &p_parameters, unsigned &p_n)
{
    m_id        = std::any_cast<unsigned>(p_parameters[p_n++]);
    m_dim       = std::any_cast<unsigned>(p_parameters[p_n++]);
    m_begin     = std::any_cast<double  >(p_parameters[p_n++]);
    m_end       = std::any_cast<double  >(p_parameters[p_n++]);
    m_is_closed = std::any_cast<bool    >(p_parameters[p_n++]);

    initializeDependentData();
}*/


void AbstractCurve::initializeDependentData()
{
    int threads_count = omp_get_max_threads();
    m_max_der_order = 2;

    m_counts.resize(threads_count * (m_max_der_order + 1));
    m_counts_array = m_counts.data();
}


AbstractCurve::~AbstractCurve()
{ }


double AbstractCurve::getBegin() const
{
    return m_begin;
}


double AbstractCurve::getEnd() const
{
    return m_end;
}


bool AbstractCurve::isClosed() const
{
    return m_is_closed;
}


unsigned AbstractCurve::getDim() const
{
    return m_dim;
}

unsigned AbstractCurve::getId() const
{
    return m_id;
}
