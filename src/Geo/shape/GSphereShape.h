#ifndef GEO_SHAPE_GSPHERESHAPE_H
#define GEO_SHAPE_GSPHERESHAPE_H

#include <Geo/shape/IShape.h>

// BV types for ComputeBV()
#include <Geo/bv/bv.h>

#include <Mal/GRandom.h> //for GSphereShapeDomainSamplerD
#include <util/GIndexedRefCountedPoolDA.h> //for GSphereShapeDomainSamplerD

namespace geo {

//TEMP
template <unsigned D> class GSphereShapeDomainSamplerD;
class GSphereShapeDomainSamplerD<2>;
class GSphereShapeDomainSamplerD<3>;

/*! 2/3d sphere shape. */
template <unsigned D>
class GSphereShape: public GIShapeD<D>
{
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;
    typedef Real sdof_type; //!< Useless but cannot be void
    //@}

    typedef GSphereShape<D> self_type;

public:
    GSphereShape( Real radius = Real(1) ) : m_Radius(radius) {}
    ~GSphereShape() {}

    //!\name Shape params
    //@{
    finline void SetRadius( Real radius ) { m_Radius = radius; }
    finline Real GetRadius() const { return m_Radius; }
    //@}

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return shape_type_of<self_type>(); }
    unsigned int GetNumSDOF() const { return 0; }
    const sdof_type *GetVecDefaultSDOF() const { return 0; }
    const Real *GetVecDefaultDOF() const { return 0; }
    GIDomainSamplerD<D> *CreateDomainSampler() const; //\note Specialized below
    //@}

    //dimension-specific
    void ComputeBVD( bv::GBoundingVolumeD<D> &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
        {
            switch( bv.GetType() )
            {
            case bv::eBV_Sphere2:
            case bv::eBV_Sphere3:
                bv.template As< bv::GSphere<D> >().SetPosRadius( transform.m_Pos, m_Radius );
                break;
            case bv::eBV_AABB2:
            case bv::eBV_AABB3:
                bv.template As< bv::GAABB<D> >().SetPosHalfSizes( transform.m_Pos, mal::GVec<Real,D>(m_Radius) );
                break;
            case bv::eBV_LSS2:
            case bv::eBV_LSS3:
                bv.template As< bv::GLineSweptSphere<D> >().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
                break;
            case bv::eBV_Void:
                break;
            case bv::eBV_Infinite:
                break;
            default:
                GEO_ASSERT( false ); //wrong type or dimension
                break;
            }
        }

protected:
    Real m_Radius;
};

//---- 2d Specialization of GSphereShapeDomainSamplerD<D>
template<>
class GSphereShapeDomainSamplerD<2>: public GIDomainSamplerD<2>
{
public:
    // !Stores minimal local data to identify a POB and compute its geometry
    struct PointOnBoundary
    {
        PointOnBoundary( Real angle_rad = Real(0) ) : m_AngleRad(angle_rad) {}
        ~PointOnBoundary() {}
        Real m_AngleRad;
    };

public:
    GSphereShapeDomainSamplerD<2>( const GSphereShape<2> *p_shape ) : m_pShape(p_shape) {}
    ~GSphereShapeDomainSamplerD<2>() {}

    unsigned int GetNumPOB() const { return m_poolPOB.Size(); }
    unsigned int GetNumAllocatedPOB() const { return m_poolPOB.NumAllocated(); }

    void UpdateDOF( const Real *vec_dof ) {}
    void ClearPOB() { m_poolPOB.Clear(); }

    //\name PointOnBoundary creation
    //@{
    point_on_boundary_id_type CreatePOB( const PointOnFeature &pof )
        {
            GEO_LOG_ASSERT(false, "GSphereShapeDomainSamplerD<2>::CreatePOB undefined!");
            return cInvalidPOB;
        }
    point_on_boundary_id_type StepPOB( point_on_boundary_id_type id,
                                       const vec_type &dir_local, Real step_length, Real eps_length )
        {
            Real angle_rad_step( step_length / m_pShape->GetRadius() );
            vec_type normal( mal::Cos( m_poolPOB[id].m_AngleRad ), mal::Sin( m_poolPOB[id].m_AngleRad ) );
            if( mal::Dot(normal,mal::PerpendicularCW(dir_local)) > 0 ) angle_rad_step = -angle_rad_step; // Select rotation sign
            return m_poolPOB.New( PointOnBoundary( m_poolPOB[id].m_AngleRad + angle_rad_step ) );
        }
    void CreateNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                             point_on_boundary_id_type id,
                             Real neighbourhood_radius, Real neighbourhood_epsilon )
        {
            GEO_ASSERT( count == 2 );
            Real angle_rad_step( neighbourhood_radius / m_pShape->GetRadius() );
            vec_id[0] = m_poolPOB.New( PointOnBoundary( m_poolPOB[id].m_AngleRad - angle_rad_step ) );
            vec_id[1] = m_poolPOB.New( PointOnBoundary( m_poolPOB[id].m_AngleRad + angle_rad_step ) );
        }
    point_on_boundary_id_type ClosestPOB( point_on_boundary_id_type id,
                                          const vec_type &pos_local,
                                          Real neighbourhood_radius, Real neighbourhood_epsilon )
        {
            GEO_ASSERT(false);
            return id;
            /*\todo
            Real neighbourhood_radius_angle_rad( neighbourhood_radius / m_pShape->GetRadius() );
            Real angle_rad( mal::Atan2( pos_local[1], pos_local[0] ) );
            if( mal::Abs(angle_rad - m_poolPOB[id].m_AngleRad) < neighbourhood_epsilon )
                return id;
            else
            {
                //\todo choose CW/CCW for min angle first, clamp afterwards!
                angle_rad = mal::Clamp( angle_rad,
                                        m_poolPOB[id].m_AngleRad - neighbourhood_radius_angle_rad,
                                        m_poolPOB[id].m_AngleRad + neighbourhood_radius_angle_rad );
                return m_poolPOB.New( PointOnBoundary( angle_rad ) );
            }
            */
        }
    point_on_boundary_id_type CreateRandomPOB()
        { return m_poolPOB.New( PointOnBoundary( mal::RandomF<Real>( 0, mal::TwoPi<Real>() ) ) ); }
    void CreateRandomPOB( point_on_boundary_id_type *vec_id, unsigned int count )
        { for( unsigned int i=0; i<count; i++ ) vec_id[i] = m_poolPOB.New( PointOnBoundary( mal::RandomF<Real>( 0, mal::TwoPi<Real>() ) ) ); }
    void CreateRandomNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                                   Real neighbourhood_radius, point_on_boundary_id_type id )
        {
            //\todo Consider normalizing resulting angle to [0,2Pi)
            Real angle_rad_range( neighbourhood_radius / m_pShape->GetRadius() );
            for( unsigned int i=0; i<count; i++ )
                vec_id[i] = m_poolPOB.New( PointOnBoundary( m_poolPOB[id].m_AngleRad + mal::RandomF<Real>( -angle_rad_range, angle_rad_range ) ) );
        }
    //@}

    //\name array-POB lifetime methods
    //@{
    void APOB_Destroy( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.Delete( vec_id[i] ); }
    void APOB_IncRef( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.IncRef( vec_id[i] ); }
    void APOB_DecRef( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.DecRef( vec_id[i] ); }
    //@}

    void AllPositionsAndNormals( vec_type *vec_pos, vec_type *vec_normal ) const
        {
            for( unsigned int it_pob=0; it_pob<m_poolPOB.NumAllocated(); it_pob++ ) //\note We iterate over ALL possible POB, but update only valid ones
                if( m_poolPOB.IsValid(it_pob) )
                {
                    vec_normal[it_pob] = Vec2( mal::Cos( m_poolPOB[it_pob].m_AngleRad ), mal::Sin( m_poolPOB[it_pob].m_AngleRad ) );
                    vec_pos[it_pob] = m_pShape->GetRadius() * vec_normal[it_pob];
                }
        }
    Vec2 POB_Position( point_on_boundary_id_type id ) const
        { return m_pShape->GetRadius() * POB_Normal(id); }
    Vec2 POB_Normal( point_on_boundary_id_type id ) const
        { return Vec2( mal::Cos( m_poolPOB[id].m_AngleRad ), mal::Sin( m_poolPOB[id].m_AngleRad ) ); }
    feature_id POB_FeatureId( point_on_boundary_id_type id ) const
        { return feature_id(); } //\note No features...
    PointOnFeature POB_PointOnFeature( point_on_boundary_id_type id ) const
        {
            GEO_LOG_ASSERT(false, "GSphereShapeDomainSamplerD<2>::POB_PointOnFeature undefined!");
            return PointOnFeature( feature_id(), Vec4::Zero() );
        }

    const GIShapeD<2> *GetShape() const { return m_pShape; }

    //TEMP
    void Trace() const { m_poolPOB.Trace(); }

private:
    const GSphereShape<2> *m_pShape;
    util::GIndexedRefCountedPoolDA<PointOnBoundary,point_on_boundary_id_type> m_poolPOB;
};

//---- 3d Specialization of GSphereShapeDomainSamplerD<D>
template<>
class GSphereShapeDomainSamplerD<3>: public GIDomainSamplerD<3>
{
public:
    // !Stores minimal local data to identify a POB and compute its geometry
    struct PointOnBoundary
    {
        PointOnBoundary( Real alpha_rad = Real(0), Real beta_rad = Real(0) ) : m_AlphaRad(alpha_rad), m_BetaRad(beta_rad) {}
        ~PointOnBoundary() {}
        Real m_AlphaRad;
        Real m_BetaRad;
    };

public:
    GSphereShapeDomainSamplerD<3>( const GSphereShape<3> *p_shape ) : m_pShape(p_shape) {}
    ~GSphereShapeDomainSamplerD<3>() {}

    unsigned int GetNumPOB() const { return m_poolPOB.Size(); }
    unsigned int GetNumAllocatedPOB() const { return m_poolPOB.NumAllocated(); }

    void UpdateDOF( const Real *vec_dof ) {}
    void ClearPOB() { m_poolPOB.Clear(); }

    //\name PointOnBoundary creation
    //@{
    point_on_boundary_id_type CreatePOB( const PointOnFeature &pof )
        {
            GEO_LOG_ASSERT(false, "GSphereShapeDomainSamplerD<2>::CreatePOB undefined!");
            return cInvalidPOB;
        }
    point_on_boundary_id_type StepPOB( point_on_boundary_id_type id,
                                       const vec_type &dir_local, Real step_length, Real eps_length )
        { GEO_ASSERT(false); return point_on_boundary_id_type(cInvalidPOB); }

    void CreateNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                             point_on_boundary_id_type id,
                             Real radius, Real eps_radius ) { GEO_ASSERT(false); }
    point_on_boundary_id_type ClosestPOB( point_on_boundary_id_type id,
                                          const vec_type &pos_local,
                                          Real neighbourhood_radius, Real neighbourhood_epsilon )
        { GEO_ASSERT(false); return id; }

    point_on_boundary_id_type CreateRandomPOB()
        { return m_poolPOB.New( PointOnBoundary( mal::RandomF<Real>( 0, mal::Pi<Real>() ), mal::RandomF<Real>( 0, mal::TwoPi<Real>() ) ) ); }
    void CreateRandomPOB( point_on_boundary_id_type *vec_id, unsigned int count )
        { GEO_ASSERT(false); }//\todo for( unsigned int i=0; i<count; i++ ) vec_id[i] = m_poolPOB.New( PointOnBoundary( mal::Random( Real(0), mal::TwoPi<Real>() ) ) ); }
    void CreateRandomNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                                   Real radius, point_on_boundary_id_type id )
        {
            GEO_ASSERT(false);
            /*\todo
            Real angle_rad_range( radius / m_pShape->GetRadius() );
            for( unsigned int i=0; i<count; i++ )
                vec_id[i] = m_poolPOB.New( PointOnBoundary( m_poolPOB[id].m_AngleRad + mal::Random( -angle_rad_range, angle_rad_range ) ) );
            */
        }
    //@}

    //\name array-POB lifetime methods
    //@{
    void APOB_Destroy( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.Delete( vec_id[i] ); }
    void APOB_IncRef( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.IncRef( vec_id[i] ); }
    void APOB_DecRef( const point_on_boundary_id_type *vec_id, unsigned int count ) { for(unsigned int i=0;i<count;i++) m_poolPOB.DecRef( vec_id[i] ); }
    //@}

    void AllPositionsAndNormals( vec_type *vec_pos, vec_type *vec_normal ) const
        {
            for( unsigned int it_pob=0; it_pob<m_poolPOB.NumAllocated(); it_pob++ ) //\note We iterate over ALL possible POB, but update only valid ones
                if( m_poolPOB.IsValid(it_pob) )
                {
                    vec_normal[it_pob] = Vec3( mal::Sin( m_poolPOB[it_pob].m_AlphaRad * mal::Cos( m_poolPOB[it_pob].m_BetaRad ) ),
                                               mal::Sin( m_poolPOB[it_pob].m_AlphaRad * mal::Sin( m_poolPOB[it_pob].m_BetaRad ) ),
                                               mal::Cos( m_poolPOB[it_pob].m_AlphaRad ) );
                    vec_pos[it_pob] = m_pShape->GetRadius() * vec_normal[it_pob];
                }
        }

    Vec3 POB_Position( point_on_boundary_id_type id ) const { return m_pShape->GetRadius() * POB_Normal( id ); }
    Vec3 POB_Normal( point_on_boundary_id_type id ) const
        {
            return Vec3( mal::Sin( m_poolPOB[id].m_AlphaRad * mal::Cos( m_poolPOB[id].m_BetaRad ) ),
                         mal::Sin( m_poolPOB[id].m_AlphaRad * mal::Sin( m_poolPOB[id].m_BetaRad ) ),
                         mal::Cos( m_poolPOB[id].m_AlphaRad ) );
        }
    feature_id POB_FeatureId( point_on_boundary_id_type id ) const
        { return feature_id(); } //\note No features...
    PointOnFeature POB_PointOnFeature( point_on_boundary_id_type id ) const
        {
            GEO_LOG_ASSERT(false, "GSphereShapeDomainSamplerD<2>::POB_PointOnFeature undefined!");
            return PointOnFeature( feature_id(), Vec4::Zero() );
        }

    const GIShapeD<3> *GetShape() const { return m_pShape; }

    //TEMP
    void Trace() const { m_poolPOB.Trace(); }

private:
    const GSphereShape<3> *m_pShape;
    util::GIndexedRefCountedPoolDA<PointOnBoundary,point_on_boundary_id_type> m_poolPOB;
};

//---- Specialization of CreateDomainSampler() to use unrelated dimension-specific implementations
template<>
inline GIDomainSamplerD<2> *GSphereShape<2>::CreateDomainSampler() const
{
    return new GSphereShapeDomainSamplerD<2>(this);
}

template<>
inline GIDomainSamplerD<3> *GSphereShape<3>::CreateDomainSampler() const
{
    return new GSphereShapeDomainSamplerD<3>(this);
}

} //namespace geo

#endif // GEO_SHAPE_GSPHERESHAPE_H
