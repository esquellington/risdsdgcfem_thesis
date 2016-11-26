#ifndef UTIL_GSTATICIAN_H
#define UTIL_GSTATICIAN_H

#include <util/Config.h>
#include <Mal/RealUtils.h>
#include <vector>
#include <algorithm>

namespace util
{

/*! Simple statician that computes basic magnitudes for given samples of type T.
  \note Samples are stored in insertion order
  \todo Cache sorted sample vector for median and histogram
*/
template <typename T>
class GStatician
{
public:
    struct histogram_entry
    {
        histogram_entry( unsigned int count, T slice_begin, T slice_end )
        : m_Count(count), m_SliceBegin(slice_begin), m_SliceEnd(slice_end) {}
        unsigned int m_Count;
        T m_SliceBegin;
        T m_SliceEnd;
    };
    typedef std::vector< histogram_entry > histogram_type;

public:
    inline GStatician() { Clear(); }
    inline ~GStatician() {}

    inline void Clear()
        {
            m_vecSamples.clear();
            m_MinIdx = 0; m_MaxIdx = 0;
            m_Min = T(0); m_Max = T(0);
            m_Avg = 0; m_Stdev = 0;
        }
    inline void Begin( unsigned int num_reserved_samples = 0 )
        {
            Clear();
            if( num_reserved_samples > 0 ) m_vecSamples.reserve(num_reserved_samples);
        }
    inline void Add( const T &sample )
        {
            m_vecSamples.push_back(sample);
            if( 1 == m_vecSamples.size() )
            {
                m_MinIdx = 0;
                m_MaxIdx = 0;
            }
            else
            {
                if( sample < m_vecSamples[m_MinIdx] ) m_MinIdx = m_vecSamples.size()-1;
                if( sample > m_vecSamples[m_MaxIdx] ) m_MaxIdx = m_vecSamples.size()-1;
            }
        }
    inline void End()
        {
            if( 0 == m_vecSamples.size() ) return;
            // Min/Max values
            m_Min = m_vecSamples[m_MinIdx];
            m_Max = m_vecSamples[m_MaxIdx];
            // Average
            T acc_sum(0);
            for( unsigned i=0; i<m_vecSamples.size(); i++ )
                acc_sum += m_vecSamples[i];
            m_Avg = (m_vecSamples.size()>0) ? double(acc_sum)/m_vecSamples.size() : 0;
            // Variance and Stdev
            double acc_variance(0);
            for( unsigned i=0; i<m_vecSamples.size(); i++ )
                acc_variance += mal::Sq( double(m_vecSamples[i])-m_Avg );
            m_Stdev = (m_vecSamples.size()>1) ? mal::Sqrt( acc_variance / (m_vecSamples.size()-1) ) : 0;
        }
    inline unsigned int GetNumSamples() const { return m_vecSamples.size(); }
    inline T GetSample( int aIdx ) const { return m_vecSamples[aIdx]; }

    //\name Sorting functionality, consider caching sorted sample vector
    //@{
    inline T ComputeMedian()
        {
            if( 0 == m_vecSamples.size() ) return 0;
            // Median
            std::vector<T> vec_sorted_samples( m_vecSamples );
            std::sort( vec_sorted_samples.begin(), vec_sorted_samples.end() ); //increasing order
            return vec_sorted_samples[ vec_sorted_samples.size() / 2 ];
        }
    inline void ComputeHistogram( int num_slices, histogram_type &histogram ) const
        {
            if( 0 == m_vecSamples.size() ) return;

            std::vector<T> vec_sorted_samples( m_vecSamples );
            std::sort( vec_sorted_samples.begin(), vec_sorted_samples.end() ); //increasing order

            histogram.clear();
            UTIL_ASSERT( num_slices > 1 );

            T slice_size( (m_Max-m_Min) / (num_slices-1) );
            T slice_begin( m_Min );
            T slice_end( slice_begin + slice_size );
            if( slice_begin == slice_end ) //This may happen if all samples are equal or too many slices are requested
            {
                histogram.push_back( histogram_entry( vec_sorted_samples.size(), m_Min, m_Max ) );
                return;
            }
            histogram.push_back( histogram_entry( 0, slice_begin, slice_end ) );
            unsigned int it_sample(0);
            unsigned int it_slice(0);
            while( it_sample < vec_sorted_samples.size() )
            {
                if( vec_sorted_samples[it_sample] < slice_end )
                {
                    histogram[it_slice].m_Count++;
                    it_sample++;
                }
                else
                {
                    slice_begin = slice_end;
                    slice_end += slice_size;
                    histogram.push_back( histogram_entry( 0, slice_begin, slice_end ) );
                    it_slice++;
                }
            }
        }
    //@}

public:
    T m_Min;
    T m_Max;
    double m_Avg;
    double m_Stdev;

private:
    std::vector<T> m_vecSamples; //\note Samples are stored in insertion order
    int m_MinIdx; //\note Invalid after SORTING
    int m_MaxIdx; //\note Invalid after SORTING
};

} //namespace util

#endif //UTIL_GSTATICIAN_H
