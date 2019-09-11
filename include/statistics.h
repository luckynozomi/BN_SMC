#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdio.h>
#include <math.h>
#include <vector>

using namespace std;
class Statistics
{
    public:
        Statistics();
        virtual ~Statistics();
        template < typename T, typename A >
        void Summary(const vector< T , A >& val)
        {
        	_mean = 0.0;
			_variance = 0.0;
			_std = 0.0;

		    // Compute the mean first
			for(int i =0; i<val.size();i++)
			{
				_mean+=val[i];
			}
			_mean/= val.size();
		    
            // Then variance
			for(int i = 0; i < val.size(); i++)
			{
				_variance += pow((val[i]-_mean),2.0);
			}
			_variance /= (val.size()-1);
		    
            // The standard deviation
			_std = sqrt(_variance);
        }

        double Get_mean()
        {
            return(_mean);
        }
        double Get_variance()
        {
            return(_variance);
        }
        double Get_std()
        {
            return(_std);
        }
    private:
        double _mean;
        double _variance;
        double _std;
};

#endif // STATISTICS_H
