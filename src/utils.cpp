#include<utils.h>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

double combine_pval(vector<double> p_vals)
{
    // Combining pvalues using fisher's method.
    // Reference: https://en.wikipedia.org/wiki/Fisher%27s_method
    int num_pvals = p_vals.size();
    double test_statistic = 0;
    for (double p_val: p_vals){
        test_statistic -= 2 * log(p_val);
    }
    int dof = 2 * num_pvals; // degrees of freedom for chi-squared distribution

    boost::math::chi_squared dist(dof);
    double pValue = cdf(complement(dist, test_statistic));
    
    return pValue;
}

