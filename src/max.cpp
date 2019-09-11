#include "max.h"

MAX::MAX()
{
    //ctor
}

MAX::~MAX()
{
    //dtor
}

void MAX::Max(vector<double>& target)
// @func MAX: function to find the maximum value among some vector entries. It is a linear search, basically compare from the first element then run through the whole vector.
// @param vector<double>: the vector contains all the entries which we would like to find the maximum value from.
{
    _MaxValue = target.front();
    for(unsigned int i=1; i< target.size(); i++ )
    {
        if(_MaxValue<target[i])
        {
            _MaxIndex = static_cast<int>(i);
            _MaxValue = target[i];
        }
        else
        {
            continue;
        }
    }
}
