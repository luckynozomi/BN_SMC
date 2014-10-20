#ifndef MAX_H
#define MAX_H
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include <stdlib.h>

using namespace std;

class MAX
{
public:
    MAX();
    virtual ~MAX();
    int Get_MaxIndex()
    {
        return _MaxIndex;
    }
    void Set_MaxIndex(int val)
    {
        _MaxIndex = val;
    }
    double Get_MaxValue()
    {
        return _MaxValue;
    }
    void Set_MaxValue(double val)
    {
        _MaxValue = val;
    }
    void Max(vector<double>& );
protected:
private:
    int _MaxIndex;
    double _MaxValue;

};

#endif // MAX_H
