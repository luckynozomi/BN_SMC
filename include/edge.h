#ifndef EDGE_H
#define EDGE_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
using namespace std;
class EDGE
{
public:
    EDGE();

    virtual ~EDGE();
    
    pair < int, int > Get_dLink()
    {
        return _dLink;
    }
    
    void Set_dLink(pair < int, int > val)
    {
        _dLink = val;
    }
    
    void Set_edge(int fromNode, int toNode)
    {
        _dLink.first = fromNode;
        _dLink.second = toNode;
    }
    
    double Get_strenth()
    {
        return _strenth;
    }
    
    void Set_strenth(double val)
    {
        _strenth = val;
    }
    
    bool Matched(int,int);

private:
    pair < int, int > _dLink;
    double _strenth;
};

#endif // EDGE_H
