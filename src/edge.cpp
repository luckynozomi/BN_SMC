#include "edge.h"

EDGE::EDGE()
{
    //ctor
}

EDGE::~EDGE()
{
    //dtor
}

bool EDGE::Matched(int fromNode,int toNode)
{
    bool val=false;
    if(_dLink.first==fromNode && _dLink.second == toNode)
    {
        val=true;
    }
    return val;
}
