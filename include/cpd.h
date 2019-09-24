/*
Define the class conditional probability distribution( so far it is probably more appropriate to call it CPT).
And we will use MLE so far, The Dirichlet prior bayesian estimator will come later.

*/


#ifndef CPD_H
#define CPD_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include "data.h"
using namespace std;

class CPD
{
public:
    CPD();

    virtual ~CPD();

    int Get_correspondingNode() const
    {
        return _correspondingNode;
    }

    void Set_correspondingNode(int val)
    {
        _correspondingNode = val;
    }

    vector< double >& Get_probability()
    {
        return _probability;
    }

    void Set_probability(vector< double > & val)
    {
        _probability = val;
    }

    vector< int >& Get_conditionCount()
    {
        return _conditionCount;
    }

    void Set_conditionCount(vector<int > & val)
    {
        _conditionCount = val;
    }

    int Get_numberOfParentConfigure() const
    {
        return _numberOfParentConfigure;
    }

    void Set_numberOfParentConfigure(int val)
    {
        _numberOfParentConfigure = val;
    }

    int Get_numberOfParameter() const
    {
        return _numberOfParameter;
    }

    void Set_numberOfParameter(int val)
    {
        _numberOfParameter = val;
    }

    double Get_BIC() const
    {
        return _BIC;
    }

    double Get_totalParam()
    {
        return _totalParam;
    }

    // utility functions:
    void UpdateCPD(const vector<int>& ,DATA&);
    //@func update the conditional probability table according to the local strucutre.
    //@param const Node: this node(or say the "toNode")
    //@param const DATA: the data
    //@param const vector<Node>: the vector of all nodes

    void UpdateBIC();

private:
    int _correspondingNode;// The corresponding node id.
    vector<int> _conditionCount;// count values of each case for the corresponding node
    vector<double> _probability;// the probability for each case of the corresponding condition table
    int _numberOfParentConfigure;
    int _numberOfParameter;// number of categories in the corresponding node.
    double _BIC;
    int _totalParam;
    
    int* _CumProduct(const int*, const int,int* );
    //@func )cumProduct: returns an array that returns the cumproduct. please see example part(in cpd.cpp) for more information
    //@param const int*: an array contains number of categories of each node involved in the cpd(parents  and this node)
    //@param const int: total number of nodes involved(parents plus node)
    //@param int* : leading by 1 array stores the result.

    int _IndexSearchTarget(vector<int>&, int*);
    //@func _indexSearchTarget: return the index of a specified searching target, used to count the conditionCount
    //@param vector<int>: the searching targtet(one observation with corresponding parent and this node values)
    //@param int*: the index array(please see example in cpd.cpp for more information
};

#endif // CPD_H
