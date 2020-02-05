/*
Defines the class which holds the data, so far only work with comma separated integer file.

*/

#ifndef DATA_H
#define DATA_H
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <boost/tokenizer.hpp>

using namespace std;

class DATA
{
public:
    // ctor/dtor/set/get functions:
    DATA();

    DATA(int,int);

    virtual ~DATA();

    int Get_NumberOfObservations() const
    {
        return _NumberOfObservations;
    }

    void Set_NumberOfObservations(int val)
    {
        _NumberOfObservations = val;
    }

    int Get_NumberOfNodes() const
    {
        return _NumberOfNodes;
    }

    void Set_NumberOfNodes(int val)
    {
        _NumberOfNodes = val;
    }

    vector<vector<int > >& Get_data()
    {
        return _data;
    }

    int Get_data(int row, int col)
    {
        return _data[row][col];
    }

    void Set_data(vector<vector<int > > & val)
    {
        _data = val;
    }

    vector<int>& Get_param()
    {
        return _param;
    }

    int Get_param(int id)
    {
        return _param[id];
    }

    double Get_prior_prob(int from_node, int to_node)
    {
        return _prior_prob[from_node][to_node];
    }

    double Get_prior_pval(int from_node, int to_node)
    {
        return _prior_pval[from_node][to_node];
    }

    bool Exists_prior(int from_node, int to_node)
    {
        return Get_prior_pval(from_node, to_node) != -1.0;
    }

    bool Exists_prior()
    {
        return _hasPrior;
    }

    // Utility function:
    void ReadData(const char*);
    // @func ReadCSV: read a csv file into _data
    // @param: const char*: file name with full path.

    void ReadParam(const char*);

    void ReadPrior(const char*);
private:
    int _NumberOfObservations; // should be set during initialization
    int _NumberOfNodes; // also should be set during initialization
    bool _hasPrior; // true if there is prior
    vector < vector < int > > _data; // vector that stores the data for further use.
    vector<int> _param; // Numbers of parameters for each node.
    vector< vector<double> > _prior_prob; // Prior knowledge of the data.
    vector< vector<double> > _prior_pval; // Prior knowledge of the data.
};

#endif // DATA_H
