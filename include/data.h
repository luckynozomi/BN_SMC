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
/*    vector<vector<int > >& Get_skeleton()
    {
        return _skeleton;
    }
    int Get_skeleton(int row, int col)
    {
        return _skeleton[row][col];
    }*/
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

    // Utility function:
    void ReadData(const char*);
    // @func ReadCSV: read a csv file into _data
    // @param: const char*: file name with full path.
    void ReadParam(const char*);
//    void ReadSkeleton(const char*);
private:
    int _NumberOfObservations;//should be set during initialization
    int _NumberOfNodes; // also should be set during initialization
    vector<vector<int > > _data; // vector that stores the data for further use.
    vector<int> _param;//Numbers of parameters for each node.
//    vector < vector< int > > _skeleton;

};

#endif // DATA_H
