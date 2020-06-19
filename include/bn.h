/*
Defines the class of BN, which basically stores the structure by link and adjacent matrix,
and the score of the network with given structure.

*/
#ifndef BN_H
#define BN_H
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>


#include "node.h"
#include "cpd.h"
#include "data.h"
#include "edge.h"
#include "constraint.h"
using namespace std;

class BN
{
public:
    BN();
    // ctor/dtor/set/get functions:

    virtual ~BN();

    double Get_score() const
    {
        return _score;
    }

    void Set_score(double val)
    {
        _score = val;    // typically we will not assign a score to a network,however it makes sense to keep the mutator here..
    }

    bool Get_growing() const
    {
        return _growing;
    }

    void Set_growing(bool val)
    {
        _growing = val;    // typically we will not assign a score to a network,however it makes sense to keep the mutator here..
    }

    vector < NODE >& Get_allNodes()
    {
        return _allNodes;
    }

    NODE& Get_aNode(int id)
    {
        return _allNodes[id];
    }

    void Set_allNodes(vector < NODE >& val)
    {
        _allNodes = val;
    }

    vector < EDGE >& Get_edges()
    {
        return _edges;
    }

    EDGE& Get_anEdge(int id)
    {
        return _edges[id];
    }

    int Get_randSeed()
    {
        return _RandSeed;
    }

    void Set_randSeed(int seed)
    {
        _RandSeed = seed;
    }

    void Set_edges(vector < EDGE >& val)
    {
        _edges = val;
    }

    vector <int> & Get_dependencies()
    {
        return _dependencies;
    }

    int Get_dependencies(int id)
    {
        return _dependencies[id];
    }

    vector <int> & Get_numberOfNeighbors()
    {
        return _numberOfNeighbors;
    }

    int Get_numberOfNeighbors(int id)
    {
        return _numberOfNeighbors[id];
    }

    void Set_name(string name){
        _name = name;
    }

    string Get_name(){
        return _name;
    }

    // utility functions:
    void Initial(DATA&,double);

    void Set_arc( int, int, bool);
    // @func Set_arc: set an arc in the BN
    // @param int: starting node(fromNode)
    // @param int: pointing node(toNode)

    void BICscore();
    // BICscore: Used to initialize score.
    //void UpdateScore(const vector<CPD>&, int, int);

    void Update(DATA&,double, ostream&);
    // Update: smc updating
    // @DATA&: the data object
    // @double: temperature
    
    void HC(DATA&);
    // HC: hill climbing to ensure reaching the local optimal
    // @DATA&: data object

    string Write_Bnlearn_modelstring(DATA&);

private:
    string _name;
    double _score; // BIC score
    vector<NODE> _allNodes;
    vector<int> _dependencies;
    vector<int> _numberOfNeighbors;
    vector < EDGE > _edges;
    bool _growing;
    bool _oneNeb;
    int _RandSeed;

    // Some parameters that should not be accessed from outside
    // this chunk of paramters are used when there is no more nodes with more than one neighbor
    vector<int> _candOneNeb;
    double _totalScoreOneNeb;
    vector<double> _scoreOneNeb;

    bool _IsDAG(int , int);
    vector<int> _SearchNoneZero(vector<int>&);
    vector< int > _MinOverOne(vector<int>&);
    double _Summation(vector<double>&);// linear addition.
    void _SelectDirection(int,int,DATA&,double,ostream&);
    void _TripletSelection(vector<int>&,DATA&,double,bool,ostream&);
    int _IndexSearch(int);
    bool _ExistArc(int,int);
    vector<bool> _changedAcc;
};

#endif // BN_H
