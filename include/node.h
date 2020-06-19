/*
Define the class of a node (vertex) , essential part of a graph.

*/

#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>

#include "cpd.h"
#include "constraint.h"
using namespace std;

class NODE
{
public:
    NODE();

    virtual ~NODE();

    string Get_nameOfNode() const
    {
        return _nameOfNode;
    }

    void Set_nameOfNode(string val)
    {
        _nameOfNode = val;
    }

    int Get_nodeID() const
    {
        return _nodeID;
    }

    void Set_nodeID(int val)
    {
        _nodeID = val;
    }

    vector<int>& Get_parent()
    {
        return _parent;
    }

    void Set_parent(vector<int>& val)
    {
        _parent = val;
    }

    vector<int>& Get_neighbors()
    {
        return _neighbors;
    }

    int Get_neighbors(int val)
    {
        return _neighbors[val];
    }

    void Set_neighbors(vector<int> val)
    {
        _neighbors= val;
    }

    set<int>& Get_child()
    {
        return _child;
    }

    void Set_child(set<int> val)
    {
        _child = val;
    }

    set<int>& Get_ancester()
    {
        return _ancester;
    }

    void Set_ancester(set<int> val)
    {
        _ancester = val;
    }

    set<int>& Get_descendant()
    {
        return _descendant;
    }

    void Set_descendant(set<int> val)
    {
        _descendant = val;
    }

    double Get_probability() const
    {
        return _probability;
    }

    void Set_probability(double val)
    {
        _probability = val;
    }

    int Get_Assignment() const
    {
        return _assignment;
    }

    void Set_Assignment(int val)
    {
        _assignment = val;
    }

    double Get_scoreContribution() const
    {
        return _scoreContribution;
    }

    void Set_scoreContribution(double val)
    {
        _scoreContribution = val;
    }

    int Get_numberOfParameter() const
    {
        return _numberOfParameter;
    }

    void Set_numberOfParameter(int val)
    {
        _numberOfParameter = val;
    }

    void Set_possibleAssignment(vector<string>& val)
    {
        _possibleAssignment = val;
    }

    vector<string>& Get_possibleAssignment()
    {
        return _possibleAssignment;
    }

    CPD& Get_correspondingCPD()
    {
        return _correspondingCPD;
    }

    void Set_correspondingCPD(CPD val)
    {
        _correspondingCPD = val;
    }

    CONSTRAINT& Get_correspondingConstraint()
    {
        return _correspondingConstraint;
    }

    vector<double>& Get_scoreOfNeighbors()
    {
        return _scoreOfNeighbors;
    }

    double Get_scoreOfNeighbors(int ind)
    {
        return _scoreOfNeighbors[ind];
    }

    void InsertAncester(int node)
    {
        _ancester.insert(node);
    }

    void InsertDescendant(int node)
    {
        _descendant.insert(node);
    }

    // utility functions:
    void UpdateAncester(int, set<int>&);
    // @func UpdateAncester: update the ancester set of this node after add an arc.
    // @param int: the node id of the new parent of this node, so called the "fromNode".
    // @param set<int>: the ancester set of the "fromNode"

    void UpdateDescendant(int,set<int>&);
    // @func UpdateDescendant: update the Descendant set of this node after add an arc.
    // @param int: the node id of the new child of this node, so called the "toNode".
    // @param set<int>: the Descendant set of the "toNode"

    void UpdateParent( int );
    // @func UpdateParent: update the parent vector of this node
    // @param int: the new parent,so called "fromNode", id.

    void UpdateChild( int );
    // @func UpdateChild: update the child set of this node
    // @param int: the new child,socalled "toNode", id.

    void UpdateCPD(DATA& data)
    {
        _correspondingCPD.UpdateCPD(_parent,data);
    }

    vector<double> UpdateBIC(DATA& data);

    void ClearRank(int,bool);

    void CMIT(vector<NODE> &,DATA&,double, ofstream&);

    void SymmetryCorrection(vector<NODE>&, DATA&);

    void SortNeighbor();

    void UpdateResult(NODE&);

    void RemoveParent(int);

    void RemoveChild(int);

private:
    string _nameOfNode; // optional, the name of the node
    int _nodeID; // The internal id for the node, should be in the same order of the one of data/information file.
    CPD _correspondingCPD; // the cpd corresponding to this node
    CONSTRAINT _correspondingConstraint; // the contraints from MIT test.
    vector <int> _neighbors;
    vector <double> _scoreOfNeighbors;
    vector<int> _parent; //records the id of the parent of this node
    set<int> _child; // record child id of this node
    set<int> _ancester; //record the ancester id of this node
    set<int> _descendant;// record the descendent id of this node
    double _scoreContribution; // The BIC score of this node.
    int _numberOfParameter; // number of categories(parameters) of this node, which in most cases is obtained from the node information(parameter) file.
    
    // The following variables are used for inference only, they are not used in learning.
    double _probability; // this probability refer to the probability of a certain assignment used in inference.
    int _assignment; // the assignment of this node ontained from the inference
    
    // the next variable is used to determine whether copy all information
    vector<string> _possibleAssignment; // optional, the names of each category.

    int _IndexSearch(int target);
    bool _CMIT(int,int,int,DATA&,double,ofstream&);
    int _IndexBiSearch(int target);
    };

#endif // NODE_H
