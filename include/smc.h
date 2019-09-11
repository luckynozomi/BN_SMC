#ifndef SMC_H
#define SMC_H
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <omp.h>
#include "bn.h"
#include "data.h"
#include "statistics.h"
#include <string.h>
class SMC
{
public:
    SMC();

    SMC(int,int);

    virtual ~SMC();

    BN& Get_initBN()
    {
        return _initBN;
    }

    void Set_initBN(BN val)
    {
        _initBN = val;
    }

    vector<BN>& Get_bestBNs()
    {
        return _bestBNs;
    }

    BN& Get_bestBNs(int id)
    {
        return _bestBNs[id];
    }

    vector<double>& Get_bestBICs()
    {
        return _bestBICs;
    }

    double Get_bestBICs(int id)
    {
        return _bestBICs[id];
    }

    void Initialize(DATA&, double,int);

    void Update(DATA&,int,int,int,double,int);
    
    void DFSummary(string,int);
    
    void Summary(string,int,int,int);

private:
    BN _initBN;
    vector<BN> _bestBNs;
    vector<double> _bestBICs;
    int _averageCount;
    vector<double> _maxBicsRound;
    vector< vector < double > > _bicScores;
    vector< vector <int> > _edgeMatrixSum; //diagnositic using
    vector< vector <int> > _edgeMatrixAverage; //diagnositic using
    vector< vector <int> > _edgeMatrixOptimal;

    int _MaximumIndSearch(vector<double>&);
    bool _MyGreatThan(double,double);

};

#endif // SMC_H
