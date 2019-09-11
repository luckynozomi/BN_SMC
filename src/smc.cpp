#include "smc.h"

SMC::SMC()
{
    //ctor
}

SMC::~SMC()
{
    //dtor
}
SMC::SMC(int nodes,int chains)
{
    vector<BN> helper(chains);
    _bestBNs = helper;
    _bestBICs.reserve(chains);

    _averageCount = 0;
    _edgeMatrixSum.reserve(nodes);
    _edgeMatrixOptimal.reserve(nodes);
    _edgeMatrixAverage.reserve(nodes);
    vector<int> thehelper(nodes);
    for(int i =0; i<nodes;i++)
    {
        _edgeMatrixOptimal.push_back(thehelper);
        _edgeMatrixSum.push_back(thehelper);
        _edgeMatrixAverage.push_back(thehelper);

    }
}


void SMC::Initialize(DATA& data, double cutoff,int chains)
{
    _initBN.Initial(data,cutoff);

    for(int i=0;i<chains;i++)
    {
        _bestBNs[i]=_initBN;
        _bestBICs.push_back(_initBN.Get_score());
    }
}


void SMC::Update(DATA& data, int numberOfIter, int numberOfChains,int NofCores, double temper,int depth)
{
    omp_set_num_threads(NofCores);
#pragma omp parallel for// parallel computing using OPENMP
    for(int i=0;i<numberOfChains;i++)
    {

        srand(i);
        BN runBN;

        runBN = _initBN;
        int k = 0;
        while( (k < numberOfIter) && (runBN.Get_growing()) )
        {

            k++;
            runBN.Update(data,temper);

            if(runBN.Get_score()>_bestBICs[i])
            {
                _bestBICs[i] = runBN.Get_score();
                _bestBNs[i] = runBN;
            }
        }
    }
    int round = 0;
    _bicScores.reserve(depth+1);
    _bicScores.push_back(_bestBICs);
    while(round < depth)
    {
        round++;

        #pragma omp parallel for
        for(int i = 0; i<numberOfChains;i++)
        {
            _bestBNs[i].HC(data);
            _bestBICs[i] = _bestBNs[i].Get_score();
        }
        _bicScores.push_back(_bestBICs);
    }
}


void SMC::DFSummary(string outpath,int numberNodes)
{
    ofstream DFresult;
    DFresult.open((outpath+"_DFsummary.R").c_str());
    DFresult<<"DFcon<-list(NULL)"<<endl;
    for(int i=0; i<numberNodes;i++)
    {
        if(_initBN.Get_numberOfNeighbors(i)>0)
        {
            DFresult<<"DFcon[["<<i+1<<"]]<-c(";
            for(int j=0;j< _initBN.Get_aNode(i).Get_neighbors().size()-1 ;j++)
            {
                DFresult<< _initBN.Get_aNode(i).Get_neighbors(j)+1<<",";
            }
            DFresult<< _initBN.Get_aNode(i).Get_neighbors(_initBN.Get_aNode(i).Get_neighbors().size()-1)+1;
            DFresult<<")"<<endl;
        }
        else
        {
            continue;
        }
    }
    DFresult.close();
}



void SMC::Summary(string outpath,int chains,int nodes,int depth)
{
    int max_ind = _MaximumIndSearch(_bestBICs);

    // compute the egdes ever sampled
    for(int i=0;i<chains;i++)
    {
        for(unsigned j=0;j<_bestBNs[i].Get_edges().size();j++)
        {
            _edgeMatrixSum[_bestBNs[i].Get_anEdge(j).Get_dLink().first][_bestBNs[i].Get_anEdge(j).Get_dLink().second]++;
        }
    }

    //compute the edsges sampled for the network with a BIC score inside first standard deviation.
    Statistics sta;
    sta.Summary(_bestBICs);
    double _cutOffBic = _bestBICs[max_ind]-sta.Get_std();
    for(int i=0;i<chains;i++)
    {
        if(_bestBICs[i]>=_cutOffBic)
        {
            _averageCount++;
            for(unsigned j=0;j<_bestBNs[i].Get_edges().size();j++)
            {
                _edgeMatrixAverage[_bestBNs[i].Get_anEdge(j).Get_dLink().first][_bestBNs[i].Get_anEdge(j).Get_dLink().second]++;
            }
        }
    }

    // outputs
    ofstream outFile,outFile_best,outFile_sum,outFile_scores,outFile_average;
    outFile.open((outpath+"_summary.txt").c_str(),fstream::app);
    outFile<<_bestBICs[max_ind]<<endl;
    outFile<<_averageCount<<endl;
    for(int i=0;i<_maxBicsRound.size();i++)
    {
        outFile<<_maxBicsRound[i]<<",";
    }
    outFile<<endl;
    outFile.close();

    outFile_scores.open((outpath+"_scores.R").c_str());
    outFile_scores<<"scores<-NULL"<<endl;
    for(int j=0;j<=depth;j++)
    {
        outFile_scores<<"scores[["<< (j+1) <<"]]<-c(";
        for(int i=0;i<chains-1;i++)
        {
            outFile_scores<<_bicScores[j][i] <<",";
        }
        outFile_scores<<_bicScores[j][chains-1]<<")"<<endl;
    }
    outFile_scores.close();

    outFile_best.open((outpath+"_best_edges.R").c_str());
    outFile_best<<"res.bn<-empty.graph(nam)"<<endl;
    for(unsigned int i=0;i<_bestBNs[max_ind].Get_edges().size();i++)
    {
        outFile_best<<"res.bn<-set.arc(res.bn,nam["<<_bestBNs[max_ind].Get_anEdge(i).Get_dLink().first+1<<"],nam["
            <<_bestBNs[max_ind].Get_anEdge(i).Get_dLink().second+1<<"])"<<endl;

    }
    outFile_best.close();

    outFile_sum.open((outpath+"_sum_edges.txt").c_str());

    outFile_sum<<"From, To, Count"<<endl;
    for(int i = 0;i<nodes;i++)
    {
        for(int j = 0;j<nodes;j++)
        {
            if(_edgeMatrixSum[i][j]>0)
            {
                outFile_sum<<i+1<<","<<j+1<<","<<_edgeMatrixSum[i][j]<<endl;
            }
        }

    }
    outFile_sum.close();

    outFile_average.open((outpath+"_average_edges.txt").c_str());
    outFile_average<<"From, To, Count"<<endl;
    for(int i = 0;i<nodes;i++)
    {
        for(int j = 0;j<nodes;j++)
        {
            if(_edgeMatrixAverage[i][j]>0)
            {
                outFile_average<<i+1<<","<<j+1<<","<< static_cast<double>(_edgeMatrixAverage[i][j])/_averageCount <<endl;
            }
        }
    }
    outFile_average.close();
}

int SMC::_MaximumIndSearch(vector<double>& target)
{
    int index = 0;
    double maximum = target[index];
    for(unsigned int i=1;i<target.size();i++)
    {
        if(target[i] >maximum)
        {
            maximum = target[i];
            index = i;
        }
        else
        {
            continue;
        }
    }
    return index;
}

bool SMC::_MyGreatThan(double a, double b)
{
    return a>b;
}
