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
    _initBN.Set_name(_name);
    _initBN.Initial(data,cutoff);

    for(int i=0;i<chains;i++)
    {
        _bestBNs[i]=_initBN;
        _bestBICs.push_back(_initBN.Get_score());
    }
}


void SMC::Update(DATA& data, int numberOfIter, int numberOfChains,int NofCores, double temper,int depth, vector<int> rand_seeds)
{
    omp_set_num_threads(NofCores);
#pragma omp parallel for// parallel computing using OPENMP
    for(int i=0;i<numberOfChains;i++)
    {

        BN runBN;
        srand(rand_seeds[i]);
        runBN.Set_randSeed(rand_seeds[i]);
        runBN = _initBN;
        int k = 0;
        int last_edge_idx = -1;
        ofstream log_file;
        log_file.open(_name + "/logs/" + to_string(rand_seeds[i]) + ".txt", ios::app);

        while( (k < numberOfIter) && (runBN.Get_growing()) )
        {
            vector<EDGE> edges = runBN.Get_edges();
            for(int edge_idx=last_edge_idx+1; edge_idx!=edges.size(); ++edge_idx){
                log_file << "Added edge: "<< data.Get_node_name(edges[edge_idx].Get_dLink().first) << " --> " << data.Get_node_name(edges[edge_idx].Get_dLink().second) << endl;
                int from_edge = edges[edge_idx].Get_dLink().first;
                int to_edge = edges[edge_idx].Get_dLink().first;
            }
            last_edge_idx = edges.size() - 1;
            k++;
            runBN.Update(data,temper, log_file);

            if(runBN.Get_score()>_bestBICs[i])
            {
                _bestBICs[i] = runBN.Get_score();
                _bestBNs[i] = runBN;
            }
        }
        log_file.close();
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



void SMC::Summary(string outpath,int chains,int nodes,int depth, DATA& data, vector<int> rand_seeds)
{
    int max_ind = _MaximumIndSearch(_bestBICs);
    double best_posterior = _bestBNs[max_ind].Get_score();
    double eps = 1e-5;
    ofstream summary_file;
    double ll_prior = 0.0;
    double ll_posterior = 0.0;
    double ll_data = 0.0;
    summary_file.open(outpath+"/summaries/bn.tsv");
    summary_file << std::boolalpha;
    summary_file << "best BN?\tindex\trand seed\tprior\tdata\tposterior\tmodel string" << endl;
    for(int chain_idx=0; chain_idx!=chains; ++chain_idx){
        BN& this_BN = _bestBNs[chain_idx];
        for(int node_idx=0; node_idx!=nodes; ++node_idx){
            NODE& this_node = this_BN.Get_aNode(node_idx);
            //prior log-likelihood, data log-likelihood and posterior log-likelihood
            vector<double> BICs = this_node.UpdateBIC(data);
            ll_prior += BICs[0];
            ll_data += BICs[1];
            ll_posterior += BICs[2];
        }
        summary_file << (abs(best_posterior - ll_posterior) <= eps) << "\t";
        summary_file << chain_idx << "\t" << setw(12) << rand_seeds[chain_idx] << "\t";
        summary_file << setw(12) << ll_prior << "\t" << setw(12) << ll_data << "\t" <<  setw(12) << ll_posterior << "\t";
        summary_file << this_BN.Write_Bnlearn_modelstring(data) << endl;
        ll_data = 0;
        ll_prior = 0;
        ll_posterior = 0;
    }
    summary_file.close();


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
    outFile.open((outpath+"/summaries/summary.txt").c_str(),fstream::app);
    outFile<<_bestBICs[max_ind]<<endl;
    outFile<<_averageCount<<endl;
    for(int i=0;i<_maxBicsRound.size();i++)
    {
        outFile<<_maxBicsRound[i]<<",";
    }
    outFile<<endl;
    outFile.close();

    outFile_scores.open((outpath+"/summaries/scores.R").c_str());
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

    outFile_best.open((outpath+"/summaries/best_edges.R").c_str());
    outFile_best<<"res.bn<-empty.graph(nam)"<<endl;
    for(unsigned int i=0;i<_bestBNs[max_ind].Get_edges().size();i++)
    {
        outFile_best<<"res.bn<-set.arc(res.bn,nam["<<_bestBNs[max_ind].Get_anEdge(i).Get_dLink().first+1<<"],nam["
            <<_bestBNs[max_ind].Get_anEdge(i).Get_dLink().second+1<<"])"<<endl;

    }
    outFile_best.close();

    ofstream outFile_bestbn;
    outFile_bestbn.open((outpath+"/summaries/best_BN.txt"));
    outFile_bestbn << _bestBNs[max_ind].Write_Bnlearn_modelstring(data) << endl;
    outFile_bestbn.close();

    outFile_sum.open((outpath+"/summaries/sum_edges.txt").c_str());

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

    outFile_average.open((outpath+"/summaries/average_edges.txt").c_str());
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
