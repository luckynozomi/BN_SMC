#include "constraint.h"
#include "utils.h"

CONSTRAINT::CONSTRAINT()
{
    //ctor
}

CONSTRAINT::CONSTRAINT(double val)
{

    _cutoff = val;

}

CONSTRAINT::~CONSTRAINT()
{
    //dtor
}

void CONSTRAINT::MIT(DATA& data, ofstream& mitfile)
{
    // Mutual information test for the discretized data is essentially Pearson's Chi-square test
    vector < int > & param = data.Get_param();
    for(int i=0; i< data.Get_NumberOfNodes(); i++)// the i-th node being tested with _correspondingNode
    {
        if(_correspondingNode!=i)
        {
            int obs = 0;
            for(int j=0; j < data.Get_NumberOfObservations(); ++j){
                if(data.Get_data(j, i) != 0 && data.Get_data(j, _correspondingNode) != 0)
                    obs += 1;
            }

            if(obs==0) continue; // If no observations in common, then we assume the nodes are independent

            int countTable[(param[i])][(param[_correspondingNode])];// construct the contingency table
            int coltotal[ (param[_correspondingNode] )];// column total(marginal)
            int rowtotal[ (param[i]) ];// row total(marginal)

            // This for loop initializes all the arrays.
            for(int col=0; col< param[_correspondingNode]; col++)
            {
                coltotal[col]=0;
                for(int row=0; row< param[i]; row++)
                {
                    rowtotal[row]=0;
                    countTable[row][col]=0;
                }
            }

            for(int j=0; j< data.Get_NumberOfObservations(); j++)
            {
                if(data.Get_data(j, i) != 0 && data.Get_data(j, _correspondingNode) != 0)
                    countTable[ (data.Get_data(j,i)-1) ][ (data.Get_data(j,_correspondingNode)-1) ]++;
                    // count the observations into the contingency table.
            }

            // column total:
            for(int col=0; col< param[_correspondingNode]; col++)
            {
                for(int row=0; row< param[i]; row++)
                {
                    coltotal[col] += countTable[row][col];
                }
            }

            // row total:
            for(int row=0; row<param[i]; row++)
            {
                for(int col=0; col<param[_correspondingNode]; col++)
                {
                    rowtotal[row] += countTable[row][col];
                }
            }

            // Conduct the MI test
            double mi=0.0;
            for(int row=0; row<param[i]; row++)
            {
                for(int col=0; col<param[_correspondingNode]; col++)
                {
                    mi += 2.0*countTable[row][col]*log( (1.0e-12) + obs*countTable[row][col]/( rowtotal[row]*coltotal[col]+(1.0e-12) ) );
                }
            }
            boost::math::chi_squared thedist( (param[i]-1)*(param[_correspondingNode]-1) );
            // define a chi-square test with degree of freedom (a-1)(b-1) where a is number of categories of the tested node and b number of categories of the corresponding node.
            double pValue = cdf( complement( thedist, mi ) );
            double data_pValue = pValue;
            double prior_pValue = 0;
            double posterior_pValue = 0;
            bool exists_prior = false;
            if(data.Exists_prior(i, _correspondingNode) || data.Exists_prior(_correspondingNode, i)) 
            {
                exists_prior = true;
                vector<double> pvals;
                pvals.push_back(pValue);
                if(data.Exists_prior(i, _correspondingNode)){
                    pvals.push_back(data.Get_prior_pval(i, _correspondingNode));
                    prior_pValue = data.Get_prior_pval(i, _correspondingNode);
                }
                if(data.Exists_prior(_correspondingNode, i)){
                    pvals.push_back(data.Get_prior_pval(_correspondingNode, i));
                    prior_pValue = data.Get_prior_pval(_correspondingNode, i);
                }
                pValue = combine_pval(pvals);
            }
            posterior_pValue = pValue;
            mitfile << data.Get_node_name(i) << "," << data.Get_node_name(_correspondingNode) << "," << prior_pValue << "," << data_pValue << "," << posterior_pValue << "," << 1.0-pValue << "," << (pValue<_cutoff) << endl;
            _testScore.push_back(1.0-pValue);
            if(pValue < _cutoff)// if the p-value is smaller than the cutoff, then the tested node is dependent with the corresponding node.
            {
                _InsertOrderedScore(i,1.0-pValue);
            }
            else
            {
                continue;
            }
        }
        else
        {
            continue;
        }
    }
}


int CONSTRAINT::_IndexSearch(int target)
{
    int L = 0;
    int U = _dependentNode.size()-1;
    int tar = 0;
    while( (_dependentNode[L]!=target) && (_dependentNode[U]!=target) )
    {
        if(target<_dependentNode[(L+U)/2])
        {
            U = (L+U)/2;
        }
        else
        {
            L = (L+U)/2;
        }
    }
    if (_dependentNode[L]==target)
    {
        tar = L;
    }
    else if(_dependentNode[U]==target)
    {
        tar = U;
    }
    else
    {
        cerr<<"Index search failed. It should not happen really."<<endl;
    }
    return tar;
}

void CONSTRAINT::_InsertOrderedScore(int node, double score)
{
    int L = 0;
    int U = _orderedTestScore.size()-1;
    if(_orderedTestScore.empty())
    {
        _orderedTestScore.push_back(score);
        _dependentNode.push_back(node);
    }
    else if(_orderedTestScore.size()==1 )
    {
        if(score >= _orderedTestScore[0])
        {
            _orderedTestScore.insert(_orderedTestScore.begin(),score);
            _dependentNode.insert(_dependentNode.begin(),node);

        }
        else
        {
            _orderedTestScore.push_back(score);
            _dependentNode.push_back(node);
        }
    }
    else
    {
        if(score>=_orderedTestScore[L])
        {
            _orderedTestScore.insert(_orderedTestScore.begin(),score);
            _dependentNode.insert(_dependentNode.begin(),node);

        }
        else if(score<=_orderedTestScore[U])
        {
            _orderedTestScore.push_back(score);
            _dependentNode.push_back(node);
        }
        else
        {
            while( U-L>1 )
            {
                if(score>_orderedTestScore[(L+U)/2])
                {
                    U = (L+U)/2;
                }
                else
                {
                    L = (L+U)/2;
                }
            }
            _orderedTestScore.insert(_orderedTestScore.begin()+L,score);
            _dependentNode.insert(_dependentNode.begin()+L,node);
        }
    }
}
