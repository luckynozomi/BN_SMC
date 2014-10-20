#include "constraint.h"

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

void CONSTRAINT::MIT(DATA& data)
{
    // Mutual information test for the discretized data is essentially Pearson's Chi-square test

    vector < int > & param = data.Get_param();
    for(int i=0; i< data.Get_NumberOfNodes(); i++)// the i-th node being tested with _correspondingNode
    {

        if(_correspondingNode!=i)
        {

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

            //vector<vector< int> > temp_data = data.Get_data();
            int obs=data.Get_NumberOfObservations();
            for(int j=0; j< obs; j++)
            {
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

            //row total:
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
//                   mi += (2.0*static_cast<double>(countTable[row][col]))*log( (1.0e-12) + static_cast<double>(obs)*countTable[row][col]/( rowtotal[row]*coltotal[col]+(1.0e-12) ) );
                    mi += 2.0*countTable[row][col]*log( (1.0e-12) + obs*countTable[row][col]/( rowtotal[row]*coltotal[col]+(1.0e-12) ) );
                }
            }
            boost::math::chi_squared thedist( (param[i]-1)*(param[_correspondingNode]-1) );
            // define a chi-square test with degree of freedom (a-1)(b-1) where a is number of categories of the tested node and b number of categories of the corresponding node.
            double pValue = cdf( complement( thedist, mi ) );
            _testScore.push_back(1.0-pValue);
            if(pValue < _cutoff)// if the p-value is smaller than the cutoff, then the tested node is dependent with the corresponding node.
            {
                _InsertOrderedScore(i,1.0-pValue);
//                _testScore.push_back(1.0-pValue);
            }
            else
            {
                continue;
            }
//                _dependentNode.push_back(i);


        }
        else
        {
            continue;
        }

    }

}

//void CONSTRAINT::clear_rank(int theNode)
//{
//    _testScore[theNode]=0.0;
//}
/*void CONSTRAINT::ClearRank(int idOrPos,bool ifdirect)
{
    if(ifdirect)
    {
        _dependentNode.erase(_dependentNode.begin()+idOrPos);
        _testScore.erase(_testScore.begin()+idOrPos);
    }
    else
    {
        int ind = _IndexSearch(idOrPos);
        _dependentNode.erase(_dependentNode.begin()+ind);
        _testScore.erase(_testScore.begin()+ind);

    }
}
*/

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
