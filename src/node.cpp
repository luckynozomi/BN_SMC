#include "node.h"
#include "utils.h"

NODE::NODE()
{
    //ctor
}

NODE::~NODE()
{
    //dtor
}

void NODE::UpdateAncester(int fromNode, set<int>& acFromNode)
{
    _ancester.insert(fromNode); // first insert the fromNode
    if(!acFromNode.empty() ) // then if the ancester of the fromNode is not empty add the ancester to the ancester of this node.
    {
        _ancester.insert( acFromNode.begin(), acFromNode.end() );
    }
}

void NODE::UpdateParent(int fromNode)
{
    _parent.push_back( fromNode ); // The algorithm guarantees that the parent node should not have been in the parent vector of this node already so we just add the parent id to the parent vector.
}

void NODE::UpdateDescendant(int toNode, set<int>& deToNode) // similar as update the ancester
{
    _descendant.insert( toNode );
    if(!deToNode.empty())
    {
        _descendant.insert( deToNode.begin(), deToNode.end() );
    }
}

void NODE::UpdateChild(int toNode) //similar as update parent.
{
    _child.insert( toNode );
}

void NODE::UpdateBIC(DATA& data)
{
    _correspondingCPD.UpdateBIC(_parent, data);
    _scoreContribution = _correspondingCPD.Get_BIC();

}
void NODE::ClearRank(int idOrPos,bool ifdirect)
{
    if(ifdirect)
    {
        _neighbors.erase(_neighbors.begin()+idOrPos);
        _scoreOfNeighbors.erase(_scoreOfNeighbors.begin()+idOrPos);
    }
    else
    {
        int ind = _IndexBiSearch(idOrPos);
        _neighbors.erase(_neighbors.begin()+ind);
        _scoreOfNeighbors.erase(_scoreOfNeighbors.begin()+ind);
    }
}

void NODE::CMIT(vector<NODE>& allnodes, DATA& data,double cutoff)
{
    vector<int> tmpNeighbor = _neighbors;
    for(int i=0; i<tmpNeighbor.size();i++)
    {
        if(tmpNeighbor[i]!=-1) // tmpNeighbor[i] is going to be conditioned on.
        {
            vector<int> currentNeighbor = tmpNeighbor;
            sort(currentNeighbor.begin(),currentNeighbor.end());
            vector<int> neighborOfTarget = allnodes[tmpNeighbor[i]].Get_neighbors();
            sort(neighborOfTarget.begin(),neighborOfTarget.end());
            vector<int> intNeighbor(currentNeighbor.size()+neighborOfTarget.size());
            vector<int> :: iterator it;
            it= set_intersection(currentNeighbor.begin(),currentNeighbor.end(),
                                neighborOfTarget.begin(),neighborOfTarget.end(),
                                intNeighbor.begin());
            intNeighbor.resize(it-intNeighbor.begin());
            if( !(intNeighbor.empty()) )
            {
                for(int j=0;j<intNeighbor.size();j++)
                {
                    if(intNeighbor[j]!=-1&& !(_CMIT(_nodeID,intNeighbor[j],tmpNeighbor[i],data,cutoff)) )
                    {
                        for(int k=0; k<tmpNeighbor.size();k++)
                        {
                            if(tmpNeighbor[k]==intNeighbor[j])
                            {
                                tmpNeighbor[k] = -1;
                            }

                        }

                    }
                    else
                    {
                        continue;
                    }
                }
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
    _neighbors.clear();
    for(int i=0;i<tmpNeighbor.size();i++)
    {
        if(tmpNeighbor[i]!=-1)
        {
            _neighbors.push_back(tmpNeighbor[i]);
            _scoreOfNeighbors.push_back(_correspondingConstraint.Get_testScore(tmpNeighbor[i]));
        }
        else
        {
            continue;
        }
    }
}

void NODE::SymmetryCorrection(vector<NODE>& allnodes)
{

    for(int i = _neighbors.size()-1; i>=0;i--)
    {
        vector<int>::iterator it;
        it = find((allnodes[ _neighbors[i] ].Get_neighbors()).begin(),(allnodes[ _neighbors[i] ].Get_neighbors()).end(), _nodeID);
        if(it == (allnodes[ _neighbors[i] ].Get_neighbors()).end() )
        {
            _neighbors.erase(_neighbors.begin()+i);

            _scoreOfNeighbors.erase(_scoreOfNeighbors.begin()+i);
        }
        else
        {
            continue;
        }
    }

}

void NODE::SortNeighbor()
{
    vector<int> tmpNeb = _neighbors;
    vector<double> tmpScore = _scoreOfNeighbors;
    sort(_neighbors.begin(),_neighbors.end());
    for(int i=0; i<_neighbors.size();i++)
    {
        int tmpInd = _IndexBiSearch(tmpNeb[i]);
        _scoreOfNeighbors[tmpInd] = tmpScore[i];
    }
}

void NODE::UpdateResult(NODE& copyNode)
{
    _correspondingCPD = copyNode.Get_correspondingCPD();
    _neighbors = copyNode.Get_neighbors();
    _scoreOfNeighbors = copyNode.Get_scoreOfNeighbors();
    _scoreContribution = copyNode.Get_scoreContribution();
}


void NODE::RemoveParent(int val)
{
    vector<int>::iterator it;
    it = find(_parent.begin(),_parent.end(),val);
    if(it != _parent.end())
    {
        _parent.erase(it);
    }
}


void NODE::RemoveChild(int val)
{
    _child.erase(val);

}

// Conditional MIT test of nodeA & nodeB, given nodeC.
// Return true if dependent (pValue < cutoff)
bool NODE::_CMIT(int nodeA, int nodeB, int nodeC, DATA& data,double cutoff)
{
    int countTable[data.Get_param(nodeC)][data.Get_param(nodeB)  ][data.Get_param(nodeA)];
    int p_xz[ data.Get_param(nodeC) ][ data.Get_param(nodeA)];
    int p_yz[data.Get_param(nodeC)][data.Get_param(nodeB)];
    int p_z[ data.Get_param(nodeC) ];

    //Intialize the counting tables.
    for(int con = 0; con< data.Get_param(nodeC);con++)
    {
        p_z[con]=0;
        for(int x=0; x< data.Get_param(nodeA);x++)
        {
            p_xz[con][x]=0;
            for(int y=0; y<data.Get_param(nodeB); y++)
            {
                p_yz[con][y]=0;
                countTable[con][y][x]=0;
            }
        }
    }

    int obs = 0;
    for(int j=0; j < data.Get_NumberOfObservations(); ++j){
        if(data.Get_data(j, nodeA) != 0 && data.Get_data(j, nodeB) != 0 && data.Get_data(j, nodeC) != 0)
            obs += 1;
    }

    if(obs==0) return false; // Return independent if no data in common

    // Fill out the countTable
    for(int row=0; row<data.Get_NumberOfObservations(); row++)
    {        
        if(data.Get_data(row, nodeA) != 0 && data.Get_data(row, nodeB) != 0 && data.Get_data(row, nodeC) != 0)
            countTable[ (data.Get_data(row,nodeC)-1) ][ (data.Get_data(row,nodeB)-1) ][ (data.Get_data(row,nodeA)-1) ]++;
    }

    // marginalization:
    // p_xz and p_z
    for(int con = 0; con< data.Get_param(nodeC);con++)
    {
        for(int x=0; x< data.Get_param(nodeA);x++)
        {
            for(int y=0; y<data.Get_param(nodeB); y++)
            {
                p_xz[con][x]+= countTable[con][y][x];
            }
            p_z[con]+=p_xz[con][x];
        }
    }

    // p_yz
    for(int con = 0; con< data.Get_param(nodeC);con++)
    {
        for(int y=0; y< data.Get_param(nodeB);y++)
        {
            for(int x=0; x<data.Get_param(nodeA); x++)
            {
                p_yz[con][y]+= countTable[con][y][x];
            }
        }
    }

    // calculate the chi-square statistic
    double cmi=0.0;
    for(int con = 0; con< data.Get_param(nodeC);con++)
    {
        for(int y=0; y< data.Get_param(nodeB);y++)
        {
            for(int x=0; x<data.Get_param(nodeA); x++)
            {
                cmi+=2.0*countTable[con][y][x]*log( (1.0e-12)+ (static_cast<double>(p_z[con])*countTable[con][y][x] ) / (static_cast<double>(p_xz[con][x])*p_yz[con][y] +(1.0e-12) ) );
            }
        }
    }

    boost::math::chi_squared thedist( (data.Get_param(nodeA)-1)*(data.Get_param(nodeB)-1)*(data.Get_param(nodeC)) );
    double pValue = cdf( complement( thedist, cmi ) );
    if(data.Exists_prior(nodeA, nodeB) || data.Exists_prior(nodeB, nodeA)){
        vector<double> pvals;
        pvals.push_back(pValue);
        if(data.Exists_prior(nodeA, nodeB))
            pvals.push_back(data.Get_prior_pval(nodeA, nodeB));
        if(data.Exists_prior(nodeB, nodeA))
            pvals.push_back(data.Get_prior_pval(nodeB, nodeA));
        pValue = combine_pval(pvals);
    }
    if(pValue<cutoff)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int NODE::_IndexSearch(int target)
{
    int tar = -1;
    for(int i=0; i<_neighbors.size();i++)
    {
        if(_neighbors[i]== target)
        {
            tar=i;
            break;
        }
    }
    return tar;
}

int NODE::_IndexBiSearch(int target)
{
    int L = 0;
    int U = _neighbors.size()-1;
    int tar = 0;
    while( (_neighbors[L]!=target) && (_neighbors[U]!=target) )
    {
        if(target<_neighbors[(L+U)/2])
        {
            U = (L+U)/2;
        }
        else
        {
            L = (L+U)/2;
        }
    }
    if (_neighbors[L]==target)
    {
        tar = L;
    }
    else if(_neighbors[U]==target)
    {
        tar = U;
    }
    else
    {
        cerr<<"Index search failed. It should not happen really."<<endl;
    }
    return tar;
}

