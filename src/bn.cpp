/*
Implementing bn.h
*/


#include "bn.h"
BN::BN()
{

}

BN::~BN()
{
    //dtor
}



void BN::Initial(DATA& data, double cutoff)
{
    _allNodes.reserve(data.Get_NumberOfNodes());
    _score = 0.0;
    _dependencies.reserve(data.Get_NumberOfNodes());
    //cout<<data.Get_NumberOfNodes()<<endl;
    _numberOfNeighbors.reserve(data.Get_NumberOfNodes());

    _growing = true;
    _oneNeb = false;

    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        NODE helpNode;
       // helpNode.Set_copyall(true);
//        helpNode.Initialize();
        helpNode.Set_nodeID(i);
        helpNode.Set_numberOfParameter(data.Get_param(i));
        helpNode.Get_correspondingCPD().Set_correspondingNode(i);
        helpNode.Get_correspondingCPD().Set_numberOfParameter(data.Get_param(i));
        helpNode.UpdateCPD(data);
        helpNode.UpdateBIC(data.Get_NumberOfObservations());
        helpNode.Get_correspondingConstraint().Set_correspondingNode(i);

        helpNode.Get_correspondingConstraint().Set_cutoff(cutoff);

        helpNode.Get_correspondingConstraint().MIT(data);


        helpNode.Set_neighbors(helpNode.Get_correspondingConstraint().Get_dependentNode());
        _dependencies.push_back(helpNode.Get_correspondingConstraint().Get_dependentNode().size());
        _allNodes.push_back(helpNode);

        _score += helpNode.Get_scoreContribution();
        /*cout<<"node: "<<i<<":";
        for(int j=0;j<helpNode.Get_correspondingConstraint().Get_dependentNode().size();j++)
        {
            cout<<helpNode.Get_correspondingConstraint().Get_dependentNode(j)<<",";
        }
        cout<<endl;
        for(int j=0;j<helpNode.Get_correspondingConstraint().Get_dependentNode().size();j++)
        {
            cout<<"score: "<<helpNode.Get_correspondingConstraint().Get_orderedTestScore(j)<<",";
        }
        cout<<endl; */
    }
    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        _allNodes[i].CMIT(_allNodes,data,cutoff);
//        cout<<"node: "<<i<<":";
//        for(int j=0;j<_allNodes[i].Get_neighbors().size();j++)
//        {
//            cout<<_allNodes[i].Get_neighbors(j)<<",";
//        }
//        cout<<endl;
    }
    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        _allNodes[i].SymmetryCorrection(_allNodes);
    }
    for(int i = 0 ; i<data.Get_NumberOfNodes();i++)
    {
        _allNodes[i].SortNeighbor();
        //_numberOfNeighbors[i]=_allNodes[i].Get_neighbors().size();
        _numberOfNeighbors.push_back(_allNodes[i].Get_neighbors().size());
       // cout<<_numberOfNeighbors.size()<<endl;
        //cout<<_numberOfNeighbors[i]<<endl;
//        cout<<"i="<<i+1<<endl;
//        for(int j=0; j<_numberOfNeighbors[i];j++)
//        {
//            cout<<_allNodes[i].Get_neighbors(j)+1<<" ";
//        }
//        cout<<endl;
    }

}




void BN::Set_arc(int fromNode, int toNode, bool setparent)
{
    // update the edge set, and update numberOfEdge.
    EDGE temp_edge;
    temp_edge.Set_edge(fromNode,toNode);
    _edges.push_back(temp_edge);
    // update the child/parent/ancester/descendent
    _allNodes[fromNode].UpdateChild(toNode);
    _allNodes[fromNode].UpdateDescendant(toNode,_allNodes[toNode].Get_descendant());
    if(setparent)
    {
        _allNodes[toNode].UpdateParent(fromNode);
    }
    _allNodes[toNode].UpdateAncester(fromNode,_allNodes[fromNode].Get_ancester());
    // pass the ancester of the "fromNode" to the descendents of the "toNode"
    // also pass the descendents of the "toNode" to the ancester of the "fromNode"
    set<int>::iterator it;
    for(it=_allNodes[fromNode].Get_ancester().begin(); it!=_allNodes[fromNode].Get_ancester().end(); it++ )
    {
        _allNodes[(*it)].UpdateDescendant(toNode,_allNodes[toNode].Get_descendant());
    }
    for(it=_allNodes[toNode].Get_descendant().begin(); it!=_allNodes[toNode].Get_descendant().end(); it++ )
    {
        _allNodes[(*it)].UpdateAncester(fromNode,_allNodes[fromNode].Get_ancester());
    }
}



void BN::BICscore()
{
    for(unsigned int i=0; i< _allNodes.size(); i++)
    {
        _score += _allNodes[i].Get_scoreContribution();
    }
}

void BN::Update(DATA& data,double temper)
{
//    cout<<"tag: "<<_oneNeb<<endl;
//    cout<<"score: "<<_score<<endl;

    if(!(_oneNeb))
    {
        //cout<<1<<endl;
        vector<int> candidates;
        vector<int> triplet(3);
        triplet[0]=-1;
        triplet[1]=-1;
        triplet[2]=-1;
        candidates = _MinOverOne(_numberOfNeighbors);

        //cout<<"overone: "<<candidates.size()<<endl;
        //cout<<candidates.size()<<endl;
        if(candidates.empty())// no possible v-structures
        {
            _candOneNeb = _SearchNoneZero(_numberOfNeighbors);
            if(_candOneNeb.empty())
            {
                _growing = false;
            }
            else if( _candOneNeb.size()>1 )
            {
                int NodeA = -1 , NodeB = -1;
                int Ind;
                double totalScore = 0.0;
                for(int i = 0; i< _candOneNeb.size();i++)
                {
                    _totalScoreOneNeb += _allNodes[_candOneNeb[i]].Get_scoreOfNeighbors(0);
                    _scoreOneNeb.push_back(_allNodes[_candOneNeb[i]].Get_scoreOfNeighbors(0));
                }
                double randomNum = static_cast<double>( rand() ) / RAND_MAX;
                while(randomNum >= 1.0)
                {
                    randomNum = static_cast<double>( rand() ) / RAND_MAX;
                }
                randomNum *= _totalScoreOneNeb;
                double cumScore = 0.0;
                for(int i = 0; i<_candOneNeb.size(); i++)
                {
                    cumScore += _scoreOneNeb[i];
                    if(randomNum < cumScore)
                    {
                        NodeA = _candOneNeb[i];
                        Ind=i;
                        break;
                    }
                }

                if(NodeA==-1)
                {
                    NodeA = _candOneNeb[_candOneNeb.size()-1];
                    Ind = _candOneNeb.size()-1;
                }

                NodeB = _allNodes[NodeA].Get_neighbors(0);
                _candOneNeb.erase(_candOneNeb.begin()+Ind);// remove node A
                _scoreOneNeb.erase(_scoreOneNeb.begin()+Ind);
                int tmpInd = _IndexSearch(NodeB);// remove node B
                _candOneNeb.erase(_candOneNeb.begin()+tmpInd);
                _scoreOneNeb.erase(_scoreOneNeb.begin()+tmpInd);
//cout<<"pair"<<endl;
                _SelectDirection(NodeA,NodeB,data,temper);
                _oneNeb=true;
            }
            else
            {
                int NodeA,NodeB;
                NodeA = _candOneNeb[0];
                NodeB = _allNodes[NodeA].Get_neighbors(0);
                //_allNodes[NodeA].ClearRank(0,true);
//cout<<"pair"<<endl;
                _SelectDirection(NodeA,NodeB,data,temper);
                _growing = false;
            }
        }
        else
        {
            if(candidates.size()==1)
            {
                triplet[1] = candidates[0];
            }
            else
            {
                double arandomNum = static_cast<double>( rand() ) / RAND_MAX;
                while(arandomNum >= 1.0)
                {
                    arandomNum = static_cast<double>( rand() ) / RAND_MAX;
                }
                //cout<<arandomNum<<endl;
                triplet[1] = candidates[ min( static_cast<int>(arandomNum*candidates.size()), static_cast<int>(candidates.size())-1 ) ];
            }

            double randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }

            double totalScore = 0.0;
            for(int i=0; i< _allNodes[triplet[1]].Get_neighbors().size();i++)
            {
                totalScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
            }

            randomNum *= totalScore;

            //cout<<randomNum<<endl;
            double cumScore = 0.0;
            int tmpSampleInde = -1;
            for(int i=0; i<_allNodes[triplet[1]].Get_neighbors().size();i++)
            {
                cumScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
                if(randomNum < cumScore)
                {
                    tmpSampleInde = i;
                    break;
                }
            }
            if(tmpSampleInde == -1)
            {
                tmpSampleInde =  _allNodes[triplet[1]].Get_neighbors().size()-1;
            }
            triplet[0] = _allNodes[triplet[1]].Get_neighbors(tmpSampleInde);
            _allNodes[triplet[1]].ClearRank(tmpSampleInde,true);
            _numberOfNeighbors[triplet[1]]--;
            _allNodes[triplet[0]].ClearRank(triplet[1],false);
            _numberOfNeighbors[triplet[0]]--;
//cout<<"4neb: ";
//
//for(int i=0;i<_allNodes[4].Get_neighbors().size();i++)
//{
//    cout<<_allNodes[4].Get_neighbors(i)<<",";
//}
//cout<<endl;
            //see if these two nodes have common neighbors, where we should first consider that cases, and this kind of case refer to complete graph, K_3
//            for(int i=0; i<_allNodes[triplet[0]].Get_neighbors().size();i++)
 //           {
   //             cout<<_allNodes[triplet[0]].Get_neighbors(i)<<",";
     //       }
       //     cout<<endl;


         //   for(int i=0; i<_allNodes[triplet[1]].Get_neighbors().size();i++)
          //  {
         //       cout<<_allNodes[triplet[1]].Get_neighbors(i)<<",";
         //   }
//cout<<endl;
            vector<int> tmpCommonNeb(_allNodes[triplet[0]].Get_neighbors().size()+_allNodes[triplet[1]].Get_neighbors().size());
            vector<int> :: iterator it;
            it = set_intersection(_allNodes[triplet[0]].Get_neighbors().begin(),_allNodes[triplet[0]].Get_neighbors().end(),
                                _allNodes[triplet[1]].Get_neighbors().begin(),_allNodes[triplet[1]].Get_neighbors().end(),
                                tmpCommonNeb.begin());
            tmpCommonNeb.resize(it-tmpCommonNeb.begin());

            randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }
            //determine if the triplet forms a complete graph K_3
            if(tmpCommonNeb.empty())//incomplete graph
            {
 //cout<<"incomp"<<endl;

                double anotherTotalScore = 0.0;
                for(int i=0; i< _allNodes[triplet[1]].Get_neighbors().size();i++)
                {
                    anotherTotalScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
                }
                randomNum *= anotherTotalScore;

                cumScore = 0.0;
                tmpSampleInde = -1;
                for(int i=0; i<_allNodes[triplet[1]].Get_neighbors().size();i++)
                {
                    cumScore += _allNodes[triplet[1]].Get_scoreOfNeighbors(i);
                    if(randomNum < cumScore)
                    {
                        tmpSampleInde = i;
                        break;
                    }
                }
                if(tmpSampleInde == -1)
                {
                    tmpSampleInde =  _allNodes[triplet[1]].Get_neighbors().size()-1;
                }
                triplet[2] = _allNodes[triplet[1]].Get_neighbors(tmpSampleInde);
                _allNodes[triplet[1]].ClearRank(tmpSampleInde,true);
                _numberOfNeighbors[triplet[1]]--;
                _allNodes[triplet[2]].ClearRank(triplet[1],false);
                _numberOfNeighbors[triplet[2]]--;
//cout<<"incomp,"<<triplet[0]<<","<<triplet[1]<<","<<triplet[2]<<endl;
//cout<<"triple"<<endl;
                _TripletSelection(triplet,data,temper,false);
            }
            else//complete graph
            {
// cout<<"compl"<<endl;
                triplet[2] = tmpCommonNeb[min(static_cast<int> (randomNum * tmpCommonNeb.size()), static_cast<int>(tmpCommonNeb.size())-1)];

                _allNodes[triplet[0]].ClearRank(triplet[2],false);

                _numberOfNeighbors[triplet[0]]--;
                _allNodes[triplet[1]].ClearRank(triplet[2],false);

                _numberOfNeighbors[triplet[1]]--;
                _allNodes[triplet[2]].ClearRank(triplet[0],false);

                _numberOfNeighbors[triplet[2]]--;
                _allNodes[triplet[2]].ClearRank(triplet[1],false);
//cout<<"4neb4 ";
//
//for(int i=0;i<_allNodes[4].Get_neighbors().size();i++)
//{
//    cout<<_allNodes[4].Get_neighbors(i)<<",";
//}
//cout<<endl;
                _numberOfNeighbors[triplet[2]]--;
//cout<<"triple"<<endl;
                _TripletSelection(triplet,data,temper,true);
//cout<<"comp,"<<triplet[0]<<","<<triplet[1]<<","<<triplet[2]<<endl;
//cout<<"4nebaftertrisel: ";
//
//for(int i=0;i<_allNodes[4].Get_neighbors().size();i++)
//{
//    cout<<_allNodes[4].Get_neighbors(i)<<",";
//}
//cout<<endl;
            }


        }
    }
    else// only nodes with one neighbor remain
    {
//cout<<"haha"<<endl;
        if(_candOneNeb.empty())
        {
            _growing = false;
        }
        else if( _candOneNeb.size()>2 )
        {
            int NodeA = -1 , NodeB = -1;
            int Ind;
            double totalScore = 0.0;
            for(int i = 0; i< _candOneNeb.size();i++)
            {
                _totalScoreOneNeb += _allNodes[_candOneNeb[i]].Get_scoreOfNeighbors(0);
            }
            double randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }
            randomNum *= _totalScoreOneNeb;
            double cumScore = 0.0;
            for(int i = 0; i<_candOneNeb.size(); i++)
            {
                cumScore += _scoreOneNeb[i];
                if(randomNum < cumScore)
                {
                    NodeA = _candOneNeb[i];
                    Ind=i;
                    break;
                }
            }

            if(NodeA==-1)
            {
                NodeA = _candOneNeb[_candOneNeb.size()-1];
                Ind = _candOneNeb.size()-1;
            }

            NodeB = _allNodes[NodeA].Get_neighbors(0);
            _candOneNeb.erase(_candOneNeb.begin()+Ind);// remove node A
            _scoreOneNeb.erase(_scoreOneNeb.begin()+Ind);
            int tmpInd = _IndexSearch(NodeB);// remove node B
            _candOneNeb.erase(_candOneNeb.begin()+tmpInd);
            _scoreOneNeb.erase(_scoreOneNeb.begin()+tmpInd);
//cout<<"pair"<<endl;
            _SelectDirection(NodeA,NodeB,data,temper);
        }
        else
        {
            int NodeA,NodeB;
            NodeA = _candOneNeb[0];
            NodeB = _allNodes[NodeA].Get_neighbors(0);
            //_allNodes[NodeA].ClearRank(0,true);
//cout<<"pair"<<endl;
            _SelectDirection(NodeA,NodeB,data,temper);
            _growing = false;
        }










/*

        if(_candOneNeb.empty())
        {
            _growing = false;
        }
        else if( _candOneNeb.size()>1 )
        {
            int NodeA,NodeB;
            double randomNum = static_cast<double>( rand() ) / RAND_MAX;
            while(randomNum >= 1.0)
            {
                randomNum = static_cast<double>( rand() ) / RAND_MAX;
            }
            randomNum *= _totalScoreOneNeb;
            double cumScore = 0.0;
            for(int i = 0; i<_candOneNeb.size(); i++)
            {
                cumScore += _scoreOneNeb[i];
                if(randomNum < cumScore)
                {
                    NodeA = _candOneNeb[i];
                    NodeB = _allNodes[NodeA].Get_neighbors(0);
                    _candOneNeb.erase(_candOneNeb.begin()+i);// remove node A
                    _totalScoreOneNeb -= _scoreOneNeb[i];
                    _scoreOneNeb.erase(_scoreOneNeb.begin()+i);
                    int tmpInd = _IndexSearch(NodeB);// remove node B
                    _candOneNeb.erase(_candOneNeb.begin()+tmpInd);
                    _totalScoreOneNeb -= _scoreOneNeb[tmpInd];
                    _scoreOneNeb.erase(_scoreOneNeb.begin()+tmpInd);
                    _SelectDirection(NodeA,NodeB,data,temper);
                    break;
                }
            }
          //  _oneNeb=true;
        }
        else
        {
            int NodeA,NodeB;
            NodeA = _candOneNeb[0];
            NodeB = _allNodes[NodeA].Get_neighbors(0);
            _SelectDirection(NodeA,NodeB,data,temper);
            _growing = false;

        }
*/

    }
//cout<<"-------------------------"<<endl;


}




void BN::HC(DATA& data,int depth)
{
//first generate a sequence of updating order, algorithm: Knuth shuffles
    vector<int> sequence;
    int No_nodes = data.Get_NumberOfNodes();
	for(int i=0;i<No_nodes;i++)
	{
	 sequence.push_back(i);
	}

	for(int i=0;i<No_nodes;i++)
	{
		int j = rand()%(i+1);
		sequence[i]=sequence[j];
		sequence[j]=i;

	}

//	int revcounter = 0 ;
//	int delcounter = 0 ;
//------------------------------------------------------------------------//
    int calibrate = 0;
    vector<bool> changedAcc;
    vector<int> changedNode;
    changedAcc.reserve(No_nodes);
    changedNode.reserve(No_nodes);
    for(int i=0;i<No_nodes;i++)
    {
        changedAcc.push_back(true);
        changedNode.push_back(0);
    }
    int round=0;

    if(depth>0)
    {

        while(calibrate == 0 && round < depth )
        {
            round++;
            calibrate = 1;
            for(int i=0; i<No_nodes;i++)
            {

                if(changedAcc[i])
                {
                    for(int j=0;j<No_nodes;j++) //sequence
                    {

                        if(i!= sequence[j])
                        {


                            if(!(_ExistArc(sequence[j],i)) && !(_ExistArc(i,sequence[j])) )//no arc
                            {

                                if(_IsDAG(sequence[j],i))
                                {

                                    NODE tmpNode = _allNodes[i];
                                    tmpNode.UpdateParent(sequence[j]);
                                    tmpNode.UpdateCPD(data);
                                    tmpNode.UpdateBIC(data.Get_NumberOfObservations());
                                    if(tmpNode.Get_scoreContribution()>_allNodes[i].Get_scoreContribution())
                                    {
                                        _score -= _allNodes[i].Get_scoreContribution();
                                        _allNodes[i] = tmpNode;
                                        Set_arc(sequence[j],i,false);
                                        _score += tmpNode.Get_scoreContribution();
                                        calibrate*=0;
                                        changedNode[i]++;
                                    }

                                }

                            }
                           /* else if(_ExistArc(sequence[j],i))
                            {
                                double scoreDifDel = 0.0, scoreDifRev = 0.0;
                                // this part calculates the score change after deleting
                                NODE delNode = _allNodes[i];
                                delNode.RemoveParent(sequence[j]);
                                delNode.UpdateCPD(data);
                                delNode.UpdateBIC(data.Get_NumberOfObservations());
                                scoreDifDel = delNode.Get_scoreContribution()-_allNodes[i].Get_scoreContribution();

                                // this part calculates the score change after reversing
                                NODE fromNode =_allNodes[sequence[j]] ,toNode=_allNodes[i];
                                fromNode.UpdateParent(i);
                                toNode.RemoveParent(sequence[j]);
                                fromNode.UpdateCPD(data);
                                toNode.UpdateCPD(data);
                                fromNode.UpdateBIC(data.Get_NumberOfObservations());
                                toNode.UpdateBIC(data.Get_NumberOfObservations());
                                scoreDifRev = fromNode.Get_scoreContribution()+toNode.Get_scoreContribution()-_allNodes[i].Get_scoreContribution()-_allNodes[sequence[j]].Get_scoreContribution();

                                if(scoreDifRev>0 && scoreDifRev>scoreDifDel)
                                {
                                    if(Rev_arc(sequence[j],i))
                                    {
                                        _allNodes[i].Set_correspondingCPD(toNode.Get_correspondingCPD());
                                        _allNodes[i].Set_scoreContribution(toNode.Get_scoreContribution());
                                        _allNodes[sequence[j]].Set_correspondingCPD(fromNode.Get_correspondingCPD());
                                        _allNodes[sequence[j]].Set_scoreContribution(fromNode.Get_scoreContribution());
                                        calibrate*=0;
                                        changedNode[i]++;


                                    }
                                }
                                else if( scoreDifDel>0 && scoreDifDel>scoreDifRev)
                                {
                                    Del_arc(sequence[j],i);
                                    _allNodes[i].Set_correspondingCPD(delNode.Get_correspondingCPD());
                                    _allNodes[i].Set_scoreContribution(delNode.Get_scoreContribution());
                                    calibrate*=0;
                                    changedNode[i]++;

                                }


                            }*/
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

                    if(changedNode[i]==0)
                    {
                        changedAcc[i] = false;
                        changedNode[i] = 0 ;
                    }
                    else
                    {
                        changedAcc[i] = true;
                        changedNode[i] = 0;
                    }
                }
            }
        }
    }
    else if(depth == -1 )
    {
        while(calibrate == 0)
        {
            calibrate = 1;
            for(int i=0; i<No_nodes;i++)
            {
                if(changedAcc[i])
                {
                    for(int j=0;j<No_nodes;j++) //sequence
                    {
                        if(i!= sequence[j])
                        {
                            if(!(_ExistArc(sequence[j],i)) && !(_ExistArc(i,sequence[j])) )//no arc
                            {
                                if(_IsDAG(sequence[j],i))
                                {
                                    NODE tmpNode = _allNodes[i];
                                    tmpNode.UpdateParent(sequence[j]);
                                    tmpNode.UpdateCPD(data);
                                    tmpNode.UpdateBIC(data.Get_NumberOfObservations());
                                    if(tmpNode.Get_scoreContribution()>_allNodes[i].Get_scoreContribution())
                                    {
                                        _score -= _allNodes[i].Get_scoreContribution();
                                        _allNodes[i] = tmpNode;
                                        Set_arc(sequence[j],i,false);
                                        _score += tmpNode.Get_scoreContribution();
                                        calibrate*=0;
                                        changedNode[i]++;

                                    }

                                }

                            }
                           /* else if(_ExistArc(sequence[j],i))
                            {
                                double scoreDifDel = 0.0, scoreDifRev = 0.0;
                                // this part calculates the score change after deleting
                                NODE delNode = _allNodes[i];
                                delNode.RemoveParent(sequence[j]);
                                delNode.UpdateCPD(data);
                                delNode.UpdateBIC(data.Get_NumberOfObservations());
                                scoreDifDel = delNode.Get_scoreContribution()-_allNodes[i].Get_scoreContribution();

                                // this part calculates the score change after reversing
                                NODE fromNode =_allNodes[sequence[j]] ,toNode=_allNodes[i];
                                fromNode.UpdateParent(i);
                                toNode.RemoveParent(sequence[j]);
                                fromNode.UpdateCPD(data);
                                toNode.UpdateCPD(data);
                                fromNode.UpdateBIC(data.Get_NumberOfObservations());
                                toNode.UpdateBIC(data.Get_NumberOfObservations());
                                scoreDifRev = fromNode.Get_scoreContribution()+toNode.Get_scoreContribution()-_allNodes[i].Get_scoreContribution()-_allNodes[sequence[j]].Get_scoreContribution();

                                if(scoreDifRev>0 && scoreDifRev>scoreDifDel)
                                {
                                    if(Rev_arc(sequence[j],i))
                                    {
                                        _allNodes[i].Set_correspondingCPD(toNode.Get_correspondingCPD());
                                        _allNodes[i].Set_scoreContribution(toNode.Get_scoreContribution());
                                        _allNodes[sequence[j]].Set_correspondingCPD(fromNode.Get_correspondingCPD());
                                        _allNodes[sequence[j]].Set_scoreContribution(fromNode.Get_scoreContribution());
                                        calibrate*=0;
                                        changedNode[i]++;
                                        revcounter++;
                                    }
                                }
                                else if( scoreDifDel>0 && scoreDifDel>scoreDifRev)
                                {
                                    Del_arc(sequence[j],i);
                                    _allNodes[i].Set_correspondingCPD(delNode.Get_correspondingCPD());
                                    _allNodes[i].Set_scoreContribution(delNode.Get_scoreContribution());
                                    calibrate*=0;
                                    changedNode[i]++;
                                    delcounter++;
                                }


                            }*/
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

                    if(changedNode[i]==0)
                    {
                        changedAcc[i] = false;
                        changedNode[i] = 0 ;
                    }
                    else
                    {
                        changedAcc[i] = true;
                        changedNode[i] = 0;
                    }
                }
            }
        }
    }

//cout<<"revcount: "<<revcounter<<", delcount: "<<delcounter<<endl;

}
            /*
            for(int j=0;j<No_nodes;j++) //sequence
            {
                if(i!= sequence[j] && !(_ExistArc(i,sequence[j])))
                {
                    if(_IsDAG(i,sequence[j]) && !(_IsDAG(sequence[j],i))) //i-> sequence[j]
                    {
                       NODE tmpNode = _allNodes[sequence[j]];
                       tmpNode.UpdateParent(i);
                       tmpNode.UpdateCPD(data);
                       tmpNode.UpdateBIC(data.Get_NumberOfObservations());
                       if(tmpNode.Get_scoreContribution()>_allNodes[sequence[j]].Get_scoreContribution())
                        {
                            _score -= _allNodes[sequence[j]].Get_scoreContribution();
                            _allNodes[sequence[j]] = tmpNode;
                            Set_arc(i,sequence[j],false);
                            _score += tmpNode.Get_scoreContribution();
                            calibrate*=0;
                        }
                    }
                    else if(!(_IsDAG(i,sequence[j])) && (_IsDAG(sequence[j],i)))// sequence[j]->i
                    {
                        NODE tmpNode = _allNodes[i];
                       tmpNode.UpdateParent(sequence[j]);
                       tmpNode.UpdateCPD(data);
                       tmpNode.UpdateBIC(data.Get_NumberOfObservations());
                       if(tmpNode.Get_scoreContribution()>_allNodes[i].Get_scoreContribution())
                        {
                            _score -= _allNodes[i].Get_scoreContribution();
                            _allNodes[i] = tmpNode;
                            Set_arc(sequence[j],i,false);
                            _score += tmpNode.Get_scoreContribution();
                            calibrate*=0;
                        }
                    }
                    else if(_IsDAG(i,sequence[j]) && _IsDAG(sequence[j],i))//both
                    {
                        NODE tmpNodeA = _allNodes[sequence[j]];
                       tmpNodeA.UpdateParent(i);
                       tmpNodeA.UpdateCPD(data);
                       tmpNodeA.UpdateBIC(data.Get_NumberOfObservations());

                       NODE tmpNodeB = _allNodes[i];
                       tmpNodeB.UpdateParent(sequence[j]);
                       tmpNodeB.UpdateCPD(data);
                       tmpNodeB.UpdateBIC(data.Get_NumberOfObservations());

                       if(tmpNodeA.Get_scoreContribution()>_allNodes[sequence[j]].Get_scoreContribution() || tmpNodeB.Get_scoreContribution()>_allNodes[i].Get_scoreContribution() )
                        {
                            if(tmpNodeA.Get_scoreContribution()>tmpNodeB.Get_scoreContribution())
                            {
                                _score -= _allNodes[sequence[j]].Get_scoreContribution();
                                _allNodes[sequence[j]] = tmpNodeA;
                                Set_arc(i,sequence[j],false);
                                _score += tmpNodeA.Get_scoreContribution();
                                calibrate*=0;
                            }
                            else
                            {
                                _score -= _allNodes[i].Get_scoreContribution();
                                _allNodes[i] = tmpNodeB;
                                Set_arc(sequence[j],i,false);
                                _score += tmpNodeB.Get_scoreContribution();
                                calibrate*=0;
                            }

                        }


                    }
                    else
                    {
                        //dummy
                    }
                }
                else
                {
                    continue;
                }
            }


        }
*/


/*void BN::Del_arc(int fromNode,int toNode)
{
    if(_ExistArc(fromNode,toNode))
    {
        _allNodes[fromNode].RemoveChild(toNode);
        _allNodes[toNode].RemoveParent(fromNode);
        _RecombAncester(toNode);
        _RecombDescedant(fromNode);
        for(int i =0;i<_edges.size();i++)
        {
            if(_edges[i].Matched(fromNode,toNode))
            {
                _edges.erase(_edges.begin()+i);

            }
        }
    }


}
*/

/*bool BN::Rev_arc(int fromNode,int toNode)// this from and to are the present edge
{
    Del_arc(fromNode,toNode);
    if(_IsDAG(toNode,fromNode))
    {
        Set_arc(toNode,fromNode,true);
        return true;
    }
    else
    {
        Set_arc(fromNode,toNode,true);
        return false;
    }


}
*/







void BN::_SelectDirection(int NodeA, int NodeB, DATA& data,double temper)
{
    //_numberOfNeighbors[NodeA]--;
    //_numberOfNeighbors[NodeB]--;

    //cout<<NodeA<<" "<<NodeB<<endl;
    //cout<<"size "<<_allNodes[NodeB].Get_neighbors().size()<<endl;
    //for(int p = 0; p<_allNodes[NodeB].Get_neighbors().size();p++)
    //{
    //    cout<<_allNodes[NodeB].Get_neighbors(p)<<",";
   // }
    //cout<<endl;
    //_allNodes[NodeB].ClearRank(NodeA,false);

    bool flagAtoB, flagBtoA;
    flagAtoB = _IsDAG(NodeA,NodeB);
    flagBtoA = _IsDAG(NodeB,NodeA);


    /* some cases:
    1) only A to B valid,
    2) only B to A valid,
    3) both are valid,
    4) neither is valid.
    */
    if(flagAtoB && (!flagBtoA) ) //only A to B valid,
    {
        double tempScore = _score;
        NODE tempNode = _allNodes[NodeB];
        tempScore -= tempNode.Get_scoreContribution();
        tempNode.UpdateParent(NodeA);
        tempNode.UpdateCPD(data);
        tempNode.UpdateBIC(data.Get_NumberOfObservations());
        tempScore += tempNode.Get_scoreContribution();
        double probSampleAtoB;
        probSampleAtoB = exp((tempScore-_score)/(temper))/(1.0+exp((tempScore-_score)/(temper)));
        double aRandomNumber = static_cast<double>(rand())/RAND_MAX;
        while(aRandomNumber >= 1.0)
        {
            aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
        }
        if(aRandomNumber <= probSampleAtoB )// which means we do accept A to B
        {
//cout<<NodeA<<"->"<<NodeB<<endl;
            _allNodes[NodeB] = tempNode;
            Set_arc(NodeA,NodeB,false);
            _score = tempScore;
        }

    }
    else if( (!flagAtoB) && flagBtoA ) //only B to A valid,
    {
        double tempScore = _score;
        NODE tempNode = _allNodes[NodeA];
        tempScore -= tempNode.Get_scoreContribution();
        tempNode.UpdateParent(NodeB);
        tempNode.UpdateCPD(data);
        tempNode.UpdateBIC(data.Get_NumberOfObservations());
        tempScore += tempNode.Get_scoreContribution();
        double probSampleBtoA;
        probSampleBtoA = exp((tempScore-_score)/(temper))/(1.0+exp((tempScore-_score)/(temper)));
        double aRandomNumber = static_cast<double>(rand())/RAND_MAX;
        while(aRandomNumber >= 1.0)
        {
            aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
        }
        if(aRandomNumber <= probSampleBtoA )// which means we do accept B to A
        {
//cout<<NodeB<<"->"<<NodeA<<endl;
            _allNodes[NodeA] = tempNode;
            Set_arc(NodeB,NodeA,false);
            _score = tempScore;
        }
    }
    else if( flagAtoB && flagBtoA ) //both are valid,
    {
        //A to B part
        double tempScoreAtoB = _score;
        NODE tempNodeAtoB = _allNodes[NodeB];
        tempScoreAtoB -= tempNodeAtoB.Get_scoreContribution();
        tempNodeAtoB.UpdateParent(NodeA);
        tempNodeAtoB.UpdateCPD(data);
        tempNodeAtoB.UpdateBIC(data.Get_NumberOfObservations());
        tempScoreAtoB += tempNodeAtoB.Get_scoreContribution();
        //B to A part
        double tempScoreBtoA = _score;
        NODE tempNodeBtoA = _allNodes[NodeA];
        tempScoreBtoA -= tempNodeBtoA.Get_scoreContribution();
        tempNodeBtoA.UpdateParent(NodeB);
        tempNodeBtoA.UpdateCPD(data);
        tempNodeBtoA.UpdateBIC(data.Get_NumberOfObservations());
        tempScoreBtoA += tempNodeBtoA.Get_scoreContribution();
        //probabilities
        double probSampleAtoB,probSampleBtoA;
        probSampleAtoB = exp((tempScoreAtoB-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
        probSampleBtoA = exp((tempScoreBtoA-_score)/(temper))/(1.0+exp((tempScoreAtoB-_score)/(temper))+exp((tempScoreBtoA-_score)/(temper)));
        //sample
        double aRandomNumber = static_cast<double>(rand())/RAND_MAX;
        while(aRandomNumber >= 1.0)
        {
            aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
        }
        if(aRandomNumber <= probSampleAtoB )// which means we do accept A to B
        {
//cout<<NodeA<<"->"<<NodeB<<endl;
            _allNodes[NodeB] = tempNodeAtoB;
            Set_arc(NodeA,NodeB,false);
            _score = tempScoreAtoB;
        }
        else if(aRandomNumber <= probSampleAtoB+probSampleBtoA) // accept B to A
        {
//cout<<NodeB<<"->"<<NodeA<<endl;
            _allNodes[NodeA] = tempNodeBtoA;
            Set_arc(NodeB,NodeA,false);
            _score = tempScoreBtoA;
        }
    }
    else
    {
        //dummy...
    }
   // cout<<_score<<endl;

}

void BN::_TripletSelection(vector<int>& candidateNodes, DATA& data, double temper,bool isCompleteGraph)
{
    double originalBIC = 0.0;
    vector<NODE> originalNodes(3);
    double TotalProb = 0.0;
    vector<bool> tmpAddedEdge(3);
    int MaxInd = 3;
    vector< vector<NODE> > candidates;
    vector<double> potentialBIC;
    vector<double> probOfSampling;
    vector<vector<int> > configurations;
    vector<int> tmpConfigurations(3);
    vector< vector <bool> > AddedEdge;
    for(int i=0; i<3;i++)
    {
//cout<<"original "<<candidateNodes[i]<<","<<_allNodes[candidateNodes[i]].Get_scoreContribution()<<endl;
        originalBIC += _allNodes[candidateNodes[i]].Get_scoreContribution();
        originalNodes[i] = _allNodes[candidateNodes[i]];
        tmpAddedEdge[i] = true;

    }
//cout<<"oBIC: "<<originalBIC<<endl;
    if(isCompleteGraph)
    {
        MaxInd = 3;
        candidates.reserve(27);
        potentialBIC.reserve(27);
        probOfSampling.reserve(27);
        configurations.reserve(27);

        AddedEdge.reserve(27);
        for(int i =0; i< 27;i++)
        {
            tmpConfigurations[0] = i%3;
            tmpConfigurations[1] = (i/3)%3;
            tmpConfigurations[2] = (i/9)%3;
            configurations.push_back(tmpConfigurations);
            candidates.push_back(originalNodes);
            potentialBIC.push_back(0.0);
            probOfSampling.push_back(0.0);
            AddedEdge.push_back(tmpAddedEdge);
        }

    }
    else
    {
        MaxInd = 2;
        candidates.reserve(9);
        potentialBIC.reserve(9);
        probOfSampling.reserve(9);
        configurations.reserve(9);
       // tmpConfigurations.reserve(2);
        AddedEdge.reserve(9);
        for(int i =0; i< 9;i++)
        {
            tmpConfigurations[0] = i%3;
            tmpConfigurations[1] = (i/3)%3;
            configurations.push_back(tmpConfigurations);
            candidates.push_back(originalNodes);
            potentialBIC.push_back(0.0);
            probOfSampling.push_back(0.0);
            AddedEdge.push_back(tmpAddedEdge);
        }
    }
    probOfSampling[0] = 1.0;//exp(0)
    potentialBIC[0] = originalBIC;
    TotalProb += probOfSampling[0];
    for(int i =1 ; i< pow(3,MaxInd);i++)
    {
        for(int j=0;j<MaxInd;j++)
        {
            switch( configurations[i][j])
            {
                case 0:
                  //  AddedEdge[i][j] = true;
                    break;
                case 1: //j -> j+1, fromNode= j%3, and toNode= (j+1)%3
                    if( candidates[i][j%3].Get_ancester().find(candidateNodes[(j+1)%3]) ==  candidates[i][j%3].Get_ancester().end() ) //if it is DAG
                    {

                        //update fromNode: descendent/ child, and if (j+2)%3 in ancester of fromNode add toNode to its descendents as well
                        candidates[i][j%3].UpdateChild(candidateNodes[(j+1)%3]);
                        candidates[i][j%3].UpdateDescendant(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_descendant() );
                        if(candidates[i][j%3].Get_ancester().find(candidateNodes[(j+2)%3]) !=  candidates[i][j%3].Get_ancester().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateDescendant(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_descendant() );
                        }
                        // update toNode
                        candidates[i][(j+1)%3].UpdateParent(candidateNodes[(j)%3]);
                        candidates[i][(j+1)%3].UpdateAncester(candidateNodes[(j)%3],candidates[i][(j)%3].Get_ancester() );
                        if(candidates[i][(j+1)%3].Get_descendant().find(candidateNodes[(j+2)%3]) !=  candidates[i][(j+1)%3].Get_descendant().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateAncester(candidateNodes[(j)%3],candidates[i][(j)%3].Get_ancester() );
                        }

                    }
                    else
                    {
                        AddedEdge[i][j] = false;
                    }
                    break;
                case 2: // j<- j+1, fromNode= (j+1)%3, and toNode= j%3
                    if( candidates[i][(j+1)%3].Get_ancester().find(candidateNodes[(j)%3]) ==  candidates[i][(j+1)%3].Get_ancester().end() ) //if it is DAG
                    {
                        //AddedEdge[i][j] = true;
                        //update fromNode: descendent/ child, and if (j+2)%3 in ancester of fromNode add toNode to its descendents as well
                        //candidates[i][(j+1)%3].UpdateChild(candidateNodes[(j)%3]);
                        candidates[i][(j+1)%3].UpdateDescendant(candidateNodes[(j)%3],candidates[i][(j)%3].Get_descendant() );
                        if(candidates[i][(j+1)%3].Get_ancester().find(candidateNodes[(j+2)%3]) !=  candidates[i][(j+1)%3].Get_ancester().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateDescendant(candidateNodes[(j)%3],candidates[i][(j)%3].Get_descendant() );
                        }
                        // update toNode
                        candidates[i][(j)%3].UpdateParent(candidateNodes[(j+1)%3]);
                        candidates[i][(j)%3].UpdateAncester(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_ancester() );
                        if(candidates[i][(j)%3].Get_descendant().find(candidateNodes[(j+2)%3]) !=  candidates[i][j%3].Get_descendant().end()) //(j+2)%3 is an ancester of j%3
                        {
                            candidates[i][(j+2)%3].UpdateAncester(candidateNodes[(j+1)%3],candidates[i][(j+1)%3].Get_ancester() );
                        }

                    }
                    else
                    {
                        AddedEdge[i][j] =false;
                    }
                    break;
                default :
                cout<<"invalid type of connections, and you should not see this message unless there is something wrong."<<endl;
                break;
            }
//cout<<"add or not:"<<AddedEdge[i][0]<<","<<AddedEdge[i][1]<<","<<AddedEdge[i][2]<<endl;
        }

        if( (AddedEdge[i][0])&&(AddedEdge[i][1])&&(AddedEdge[i][2]) ) // If any edge added to current configuration then update all the candidate Nodes.
        {
            for(int j=0;j<3;j++)
            {

                candidates[i][j].UpdateCPD(data);
                candidates[i][j].UpdateBIC(data.Get_NumberOfObservations());
                potentialBIC[i] += candidates[i][j].Get_scoreContribution();

            }
//cout<<"pBIC: "<<i<<"th "<<potentialBIC[i]<<endl;
            probOfSampling[i] += exp( (potentialBIC[i]-originalBIC)/temper );
        }
        else
        {
            probOfSampling[i] = 0.0;
        }
        TotalProb += probOfSampling[i];
    }

    for(int i=0; i<pow(3,MaxInd);i++)
    {
        probOfSampling[i] /= TotalProb;
    }
//Sample a random number
    double aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
    while(aRandomNumber > 1.0)
    {
        aRandomNumber = static_cast<double>( rand() ) / RAND_MAX;
    }
    double cumulativeProb = 0.0;
    int SampledEdge = -1;
    for(int i=0;i<pow(3,MaxInd);i++)
    {
        cumulativeProb += probOfSampling[i];
        if(aRandomNumber < cumulativeProb)
        {
            SampledEdge = i;
            break;
        }
    }

    // this part is to avoid numerical issue.
    if(SampledEdge == -1)
    {
        for(int i=probOfSampling.size()-1;i>=0;i--)
        {
            if(probOfSampling[i]!=0.0)
            {
                SampledEdge = pow(3,MaxInd)-1;
            }
        }

    }

    _score -= originalBIC;

    _score += potentialBIC[SampledEdge];

    //copy the nodes from temp to allnodes
    for(int j=0; j<3 ; j++)
    {
       // candidates[SampledEdge][j].Set_copyall(false);
        _allNodes[candidateNodes[j]].UpdateResult(candidates[SampledEdge][j]);
    }



    for(int j=0;j<MaxInd;j++)
        {
            switch( configurations[SampledEdge][j])
            {
                case 0:

                    break;
                case 1: //j -> j+1, fromNode= j%3, and toNode= (j+1)%3
                    if( _IsDAG(candidateNodes[j%3],candidateNodes[(j+1)%3]) ) //if it is DAG
                    {
                        Set_arc(candidateNodes[j%3],candidateNodes[(j+1)%3],true);
                        //cout<<candidateNodes[j%3]+1<<"->"<<candidateNodes[(j+1)%3]+1<<endl;
                    }
                    break;
                case 2: // j<- j+1, fromNode= (j+1)%3, and toNode= j%3
                    if(_IsDAG(candidateNodes[(j+1)%3],candidateNodes[(j)%3])) //if it is DAG
                    {
                        Set_arc(candidateNodes[(j+1)%3],candidateNodes[(j)%3],true);
                        //cout<<candidateNodes[(j+1)%3]+1<<"->"<<candidateNodes[j%3]+1<<endl;
                    }
                    break;
                default :
                    cout<<"invalid type of connections, and you should not see this message unless there is something wrong."<<endl;
                    break;
            }
        }

}



bool BN::_IsDAG(int fromNode, int toNode)
{
    // The method here is to see whether the "toNode" is an element of the ancester set of the "fromNode".
    bool DAG = true;
    if ( _allNodes[fromNode].Get_ancester().find( toNode ) != _allNodes[fromNode].Get_ancester().end() ) // if the iterator is different from the end means the "toNode" is found
        //from the ancester set of the "fromNode", then it is not a DAG.
    {
        DAG = false;
    }

    return DAG;
}


vector<int> BN::_SearchNoneZero(vector<int>& target)
{
    vector<int> index;
    for(unsigned int j=0; j<target.size();j++)
    {
        if(target[j] != 0)
        {
            index.push_back(j);
        }
    }
    return index;
}



vector< int > BN::_MinOverOne(vector<int>& target)// The algorithm is linear search.
{

    vector< int > index;
    int minimum;
    int ind;
    // the first loop is intended to detected if there is an at-least-two entry and start the search from the first at-least-two entry.
    for(unsigned int j=0; j<target.size();j++)
    {

        if(target[j]>1)
        {
            index.push_back(j);
            minimum = target[j];
            ind = j;
            break;
        }
        else
        {
            continue;
        }

    }
    // start the linear search.
    if(!index.empty())
    {
        for( unsigned int i=ind+1; i<target.size();i++)
        {
            if(target[i] == minimum)
            {
                index.push_back(i);
            }
            else if(target[i] < minimum && target[i]>1)
            {
                index.clear();// in case a smaller value came, clear the exisiting vector.
                index.push_back(i);
                minimum = target[i];
            }
            else
            {
                continue;
            }
        }
    }
    return index;

}


double BN::_Summation(vector<double>& target)// linear addition..
{
    double result = 0.0;

    for(unsigned int i=0; i< target.size();i++)
    {
        result += target[i];
    }
    return result;
}

int BN::_IndexSearch(int target)
{
    int L = 0;
    int U = _candOneNeb.size()-1;
    int tar = 0;
    while( (_candOneNeb[L]!=target) && (_candOneNeb[U]!=target) )
    {
        if(target<_candOneNeb[(L+U)/2])
        {
            U = (L+U)/2;
        }
        else
        {
            L = (L+U)/2;
        }
    }
    if (_candOneNeb[L]==target)
    {
        tar = L;
    }
    else if(_candOneNeb[U]==target)
    {
        tar = U;
    }
    else
    {
        cerr<<"Index search failed. It should not happen really."<<endl;
    }
    return tar;

}

bool BN::_ExistArc(int fromNode,int toNode)
{
    bool edge = false;

    if(_allNodes[fromNode].Get_child().find(toNode)!=_allNodes[fromNode].Get_child().end())
        //|| _allNodes[toNode].Get_child().find(fromNode)!=_allNodes[toNode].Get_child().end())
    {
        edge=true;
    }
    return edge;
}

/*void BN::_RecombAncester(int toNode)
{
    set<int> children_new,children_old;
    children_new.insert(toNode);
    children_old.insert(toNode);
    set<int>::iterator toNodes;
    while(!(children_new.empty()))
    {
        children_old = children_new;
        children_new.clear();
        for(toNodes = children_old.begin();toNodes != children_old.end();toNodes++)
        {
            if(!(_allNodes[*toNodes].Get_parent().empty()))
            {
                vector<int> tmpParent = _allNodes[*toNodes].Get_parent();
                set<int> tmpAncester;
                tmpAncester.insert(tmpParent.begin(),tmpParent.end());
                for(int j=0;j<tmpParent.size();j++)
                {
                    tmpAncester.insert(_allNodes[ tmpParent[j] ].Get_ancester().begin(),_allNodes[ tmpParent[j] ].Get_ancester().end());
                }
                _allNodes[*toNodes].Set_ancester(tmpAncester);
            }
            if(!(_allNodes[*toNodes].Get_child().empty()))
            {
                children_new.insert(_allNodes[*toNodes].Get_child().begin(),_allNodes[*toNodes].Get_child().end());
            }
        }
    }
}

void BN::_RecombDescedant(int fromNode)
{
    set<int> parent_new,parent_old;
    parent_new.insert(fromNode);
    parent_old.insert(fromNode);
    set<int>::iterator fromNodes;
    while(!(parent_new.empty()))
    {
        parent_old = parent_new;
        parent_new.clear();
        for(fromNodes = parent_old.begin();fromNodes != parent_old.end();fromNodes++)
        {
            if(!(_allNodes[*fromNodes].Get_child().empty()))
            {
                set<int> tmpChild = _allNodes[*fromNodes].Get_child();
                set<int>::iterator ch;
                set<int> tmpDescandants;
                tmpDescandants.insert(tmpChild.begin(),tmpChild.end());
                for(ch = tmpChild.begin();ch != tmpChild.end();ch++)
                {
                    tmpDescandants.insert(_allNodes[ *ch ].Get_descendant().begin(),_allNodes[ *ch ].Get_descendant().end());
                }
                _allNodes[*fromNodes].Set_descendant(tmpDescandants);
            }
            if(!(_allNodes[*fromNodes].Get_parent().empty()))
            {
                parent_new.insert(_allNodes[*fromNodes].Get_parent().begin(),_allNodes[*fromNodes].Get_parent().end());
            }
        }
    }


}
*/
