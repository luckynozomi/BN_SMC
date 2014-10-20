#include "cpd.h"

CPD::CPD()
{
    //ctor
}

CPD::~CPD()
{
    //dtor
}


int* CPD::_CumProduct(const int* leng,const int length,int* res)
{
    for(int i=0; i< length; i++)
    {
        res[0]*=leng[i];
    }
    for(int i=0; i<length; i++)
    {
        res[i+1]=res[i]/leng[i];
    }
    return res;
}
// @example: say this node id is 1, parents are 2,3. Id 1 has 2 categories, 2,3 have 3,4 categories respectively
// so leng is (3,4,2)
// length is 2+1=3
// res should be (1,x,x,x) where x is anything, not necessary initialized.
// the returned array is (3*4*2,4*2,2,1)=(24,8,2,1)

int CPD::_IndexSearchTarget(vector<int>& searchTarget, int* param)
{
    int index = 0;
    for(unsigned int i = 0; i< searchTarget.size(); i++)
    {
        index += (param[i+1]) * (searchTarget[i]-1);
    }

    return index;

}
//@algorithm: The idea is really like binary-> dec conversion,
// let us still use the example from previous function,
// let the vector A=(a,b,c) denotes one instance of the observations, where a is the value for node 2which is the parent of node 1,
// so a takes values 0,1,2; similarly b takes values 0,1,2,3 and c takes value 0,1
// let us index A=(0,0,0) as 0
// and (0,0,1) as 1
// so we have the following:
/*
realization index
(0,0,0)     0=0*(4*2)+0*(2)+0
(0,0,1)     1=0*(4*2)+0*(2)+1
(0,1,0)     2=0*(4*2)+1*(2)+0
(0,1,1)     3=0*(4*2)+1*(2)+1
(0,2,0)     4=0*(4*2)+2*(2)+0
...         ...
(1,0,0)     8=1*(4*2)+0*(2)+0
(1,1,0)     10=1*(4*2)+1*(2)+0
*/
// in general, say A=(a1,a2,a3,...,an) is a realization of the local structure
// and supposed that each entry has categories (b1,b2,b3,...,bn)
// then its index is I=a1*(b2*b3*...*bn)+a2*(b3*...*bn)+...+a_{n-1}*(bn)+an.
// Here we should notice that the param parameter in the function is exactly the one obtained
// by _cumProduct( start from the second entry).


void CPD::UpdateCPD(const vector<int>& parent, DATA& data)
{

    //Since BIC is locally consistant, so at each iteration we could update only the "toNode"
    //vector<int> temp_par=targetNode.Get_parent();// retrive the parent information.
    int temp_size=parent.size(); // number of parents
    int numberOfParaEach[temp_size+1];// the array restores the number of categories for each parent and the "toNode"

    //This loop is intended to give the array numberOfParaEach values for the first n-1 entries which are the parent categories.
    for(int i=0; i< temp_size; i++)
    {
        numberOfParaEach[i]=data.Get_param(parent[i]);
    }
    // this part is used to give the last entry of the array, which is the number of category
    //_numberOfParameter = numberOfParam;
    numberOfParaEach[temp_size]=_numberOfParameter;
    // c1,c2,c3,...,c_n,T

    // This chuck of code is used to generate the array that used to index each observation and the number of parent configurates, as well as the total configurations.
    int temp_param[temp_size+2];
    temp_param[0]=1;
    int* param;
    param = _CumProduct(numberOfParaEach,temp_size+1, temp_param);

    _numberOfParentConfigure = param[0]/_numberOfParameter;
    _totalParam = param[0];
    // this part is used to count the occurences of each incident
    //vector<int> temp_conditionCount;
    //vector<double> temp_probability;
    _conditionCount.clear();
    _probability.clear();
    _conditionCount.reserve(_totalParam);
    _probability.reserve(_totalParam);
    // fill out the _conditionCount
    for(int i=0; i<_totalParam; i++)//initialize, param[0] saves the total number of configurations.
    {
        _probability.push_back(1.0e-14);// avoid 0 probability which may cause the BIC gets to 0.
        _conditionCount.push_back(0);
    }

    //_conditionCount=temp_conditionCount;
    //_probability = temp_probability;

    //temp_par.push_back( targetNode.Get_nodeID() );// in addition to the parent IDs, add the "toNode" id so we could get the corresponding columns of the data.
    //int data_obs=data.Get_NumberOfObservations();// number of observations.
    //vector< vector <int> > temp_data;
    //temp_data=data.Get_data();// store data for counting use.

    for(int i=0; i< data.Get_NumberOfObservations(); i++)
    {
        vector<int> searchTarget(temp_size+1);
        // create the searching target which is one observation(incident) of the local strucutre.
        for(int j=0; j< temp_size; j++)
        {
            //searchTarget.push_back(temp_data[i][ temp_par[j] ]);
            searchTarget[j] = data.Get_data(i, parent[j]);
        }
        searchTarget[temp_size] = data.Get_data(i, _correspondingNode);
        //int index;
        //index = _IndexSearchTarget(searchTarget,param);// index the target
        _conditionCount[ _IndexSearchTarget(searchTarget,param) ]++;// set the count to be one more.
    }





    //calculate _probability

    int numberOfconditions=_totalParam/_numberOfParameter;


    vector<int> cumsum(numberOfconditions);// normalizer

/*    for(int i=0; i<numberOfconditions; i++)
    {
        cumsum[i]=0;
    }
*/
    for(int i=0; i<_totalParam; i++)
    {
        cumsum[ ( i / (_numberOfParameter) ) ]+= _conditionCount[i];//numberOfParaEach[temp_size] is the number of categories of the "toNode"
        // if there are 3 categories for the "toNode", then 0/3,1/3,2/3 will all lead to 0, and 3/3,4/3,5/3 will lead to 1.
    }

    for(int i=0; i<_totalParam; i++)
    {

        _probability[i]+= _conditionCount[i] /( (1.0e-12)+cumsum[ ( i / (_numberOfParameter) ) ]);
        // the 1.0e-12 is used to avoid 0/0.

    }

}

void CPD::UpdateBIC(int obs)
{
    _BIC = - 0.5*_numberOfParentConfigure*(_numberOfParameter-1)*log( 1.0e-12+obs );
    for(unsigned int i = 0; i< _probability.size(); i++)
    {
        _BIC += _conditionCount[i]*log(_probability[i]);
    }
}



