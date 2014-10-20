#include "main.h"

int main(int argc, char* argv[])
{
    //read in the command line argument
    if (argc < 11)
    {
        cerr<<"Not enough input, the format should be:"<<argv[0]<< " #NODES #OBS #iter #chains cutoff #cores temperature depth DATA_path PARAM_path output_path"<<endl;
    }
    // Number of nodes (parameters) in the data
    string Node_temp = argv[1];
    int NODES = atoi(Node_temp.c_str());
    // Number of observations in the data
    string obs_temp = argv[2];
    int OBSERVATIONS = atoi(obs_temp.c_str());
    // number of maximum iteration in each chain
    string iter_temp = argv[3];
    int NUMBER_OF_ITERATION = atoi(iter_temp.c_str());
    // Number of chains to grow
    string chain_temp = argv[4];
    int NUMBER_OF_CHAIN = atoi(chain_temp.c_str());
    // The cut off point where MIT test determine dependency
    string cut_temp = argv[5];
    double CUT_OFF = atof(cut_temp.c_str());
    // this is somehow an obsolete parameter, may remove it in later version.
//    string independent_temp = argv[6];
//    double PROBABILITY_OF_INDEPENDENT = atof(independent_temp.c_str());
    // Number of cores to use, it gets replaced by number of cores available if the later one is smaller.
    string cores_temp = argv[6];
    int NUMBER_OF_CORES = atoi(cores_temp.c_str());
    // temperature, used to control the greedy level.
    string temperature_temp = argv[7];
    double TEMPERATURE = atof(temperature_temp.c_str());

    string depth_temp = argv[8];
    int DEPTH = atof(depth_temp.c_str());

    // path to the data file
    char* DATA_NAME = argv[9];
    // path to the parameter file
    char* PARA_RANGE = argv[10];

//    char* SKELTON_FILE = argv[10];
    // path to the output file.
    string OUTPUT_FILE = argv[11];

// First get the number of cores in current system
    int iCPU = omp_get_num_procs();
    // if number of available cores smaller than required then use maximum cores
    if (iCPU < NUMBER_OF_CORES )
    {
        NUMBER_OF_CORES = iCPU;
    }
    time_t iniStart,iniEnd,procStart,procEnd;

    DATA data1(OBSERVATIONS,NODES);
    data1.ReadData(DATA_NAME);
    data1.ReadParam(PARA_RANGE);
 //   data1.ReadSkeleton(SKELTON_FILE);
    TEMPERATURE*=OBSERVATIONS;

    time(&iniStart);
    SMC smc1(NODES,NUMBER_OF_CHAIN);
    //cout<<data1.Get_param().size()<<endl;
    smc1.Initialize(data1, CUT_OFF,NUMBER_OF_CHAIN);
    smc1.DFSummary(OUTPUT_FILE,NODES);
    time(&iniEnd);
    time(&procStart);
    smc1.Update(data1,NUMBER_OF_ITERATION,NUMBER_OF_CHAIN,NUMBER_OF_CORES,TEMPERATURE,DEPTH);
    time(&procEnd);
    ofstream outFile;
    outFile.open((OUTPUT_FILE+"_summary.txt").c_str());
    outFile<<"Elasped time for initialization is "<< difftime(iniEnd,iniStart)<<endl;
    outFile<< "Time elasped is "<< difftime(procEnd,procStart)<<endl;
    outFile.close();

    smc1.Summary(OUTPUT_FILE,NUMBER_OF_CHAIN,NODES);


    return 0;
}
