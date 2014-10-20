#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <boost/tokenizer.hpp>
#include <math.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <omp.h>

#include "node.h"
#include "bn.h"
#include "cpd.h"
#include "data.h"
#include "smc.h"
#include "constraint.h"
#include "edge.h"

#include <time.h>

using namespace std;
