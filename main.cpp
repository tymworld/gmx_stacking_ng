#include <iostream>
#include "stacking.h"

int main(int argc, char * argv[]) {
    // insert code here...
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<stacking>(argc, argv);

}

