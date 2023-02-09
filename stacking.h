//
// Created by Yiming Tang on 2020/9/24.
//

#ifndef GMX_STACKING_STACKING_H
#define GMX_STACKING_STACKING_H

#include "gromacs/trajectoryanalysis.h"
#include <string>
#include <vector>

using namespace std;
using namespace gmx;

struct coordinate
{
    float x;
    float y;
    float z;
};

class stacking: public TrajectoryAnalysisModule
{
public:
    stacking();

    virtual void initOptions(IOptionsContainer          *options,
                             TrajectoryAnalysisSettings *settings);

    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation        &top);

    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);

    virtual void finishAnalysis(int nframes);

    virtual void writeOutput();

private:

    static float calDistance(coordinate atom11, coordinate atom12, coordinate atom13,
                             coordinate atom21, coordinate atom22, coordinate atom23);

    static float calAngle(coordinate atom11, coordinate atom12, coordinate atom13,
                          coordinate atom21, coordinate atom22, coordinate atom33);

    bool identical_selection;


    gmx::AnalysisData             data_probability_;
    gmx::Selection                ref_;
    gmx::Selection                sel_;

    gmx::AnalysisNeighborhood     nb_;

    double cutoff_;

    unsigned long **probability_;

    std::string     fnEnergySurface_;
    std::string     fnEnergySurfaceRaw_;
    std::string     fnPropability_;

    double          maxDistance_;
    double          stepDistance_;
    double          maxAngle_;
    double          stepAngle_;

    double          temperature_;

    // The following parameters are only valid if the sel_ and ref_ selection groups are identical.
    // And, if the inter_molecule_ is true, the ring_number_in_molecule must be specified.
    // Vice versa, if the ring_number_in_molecule_ is specified, the inter_molecule_ must be true.

    bool            inter_molecule_;
    bool            intra_molecule_;
    int             ring_number_in_molecule_;

    //t_topology     *top_;
    const t_atoms         *atoms_;

    bool **ring_has_contact_;

};

#endif //GMX_STACKING_STACKING_H
