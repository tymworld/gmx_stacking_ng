//
// Created by Yiming Tang on 2020/9/24.
//

#include "stacking.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <gromacs/fileio/matio.h>
#include <gromacs/fileio/gmxfio.h>

float stacking::calDistance(const coordinate atom11, const coordinate atom12, const coordinate atom13,
                            const coordinate atom21, const coordinate atom22, const coordinate atom23)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x + atom13.x) / 3;
    y_1 = (atom11.y + atom12.y + atom13.y) / 3;
    z_1 = (atom11.z + atom12.z + atom13.z) / 3;

    x_2 = (atom21.x + atom22.x + atom23.x) / 3;
    y_2 = (atom21.y + atom22.y + atom23.y) / 3;
    z_2 = (atom21.z + atom22.z + atom23.z) / 3;

    return (float)sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
}

float stacking::calAngle(coordinate atom11, coordinate atom12, coordinate atom13,
                         coordinate atom21, coordinate atom22, coordinate atom23) {
    float a1 = ((atom12.y - atom11.y) * (atom13.z - atom11.z) - (atom12.z - atom11.z) * (atom13.y - atom11.y));
    float b1 = ((atom12.z - atom11.z) * (atom13.x - atom11.x) - (atom12.x - atom11.x) * (atom13.z - atom11.z));
    float c1 = ((atom12.x - atom11.x) * (atom13.y - atom11.y) - (atom12.y - atom11.y) * (atom13.x - atom11.x));

    float a2 = ((atom22.y - atom21.y) * (atom23.z - atom21.z) - (atom22.z - atom21.z) * (atom23.y - atom21.y));
    float b2 = ((atom22.z - atom21.z) * (atom23.x - atom21.x) - (atom22.x - atom21.x) * (atom23.z - atom21.z));
    float c2 = ((atom22.x - atom21.x) * (atom23.y - atom21.y) - (atom22.y - atom21.y) * (atom23.x - atom21.x));


    float angle = (float)(acos((a1 * a2 + b1 * b2 + c1 * c2) / sqrt(a1 * a1 + b1 * b1 + c1 * c1) /
                               sqrt(a2 * a2 + b2 * b2 + c2 * c2)) / (M_PI) * 180);

    return min(abs(angle - 0), abs(180 - angle));

}


stacking::stacking(): cutoff_(1.5), maxDistance_(1.2), maxAngle_(90.0), stepDistance_(0.01), stepAngle_(1.0)
{
    registerAnalysisDataset(&data_probability_, "probability");
}

void stacking::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
            {
                    "Next Generation Stacking Free Energy Surface Calculator by Yiming Tang @ Fudan. \n",
                    "This Tool calculates the Joint Probabilit Density (JPD) of two molecules",
                    "as a function of angle (0-90) - distance (0-1.2). The tool outputs the free",
                    "energy surface defined by E = - R * T * log(H), where H is the probability.",
                    "This tool takes two selection group named ref and sel, each should contains",
                    "rings where each ring contains specifically three atoms.",
                    "Angles and Distances will be calculated between each ring in ref and that in sel. ",
                    "If the two groups are identical, the program automatically strips out same ring. \n",
                    "This program is a totally rewritten version of the original gmx_stacking program.",
                    "As the original program can perform a bit more task than this program",
                    "(non-cg benzene rings with >3 atoms), we choose not to update that program,",
                    "but write a new program (this one) instead."
            };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("energy").filetype(eftUnknown).legacyType(efXPM).outputFile()
                               .store(&fnEnergySurface_).defaultBasename("energy")
                               .description("Energy Surface xpm map").required());

    options->addOption(FileNameOption("raw-energy").filetype(eftGenericData).outputFile()
                               .store(&fnEnergySurfaceRaw_)
                               .description("Ras Energy Surface Data for further analysis"));

    options->addOption(FileNameOption("raw-probability").filetype(eftGenericData).outputFile()
                               .store(&fnPropability_)
                               .description("Raw Probability Data for further analysis"));

    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                               .description("Rings whose maximum atom-wise distance are beyond this cutoff is not considered contacted"));

    options->addOption(DoubleOption("maxD").store(&maxDistance_).defaultValue(1.2)
                               .description("Maxmimum Centroid Distance to output."));

    options->addOption(DoubleOption("stepD").store(&stepDistance_).defaultValue(0.01)
                               .description("Step of Centroid Distance to output."));

    options->addOption(DoubleOption("maxA").store(&maxAngle_).defaultValue(90)
                               .description("Maxmimum Centroid Angle to output, default is 90 for plane and 180 for vector."));

    options->addOption(DoubleOption("stepA").store(&stepAngle_).defaultValue(1)
                               .description("Step of Centroid Angle to output."));

    options->addOption(SelectionOption("ref").store(&ref_).required()
                       .description("Groups containing rings to calculate angle and distance from. Three atoms per ring."));

    options->addOption(SelectionOption("sel").store(&sel_).required()
                               .description("Groups containing rings to calculate angle and distance to. Three atoms per ring."));

    options->addOption(DoubleOption("temperature").store(&temperature_).defaultValue(298)
                               .description("Temperature in K for energy calculation"));

    options->addOption(BooleanOption("inter_molecule").store(&inter_molecule_).defaultValue(false)
                       .description("Whether only contact between different molecule is calculated. MUST SPECIFY num_ring!!!"));

    options->addOption(BooleanOption("intra_molecule").store(&intra_molecule_).defaultValue(false)
                               .description("Whether only contact between the same molecule is calculated. MUST SPECIFY num_ring!!!"));

    options->addOption(IntegerOption("num_ring").store(&ring_number_in_molecule_).defaultValue(1)
                       .description("Number of rings in each molecule. Specify this without enabling inter_molecule is useless."));

}

void stacking::initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top)
{

    probability_ = new unsigned long*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        probability_[i] = new unsigned long[(int)(maxAngle_ / stepAngle_)];
    }

    for(int i = 0; i< ((int)(maxDistance_ / stepDistance_)); i++)
    {
        for(int j = 0; j < (int)(maxDistance_ / stepDistance_); j++)
        {
            probability_[i][j] = 0;
        }
    }

    if(ref_.atomCount() % 3 or sel_.atomCount() % 3)
    {
        cerr << "ERROR: At least one of our selection contains atoms not a multiple of three." << endl;
        exit(1);
    }

    ring_has_contact_ = new bool*[ref_.atomCount() / 3];
    for(int i = 0; i < ref_.atomCount() / 3; i++)
    {
        ring_has_contact_[i] = new bool[sel_.atomCount() / 3];
    }

    data_probability_.setDataSetCount(1);
    data_probability_.setColumnCount(0,1);

    if((not inter_molecule_ and not intra_molecule_) and ring_number_in_molecule_ > 1)
    {
        cerr << "ERROR: Specify num_ring without -inter_molecule or -intra_molecule" << endl;
        exit(1);
    }

    identical_selection = true;
    if(ref_.atomCount() != sel_.atomCount())
    {
        identical_selection = false;
    }
    else
    {
        for(int i = 0; i < ref_.atomCount(); i++)
        {
            if(ref_.atomIndices()[i] != sel_.atomIndices()[i])
            {
                identical_selection = false;
            }
        }
    }

    if((inter_molecule_ or intra_molecule_) and not identical_selection)
    {
        cerr << "Do not use -inter/intra_molecule if the two selections ref_ and sel_ are different. " << endl;
        cerr << "This will produce artificial results." << endl;
        cerr << "The program will now end." << endl;
        exit(1);
    }

    if(inter_molecule_)
    {
        cout << "*****************************************" << endl;
        cout << "**        VERY IMPORTANT NOTICE        **" << endl;
        cout << "** You have specified -inter-molecule_ **" << endl;
        cout << "** You should also specify -num_ring   **" << endl;
        cout << "** Now a molecule contains " << setw(3) << ring_number_in_molecule_ << " rings!! **" << endl;
        cout << "*****************************************" << endl;
    }

    if(intra_molecule_)
    {
        cout << "*****************************************" << endl;
        cout << "**        VERY IMPORTANT NOTICE        **" << endl;
        cout << "** You have specified -intra-molecule_ **" << endl;
        cout << "** You should also specify -num_ring   **" << endl;
        cout << "** Now a molecule contains " << setw(3) << ring_number_in_molecule_ << " rings!! **" << endl;
        cout << "*****************************************" << endl;
    }

    top_   = top.topology();
    atoms_ = top.topology()->atoms;

    nb_.setCutoff(cutoff_);

}

void stacking::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata)
{

    AnalysisDataHandle dhProbability = pdata->dataHandle(data_probability_);
    dhProbability.startFrame(frnr, fr.time);

    // We first clear all true values in ring_has_contact_ parameter.

    for(int i = 0; i < ref_.atomCount() / 3; i++)
    {
        for(int j = 0; j < sel_.atomCount() / 3; j++)
        {
            ring_has_contact_[i][j] = false;
        }
    }

    // Begin analysis

    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, ref_);
    gmx::ArrayRef<float const[3]>::iterator iter_coordinate_sel;
    for(iter_coordinate_sel = sel_.coordinates().begin(); iter_coordinate_sel != sel_.coordinates().end(); iter_coordinate_sel++)
    {
        AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(*iter_coordinate_sel);
        AnalysisNeighborhoodPair pair;

        while(pairSearch.findNextPair(&pair))
        {
            int atom_index_in_ref = pair.refIndex();
            int atom_index_in_sel = pair.testIndex();
            int ring_index_in_ref = floor(atom_index_in_ref / 3.0);
            int ring_index_in_sel = floor(atom_index_in_sel / 3.0);
            ring_has_contact_[ring_index_in_ref][ring_index_in_sel] = true;
        }
    }

    // We now iterate through contacted rings to calculate their centroid distances and angles.

    for(int ring_index_in_ref = 0; ring_index_in_ref < ref_.atomCount() / 3; ring_index_in_ref++)
    {
        for(int ring_index_in_sel = 0; ring_index_in_sel < sel_.atomCount() / 3; ring_index_in_sel++)
        {
            if(inter_molecule_
            and floor(ring_index_in_ref / float(ring_number_in_molecule_))
             == floor(ring_index_in_sel / float(ring_number_in_molecule_)))
            {
                continue;
            }

            if(intra_molecule_
               and floor(ring_index_in_ref / float(ring_number_in_molecule_))
                   != floor(ring_index_in_sel / float(ring_number_in_molecule_)))
            {
                continue;
            }

            if(identical_selection and ring_index_in_ref == ring_index_in_sel)
            {
                continue;
            }

            coordinate molecule1[3];
            coordinate molecule2[3];
            for(int i = 0; i < 3; i++)
            {
                molecule1[i].x = ref_.coordinates()[ring_index_in_ref * 3 + i][0];
                molecule1[i].y = ref_.coordinates()[ring_index_in_ref * 3 + i][1];
                molecule1[i].z = ref_.coordinates()[ring_index_in_ref * 3 + i][2];
                molecule2[i].x = sel_.coordinates()[ring_index_in_sel * 3 + i][0];
                molecule2[i].y = sel_.coordinates()[ring_index_in_sel * 3 + i][1];
                molecule2[i].z = sel_.coordinates()[ring_index_in_sel * 3 + i][2];
            }

            float tempDistance = calDistance(molecule1[0], molecule1[1], molecule1[2],
                                             molecule2[0], molecule2[1], molecule2[2]);

            float tempAngle = calAngle(molecule1[0], molecule1[1], molecule1[2],
                                       molecule2[0], molecule2[1], molecule2[2]);

            if (tempDistance <= maxDistance_ && tempAngle <= maxAngle_ ) {
                probability_[(int) floor(tempDistance / stepDistance_)][(int) floor(tempAngle / stepAngle_)] += 1;
            }
        }
    }

    dhProbability.finishFrame();
}

void stacking::finishAnalysis(int /*nframes*/) {}


void stacking::writeOutput()
{

    real DistanceVector[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_) ; i++)
    {
        DistanceVector[i] = real(i * stepDistance_);
    }

    real AngleVector[(int)(maxAngle_ / stepAngle_)];

    for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
    {
        AngleVector[j] = real(j * stepAngle_);
    }


    if       (fnEnergySurface_.empty())         {fnEnergySurface_ = "energy.xpm";}
    else if  (fnEnergySurface_.compare(".xpm")) {}
    else                                        {fnEnergySurface_ += ".xpm";}

    // Construct matrix probability
    real **matProbability = new real*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matProbability[i] = new real[(int)(maxAngle_ / stepAngle_)];
    }

    real **matEnergy = new real*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matEnergy[i] = new real[(int)(maxAngle_ / stepAngle_)];
    }

    int scale = (int)(sel_.atomCount() / 3) * (int)(ref_.atomCount() / 3);

    for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
        {
            matProbability[i][j] = probability_[i][j] / (real)data_probability_.frameCount();
            //matProbability[i][j] = avem_probability_->average(i, j) ;
            //cout << i << '\t' << j << '\t' << avem_probability_->average(i, j) << '\t' << avem_probability_->average(i, j) / (float)scale << '\t' << matProbability[i][j] << endl;
            matEnergy[i][j]      = real(- 8.31447 * temperature_ * log(matProbability[i][j] / (float)scale) * 0.0002389);
        }
    }

    t_rgb rlo, rhi;
    rlo.r = 0.0; rlo.g = 0.0; rlo.b = 1.0;
    rhi.r = 1.0; rhi.g = 0.0; rhi.b = 0.0;
    int nlevels  = 400;

    FILE *fpEnergySurface;
    fpEnergySurface = fopen(fnEnergySurface_.c_str(), "w");



    write_xpm(fpEnergySurface, 0, "Free Energy Surface"
            , "Contact Probability", "Distance", "angle"
            , (int)(maxDistance_ / stepDistance_), (int)(maxAngle_ / stepAngle_), DistanceVector, AngleVector
            , matEnergy, 0, 20, rlo, rhi, &nlevels);

    fclose(fpEnergySurface);

    if (!fnEnergySurfaceRaw_.empty())
    {
        if (fnEnergySurfaceRaw_.compare(".dat"))    {}
        else                                        {fnEnergySurfaceRaw_ += ".dat";}

        FILE *fpEnergySurfaceRaw;
        fpEnergySurfaceRaw = fopen(fnEnergySurfaceRaw_.c_str(), "w" );

        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
            {
                fprintf(fpEnergySurfaceRaw, "%f ", matEnergy[i][j]);
            }
            fprintf(fpEnergySurfaceRaw, "\n");
        }
        fclose(fpEnergySurfaceRaw);


    }

    if (!fnPropability_.empty())
    {
        if (fnPropability_.compare(".dat"))    {}
        else                                        {fnPropability_ += ".dat";}

        FILE *fpProbability;
        fpProbability = fopen(fnPropability_.c_str(), "w" );

        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
            {
                fprintf(fpProbability, "%f ", matProbability[i][j]);
            }
            fprintf(fpProbability, "\n");
        }
        fclose(fpProbability);

    }

}