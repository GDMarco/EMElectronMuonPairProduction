#include "crpropa/Module.h"
#include "crpropa/PhotonBackground.h"

#include <fstream>
#include <cmath>

namespace crpropa {
/**
 * \addtogroup EnergyLosses
 * @{
 */
 
/**
 @class EMChargedPionPairProduction
 @brief Charged pion-pair production of photons with background photons.

 This module simulates charged pion-pair production of cosmic ray photons with background photons:
 gamma + gamma_b --> pi+ + pi- (Breit-Wheeler process).
 The resulting charged pion pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked.
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class EMChargedPionPairProduction : public Module
{
private:
    ref_ptr<PhotonField> photonField;     // target photon field
    bool havePions;                    // add secondary pions to simulation
    double limit;                        // limit the step to a fraction of the mean free path
    double thinning;                    // factor of the thinning (0: no thinning, 1: maximum thinning)
    std::string interactionTag = "EMCPPP";

    // tabulated interaction rate 1/lambda(E)
    std::vector<double> tabEnergy;  //!< pion energy in [J]
    std::vector<double> tabRate;  //!< interaction rate in [1/m]
    
    // tabulated CDF(s_kin, E) = cumulative differential interaction rate
    std::vector<double> tabE;  //!< pion energy in [J]
    std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
    std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate
    
    std::vector<double> tabsIn; //!< inelasticity s table, given in [eV**2]
    std::vector<double> tabInelasticity; //!< inelasticity values
public:
    /** Constructor
     @param photonField        target photon field
     @param havePions   if true, add secondary pions as candidates
     @param thinning        weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
     @param limit            step size limit as fraction of mean free path
     */
    EMChargedPionPairProduction(ref_ptr<PhotonField> photonField, bool havePions = false, double thinning = 0, double limit = 0.1);
    
    // set the target photon field
    void setPhotonField(ref_ptr<PhotonField> photonField);

    // decide if secondary pions are added to the simulation
    void setHavePions(bool havePions);

    /** Limit the propagation step to a fraction of the mean free path
     * @param limit fraction of the mean free path
     */
    void setLimit(double limit);

    /** Apply thinning with a given thinning factor
     * @param thinning factor of thinning (0: no thinning, 1: maximum thinning)
     */
    void setThinning(double thinning);

    /** set a custom interaction tag to trace back this interaction
     * @param tag string that will be added to the candidate and output
     */
    void setInteractionTag(std::string tag);
    std::string getInteractionTag() const;

    void initRate(std::string filename);
    void initCumulativeRate(std::string filename);
    void initInelasticity(std::string filename);

    void performInteraction(Candidate *candidate) const;
    void process(Candidate *candidate) const;
};
/** @}*/

} // end namespace crpropa
