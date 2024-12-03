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
 @class EMElectronMuonPairProduction
 @brief Charged pion-pair production of photons with background photons.

 This module simulates electron muon pair production of cosmic ray electrons with background photons:
 electron + gamma_b -->  electron + mu+ + mu- (Breit-Wheeler process).
 The resulting charged muon pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked.
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class EMElectronMuonPairProduction : public Module
{
private:
    ref_ptr<PhotonField> photonField;     // target photon field
    bool haveMuons;                    // add secondary muons to simulation
    bool haveElectron;               // add secondary electron to simulation
    double limit;                // limit the step to a fraction of the mean free path
    double thinning;                    // factor of the thinning (0: no thinning, 1: maximum thinning)
    std::string interactionTag = "ElMuPP";

    // tabulated interaction rate 1/lambda(E)
    std::vector<double> tabEnergy;  //!< electron energy in [J]
    std::vector<double> tabRate;  //!< interaction rate in [1/m]
    
    // tabulated CDF(s_kin, E) = cumulative differential interaction rate
    std::vector<double> tabE;  //!< electron energy in [J]
    std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
    std::vector< std::vector<double> > tabCDF;  //!< cumulative interaction rate
    
    std::vector<double> tabsIn; //!< inelasticity s table, given in [eV**2]
    std::vector<double> tab2Mu; //!< inelasticity table for the energy fraction of the two muons
    std::vector<double> tabHeMu; //!< inelasticity table for the leading muon
    std::vector<double> tabLeMu; //!< inelasticity table for the low-energy muon

public:
    /** Constructor
     @param photonField        target photon field
     @param haveMuons   if true, add secondary muons as candidates
     @param haveElectron   if true, add the secondary electron as candidate
     @param thinning        weighted sampling of secondaries (0: all particles are tracked; 1: maximum thinning)
     @param limit            step size limit as fraction of mean free path
     */
    EMElectronMuonPairProduction(ref_ptr<PhotonField> photonField, bool haveMuons = false, bool haveElectron = false, double thinning = 0, double limit = 0.1);
    
    // set the target photon field
    void setPhotonField(ref_ptr<PhotonField> photonField);

    // decide if secondary pions are added to the simulation
    void setHaveMuons(bool haveMuons);
    
    // decide if secondary pions are added to the simulation
    void setHaveElectron(bool haveElectron);

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
