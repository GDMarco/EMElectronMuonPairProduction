#include "EMElectronMuonPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

static const double mass_muon = 1.883531627e-28 * kilogram;
static const double mmc2 = mass_muon * c_squared;
static const double mec2 = mass_electron * c_squared;

EMElectronMuonPairProduction::EMElectronMuonPairProduction(ref_ptr<PhotonField> photonField, bool haveMuons, bool haveElectron, double thinning, double limit) {
    setPhotonField(photonField);
    setThinning(thinning);
    setLimit(limit);
    setHaveMuons(haveMuons);
    setHaveElectron(haveElectron);
}

void EMElectronMuonPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
    this->photonField = photonField;
    std::string fname = photonField->getFieldName();
    setDescription("EMElectronMuonPairProduction: " + fname);
    initRate(getDataPath("EMElectronMuonPairProduction/rate_" + fname + ".txt"));
    initCumulativeRate(getDataPath("EMElectronMuonPairProduction/cdf_" + fname + ".txt"));
    initInelasticity(getDataPath("EMElectronMuonPairProduction/inelasticity.txt")); // data from MUNHECA code.
}

void EMElectronMuonPairProduction::setHaveMuons(bool haveMuons) {
    this->haveMuons = haveMuons;
}

void EMElectronMuonPairProduction::setHaveElectron(bool haveElectron) {
    this->haveElectron = haveElectron;
}

void EMElectronMuonPairProduction::setLimit(double limit) {
    this->limit = limit;
}

void EMElectronMuonPairProduction::setThinning(double thinning) {
    this->thinning = thinning;
}

void EMElectronMuonPairProduction::initRate(std::string filename) {
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("EMElectronMuonPairProduction: could not open file " + filename);

    // clear previously loaded interaction rates
    tabEnergy.clear();
    tabRate.clear();

    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infile) {
                tabEnergy.push_back(pow(10, a) * eV);
                tabRate.push_back(b / Mpc);
            }
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infile.close();
}

void EMElectronMuonPairProduction::initCumulativeRate(std::string filename) {
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("EMElectronMuonPairProduction: could not open file " + filename);

    // clear previously loaded tables
    tabE.clear();
    tabs.clear();
    tabCDF.clear();
    
    // skip header
    while (infile.peek() == '#')
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');

    // read s values in first line
    double a;
    infile >> a; // skip first value
    while (infile.good() and (infile.peek() != '\n')) {
        infile >> a;
        tabs.push_back(pow(10, a) * eV * eV);
    }

    // read all following lines: E, cdf values
    while (infile.good()) {
        infile >> a;
        if (!infile)
            break;  // end of file
        tabE.push_back(pow(10, a) * eV);
        std::vector<double> cdf;
        for (int i = 0; i < tabs.size(); i++) {
            infile >> a;
            cdf.push_back(a / Mpc);
        }
        tabCDF.push_back(cdf);
    }
    infile.close();
}

void EMElectronMuonPairProduction::initInelasticity(std::string filename) {
    std::ifstream infile(filename.c_str());
    
    if (!infile.good())
        throw std::runtime_error("EMElectronMuonPairProduction: could not open file " + filename);
    
    tabsIn.clear();
    tab2Mu.clear();
    tabHeMu.clear();
    tabLeMu.clear();
    
    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b, c, d;
            infile >> a >> b >> c >> d;
            if (infile) {
                this->tabsIn.push_back(a * eV * eV);
                this->tab2Mu.push_back(b);
                this->tabHeMu.push_back(c);
                this->tabLeMu.push_back(d);
            }
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infile.close();
}
   
void EMElectronMuonPairProduction::performInteraction(Candidate *candidate) const {
    // scale particle energy instead of background photon energy
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);
    
    // it is only an energy loss, so see TPP
    // candidate->setActive(false);
    
    // check if secondary charged muons pair needs to be produced
    if (not haveMuons and not haveElectron)
        return;
    
    // check if in tabulated energy range
    if (E < tabE.front() or (E > tabE.back()))
        return;
    
    // sample the value of s
    Random &random = Random::instance();
    size_t i = closestIndex(E, tabE);  // find closest tabulation point
    size_t j = random.randBin(tabCDF[i]);
    double lo = std::max((2 * mmc2 + mec2) * (2 * mmc2 + mec2), tabs[j-1]);  // first s-tabulation point below min(s_kin) = (2 mm c^2 + me c^2)^2; ensure physical value
    double hi = tabs[j];
    double s = lo + random.rand() * (hi - lo);
    
    if (haveMuons) {
        // sample muon+ / muon- energy
        double inelasticityMuHe = interpolate(s, this->tabsIn, this->tabHeMu);
        double inelasticityMuLe = interpolate(s, this->tabsIn, this->tabHeMu);
        
        double EmuHe = inelasticityMuHe * E;
        double EmuLe = inelasticityMuLe * E;
        
        // for some backgrounds Ee=nan due to precision limitations.
        if (not std::isfinite(EmuLe) || not std::isfinite(EmuHe))
            return;
        
        // sample random position along current step
        Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
        
        // create a random number + or - 1 to randomly select the leading particle
        int randIntPM = (random.randUniform(0, 1) < 0.5) ? -1 : 1;
        
        // apply sampling
        if (random.rand() < pow(inelasticityMuHe, thinning)) {
            double w = 1. / pow(inelasticityMuHe, thinning);
            candidate->addSecondary(randIntPM * 13, EmuHe / (1 + z), pos, w, interactionTag); // positively charged muon
        }
        if (random.rand() < pow(1 - inelasticityMuLe, thinning)){
            if (EmuLe >= mmc2) {
                double w = 1. / pow(1 - inelasticityMuLe, thinning);
                candidate->addSecondary(randIntPM * -13, EmuLe / (1 + z), pos, w, interactionTag);
            }
        }
    }
    
    if (haveElectron) {
        
        double inelasticityEl = interpolate(s, this->tabsIn, this->tab2Mu);
        double Eel = E * (1 - inelasticityEl);
        
        if (not std::isfinite(Eel))
            return;
        
        candidate->current.setEnergy(Eel / (1 + z));
        
    } else {
        
        candidate->setActive(false);
        return;
        
    }
}

void EMElectronMuonPairProduction::process(Candidate *candidate) const {
    
    // check if electron
    if (std::abs(candidate->current.getId()) != 11)
        return;
    
    // scale particle energy instead of background photon energy
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);

    // check if in tabulated energy range
    if ((E < tabEnergy.front()) or (E > tabEnergy.back()))
        return;
    
    // interaction rate
    double rate = interpolate(E, tabEnergy, tabRate);
    rate *= pow_integer<2>(1 + z) * photonField->getRedshiftScaling(z);
    
    // run this loop at least once to limit the step size
    double step = candidate->getCurrentStep();
    Random &random = Random::instance();
    do {
        double randDistance = -log(random.rand()) / rate;
        // check for interaction; if it doesn't ocurr, limit next step
        if (step < randDistance) {
            candidate->limitNextStep(limit / rate);
        } else {
            performInteraction(candidate);
            return;
        }
        step -= randDistance;
    } while (step > 0.);

}

void EMElectronMuonPairProduction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string EMElectronMuonPairProduction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace crpropa
