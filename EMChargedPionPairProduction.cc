#include "EMChargedPionPairProduction.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <string>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace crpropa {

// mass W boson
static const double massWGeV = 80.379;
static const double mWkg = 1.43288582 * 1e-25; // kg
// static const double mW2 = (mWkg*c_light**2.) ** 2; //squared W mass [J^2/c^4]

// mass charged pion
static const double massPionGeV = 139.57039e-3; // GeV/c^2
static const double mass_cpion = mWkg * massPionGeV / massWGeV * kilogram;
static const double mcpic2 = mass_cpion * c_squared;

EMChargedPionPairProduction::EMChargedPionPairProduction(ref_ptr<PhotonField> photonField, bool havePions, double thinning, double limit) {
    setPhotonField(photonField);
    setThinning(thinning);
    setLimit(limit);
    setHavePions(havePions);
}

void EMChargedPionPairProduction::setPhotonField(ref_ptr<PhotonField> photonField) {
    this->photonField = photonField;
    std::string fname = photonField->getFieldName();
    setDescription("EMChargedPionPairProduction: " + fname);
    initRate(getDataPath("EMChargedPionPairProduction/rate_" + fname + ".txt"));
    initCumulativeRate(getDataPath("EMChargedPionPairProduction/cdf_" + fname + ".txt"));
    initInelasticity(getDataPath("EMChargedPionPairProduction/inelasticity.txt")); // data from MUNHECA code.
}

void EMChargedPionPairProduction::setHavePions(bool havePions) {
    this->havePions = havePions;
}

void EMChargedPionPairProduction::setLimit(double limit) {
    this->limit = limit;
}

void EMChargedPionPairProduction::setThinning(double thinning) {
    this->thinning = thinning;
}

void EMChargedPionPairProduction::initRate(std::string filename) {
    std::ifstream infile(filename.c_str());

    if (!infile.good())
        throw std::runtime_error("EMChargedPionPairProduction: could not open file " + filename);

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

void EMChargedPionPairProduction::initCumulativeRate(std::string filename) {
    std::ifstream infile(filename.c_str());

    if (!infile.good())
        throw std::runtime_error("EMChargedPionPairProduction: could not open file " + filename);

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

void EMChargedPionPairProduction::initInelasticity(std::string filename) {
    std::ifstream infile(filename.c_str());

    if (!infile.good())
        throw std::runtime_error("EMChargedPionPairProduction: could not open file " + filename);
    
    tabsIn.clear();
    tabInelasticity.clear();
    
    while (infile.good()) {
        if (infile.peek() != '#') {
            double a, b;
            infile >> a >> b;
            if (infile) {
                tabsIn.push_back(a * eV * eV);
                tabInelasticity.push_back(b); // inelasticity is dimensionless
            }
        }
        infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
    }
    infile.close();
}

// maybe to rewrite according to the inelasticity...
// class PPSecondariesEnergyDistribution
   
void EMChargedPionPairProduction::performInteraction(Candidate *candidate) const {
    // scale particle energy instead of background photon energy
    double z = candidate->getRedshift();
    double E = candidate->current.getEnergy() * (1 + z);

    // cosmic ray photon is lost after interacting
    candidate->setActive(false);

    // check if secondary charged pions pair needs to be produced
    if (not havePions)
        return;

    // check if in tabulated energy range
    if (E < tabE.front() or (E > tabE.back()))
        return;

    // sample the value of s
    Random &random = Random::instance();
    size_t i = closestIndex(E, tabE);  // find closest tabulation point
    size_t j = random.randBin(tabCDF[i]);
    double lo = std::max(4 * mcpic2 * mcpic2, tabs[j-1]);  // first s-tabulation point below min(s_kin) = (2 me c^2)^2; ensure physical value
    double hi = tabs[j];
    double s = lo + random.rand() * (hi - lo);

    // sample pion+ / pion- energy
    double inelasticity = interpolate(s, tabsIn, tabInelasticity); // it should come from the tables!
    double EpiP = inelasticity * E; // it is assumed that the leading particle is the positively charged pion
    
    double EpiM = E - EpiP;
    // double f = EpiP / E; // this is inelasticity itself

    // for some backgrounds Ee=nan due to precision limitations.
    if (not std::isfinite(EpiP) || not std::isfinite(EpiM))
        return;

    // sample random position along current step
    Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
    // apply sampling
    if (random.rand() < pow(inelasticity, thinning)) {
        double w = 1. / pow(inelasticity, thinning);
        candidate->addSecondary(211, EpiP / (1 + z), pos, w, interactionTag); // positively charged pion
    }
    if (random.rand() < pow(1 - inelasticity, thinning)){
        double w = 1. / pow(1 - inelasticity, thinning);
        candidate->addSecondary(-211, EpiM / (1 + z), pos, w, interactionTag);
    }
}

void EMChargedPionPairProduction::process(Candidate *candidate) const {
    // check if photon
    if (candidate->current.getId() != 22)
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

void EMChargedPionPairProduction::setInteractionTag(std::string tag) {
    interactionTag = tag;
}

std::string EMChargedPionPairProduction::getInteractionTag() const {
    return interactionTag;
}

} // end namespace crpropa
