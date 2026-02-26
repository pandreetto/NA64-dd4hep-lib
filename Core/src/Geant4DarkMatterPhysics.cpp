#include "Geant4DarkMatterPhysics.h"

#include <G4ProcessManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4PhysicsListHelper.hh>

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

#include "DarkMatterAnnihilation.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "DarkZAnnihilation.hh"
#include "DarkMuPhilicScalars.hh"
#include "DarkMuPhilicPseudoScalars.hh"
#include "ALP.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkPseudoScalarsAnnihilation.hh"
#include "DarkAxialsAnnihilation.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkMassSpin2.hh"
#include "DarkMassSpin2Annihilation.hh"
#include "DarkAxials.hh"

#include "DMProcessDMBrem.hh"
#include "DMProcessPrimakoffALP.hh"
#include "DMProcessAnnihilation.hh"
#include "AnnihilationStepLimiter.hh"

#include "DMParticleChi.hh"
#include "DMParticleChiScalar.hh"
#include "DMParticleChi1.hh"
#include "DMParticleChi2.hh"

#include "DMG4/DMParticles/DMParticleAPrime.hh"
#include "DMG4/DMParticles/DMParticleZPrime.hh"
#include "DMG4/DMParticles/DMParticleALP.hh"
#include "DMG4/DMParticles/DMParticleScalar.hh"
#include "DMG4/DMParticles/DMParticlePseudoScalar.hh"
#include "DMG4/DMParticles/DMParticleAxial.hh"

Geant4DarkMatterPhysics::Geant4DarkMatterPhysics(Geant4Context* ctxt, const std::string& name) :
    Geant4PhysicsList(ctxt, name)
{
    declareProperty("DMProcessType", DMProcessType = 0);
    declareProperty("EThresh", EThresh = 0);
    declareProperty("DMMass", DMMass = 0);
    declareProperty("Epsilon", Epsilon = 0);
    declareProperty("ANucl", ANucl = 0);
    declareProperty("ZNucl", ZNucl = 0);
    declareProperty("Density", Density = 0);
    declareProperty("DecayType", DecayType = 0);
    declareProperty("RDM", RDM = 1. / 3);
    declareProperty("fFactor", fFactor = 0.1);
    declareProperty("BranchingType", BranchingType = 0);
    declareProperty("AlphaD", AlphaD = 0.5);
    declareProperty("BiasSigmaFactor0", BiasSigmaFactor0 = 1);
    declareProperty("AnnihilationStepLimiterFactor", AnnilStepLimiterFactor = 5.);
    declareProperty("VerboseLevel", m_verbosity = 0);

    EThresh /= GeV;
    DMMass /= GeV;
    Density /= (g / cm3);

    switch(DMProcessType)
    {
    case 1:
        myDarkMatter = new DarkPhotons(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType);
        break;
    case 2:
        myDarkMatter = new DarkScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType);
        break;
    case 3:
        myDarkMatter = new DarkAxials(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
        break;
    case 4:
        myDarkMatter = new DarkPseudoScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
        break;
    case 5:
        if(DecayType) throw "DarkMassSpin2 with decays is not yet implemented";
        myDarkMatter = new DarkMassSpin2(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
        break;
    case 21:
        myDarkMatter = new ALP(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
        break;
    case 31:
        myDarkMatter = new DarkZ(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
        break;
    case 32:
        if(DecayType) throw "DarkMuPhilicScalar with decays is not yet implemented";
        myDarkMatter = new DarkMuPhilicScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType);
        break;
    case 34:
        if(DecayType) throw "DarkMuPhilicPseudoScalar with decays is not yet implemented";
        myDarkMatter = new DarkMuPhilicPseudoScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
        break;
    case 11:
        myDarkMatter = new DarkPhotonsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, 
                                                   DecayType, RDM, AlphaD, BranchingType, fFactor);
    break;
    case 12:
        myDarkMatter = new DarkScalarsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon,
                                                   DecayType, RDM, AlphaD, BranchingType, fFactor );
    break;
    case 13:
        myDarkMatter = new DarkAxialsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon,
                                                  DecayType, RDM, AlphaD, BranchingType, fFactor );
    break;
    case 14:
        myDarkMatter = new DarkPseudoScalarsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon,
                                                         DecayType, RDM, AlphaD, BranchingType, fFactor );
        break;
    case 15:
        myDarkMatter = new DarkMassSpin2Annihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon,
                                                     DecayType, RDM, AlphaD, BranchingType, fFactor );
        break;
    case 16:
        myDarkMatter = new DarkZAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon,
                                             DecayType, RDM, AlphaD, BranchingType,fFactor);
        break;
    default:
        throw "Wrong DM process type specified";
    }

    BiasSigmaFactor = BiasSigmaFactor0
    * (myDarkMatter->GetepsilBench() * myDarkMatter->GetepsilBench())
    / (myDarkMatter->Getepsil() * myDarkMatter->Getepsil());

    //For the e+ e- --> Z' --> ff process, we compute the cross section using epsil, so the code above has to be changed
    if (DMProcessType == 16) BiasSigmaFactor = BiasSigmaFactor0;
}


void Geant4DarkMatterPhysics::constructParticles(G4VUserPhysicsList* physics)
{
    Geant4PhysicsList::constructParticles(physics);

    switch(DMProcessType)
    {
    case 1:  //dark-photon bremmstrahlung
        DMParticleAPrime::Definition();
        if (DecayType != 0)
        {
            if (BranchingType >= 2)
            {
                DMParticleChi1::Definition();
                DMParticleChi2::Definition();
            }
        }
        break;
    case 2:
        DMParticleScalar::Definition();
        break;
    case 3:
        DMParticleAxial::Definition();
        break;
    case 4:
        DMParticlePseudoScalar::Definition();
        break;
    case 5:
        DMParticleAPrime::Definition();
        break;
    case 11: //annihilation processes
        DMParticleAPrime::Definition();
        break;
    case 12:
        DMParticleScalar::Definition();
        break;
    case 13:
        DMParticleAxial::Definition();
        break;
    case 14:
        DMParticlePseudoScalar::Definition();
        break;
    case 15:
        DMParticleAPrime::Definition(); // A' for the moment, the spin 2 particle not yet implemented
        if (DecayType != 0)
        {
            if ((BranchingType == 0) || (BranchingType == 1))
            {
                DMParticleChi::Definition();
            } else {
                if (BranchingType == 3)
                {
                    throw "Several decay channels according to BranchingType are not allowed in annihilation";
                }
                DMParticleChi1::Definition();
                DMParticleChi2::Definition();
            }
        }
        break;
    case 16:
    DMParticleZPrime::Definition();
    if (DecayType != 0)
    {
        if ((BranchingType == 0) || (BranchingType == 10))
        {
            // neutrinos final state
        }
        else if ((BranchingType == 1) || (BranchingType == 11))
        {
            // DM final state
            DMParticleChiScalar::Definition();
        }
        else throw "BranchingType not implemented";
    }
    break;
    case 21:
        DMParticleALP::Definition();
        break;
    case 31:
        DMParticleZPrime::Definition();
        break;
    case 32:
        DMParticleZPrime::Definition();
        break;
    case 34:
        DMParticleZPrime::Definition();
    }
}


void Geant4DarkMatterPhysics::constructProcesses(G4VUserPhysicsList* physics_list)
{
    Geant4PhysicsList::constructProcesses(physics_list);

    G4ParticleDefinition* theDMParticlePtr = nullptr;
    if (myDarkMatter->GetParentPDGID() == 11)
    {
        if (myDarkMatter->GetDMType() == 1) theDMParticlePtr = DMParticleAPrime::Definition();

        if (myDarkMatter->GetDMType() == 2) theDMParticlePtr = DMParticleScalar::Definition();

        if (myDarkMatter->GetDMType() == 3) theDMParticlePtr = DMParticleAxial::Definition();

        if (myDarkMatter->GetDMType() == 4) theDMParticlePtr = DMParticlePseudoScalar::Definition();

        // A' for the moment, the spin 2 particle not yet implemented
        if (myDarkMatter->GetDMType() == 5) theDMParticlePtr = DMParticleAPrime::Definition();
    }
    else if (myDarkMatter->GetParentPDGID() == -11)
    {
        // Annihilation
        if (myDarkMatter->GetDMType() == 1) theDMParticlePtr = DMParticleAPrime::Definition();

        if (myDarkMatter->GetDMType() == 2) theDMParticlePtr = DMParticleScalar::Definition();

        if (myDarkMatter->GetDMType() == 3) theDMParticlePtr = DMParticleAxial::Definition();

        if (myDarkMatter->GetDMType() == 4) theDMParticlePtr = DMParticlePseudoScalar::Definition();

        // Annihilation through spin 2 DM, A' for the moment
        if (myDarkMatter->GetDMType() == 5) theDMParticlePtr = DMParticleAPrime::Definition();

        // Annihilation through Z' (Lmu-Ltau or B-L models)
        if(myDarkMatter->GetDMType() == 11) theDMParticlePtr = DMParticleZPrime::Definition();
    }
    // Always Z' for the moment, scalar etc. particles from muons not yet implemented
    else if (myDarkMatter->GetParentPDGID() == 13) theDMParticlePtr = DMParticleZPrime::Definition(); 

    else if (myDarkMatter->GetParentPDGID() == 22) theDMParticlePtr = DMParticleALP::Definition();

    if (theDMParticlePtr == nullptr) throw "ConstructProcess: cannot determine the DM particle type";

    myDarkMatter->SetMA(theDMParticlePtr->GetPDGMass() / GeV);
    myDarkMatter->SetDMPDGID(theDMParticlePtr->GetPDGEncoding());
    myDarkMatter->PrepareTable();

    G4PhysicsListHelper* phLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
    phLHelper->DumpOrdingParameterTable();

    if (myDarkMatter->GetParentPDGID() == 11)
    {
        DMProcessDMBrem* DMBremPointer = new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor);

        G4ProcessManager* processManager = (G4Electron::ElectronDefinition())->GetProcessManager();
        processManager->AddDiscreteProcess(DMBremPointer);
        processManager = (G4Positron::PositronDefinition())->GetProcessManager();
        processManager->AddDiscreteProcess(DMBremPointer);
    }
    else if (myDarkMatter->GetParentPDGID() == -11)
    {
        DarkMatterAnnihilation *dmAnnihil=dynamic_cast<DarkMatterAnnihilation*>(myDarkMatter);
        AnnihilationStepLimiter *dmLimiterProc=new AnnihilationStepLimiter(dmAnnihil, "StepLimiterAnnihilation");
        dmLimiterProc->SetFactor(AnnilStepLimiterFactor);

        DMProcessAnnihilation* dmAnnihilProc = new DMProcessAnnihilation(dmAnnihil, theDMParticlePtr, 
                                                                         BiasSigmaFactor, dmLimiterProc);

        G4ParticleDefinition* posi=G4Positron::Definition();
        G4ProcessManager* processManager=posi->GetProcessManager();
        processManager->AddDiscreteProcess(dmLimiterProc);
        processManager->AddDiscreteProcess(dmAnnihilProc);
    }
    else if(myDarkMatter->GetParentPDGID() == 13)
    {
        DMProcessDMBrem* DMBremPointer = new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor);

        G4ProcessManager* processManager = (G4MuonMinus::MuonMinusDefinition())->GetProcessManager();
        processManager->AddDiscreteProcess(DMBremPointer);
        processManager = (G4MuonPlus::MuonPlusDefinition())->GetProcessManager();
        processManager->AddDiscreteProcess(DMBremPointer);
    }
    else if(myDarkMatter->GetParentPDGID() == 22)
    {
        DMProcessPrimakoffALP* DMPrimakoffALPPointer = new DMProcessPrimakoffALP(myDarkMatter, theDMParticlePtr, BiasSigmaFactor);

        G4ProcessManager* processManager = (G4Gamma::GammaDefinition())->GetProcessManager();
        processManager->AddDiscreteProcess(DMPrimakoffALPPointer);
    }
}

