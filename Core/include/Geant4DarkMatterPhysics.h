#pragma once

#include <DDG4/Geant4PhysicsList.h>
#include <DDG4/Factories.h>

#include <G4VUserPhysicsList.hh>

#include "DarkMatter.hh"

namespace dd4hep
{
    namespace sim
    {
        class Geant4DarkMatterPhysics : public Geant4PhysicsList
        {
        public:
            Geant4DarkMatterPhysics() = delete;
            Geant4DarkMatterPhysics(const Geant4DarkMatterPhysics&) = delete;
            Geant4DarkMatterPhysics(Geant4Context* ctxt, const std::string& name);
            virtual ~Geant4DarkMatterPhysics() = default;
            virtual void constructParticles(G4VUserPhysicsList* physics) override;
            virtual void constructProcesses(G4VUserPhysicsList* physics_list) override;

        private:
            G4int DMProcessType;
            G4double EThresh;
            double DMMass;
            double Epsilon;
            G4double ANucl;
            G4double ZNucl;
            G4double Density;
            G4int DecayType;
            G4double RDM;
            G4double fFactor;
            G4int BranchingType;
            G4double AlphaD;
            G4double BiasSigmaFactor0;
            G4double AnnilStepLimiterFactor;
            int m_verbosity;

            DarkMatter* myDarkMatter;
            G4double BiasSigmaFactor;
        };
    }
}

using namespace dd4hep::sim;
DECLARE_GEANT4ACTION(Geant4DarkMatterPhysics)

/*
- The definition for the steering file is:

  def setupDarkMatter(kernel):
    darkmatter = PhysicsList(kernel, 'Geant4DarkMatterPhysics/DMPhys')
    darkmatter.VerboseLevel = 2
    darkmatter.enableUI()

    seq = kernel.physicsList()
    seq.adopt(darkmatter)
  SIM.physics.setupUserPhysics(setupCerenkov)
 
 - reference DD4hep/DDG4/python/DDSim/Helper/Physics.py
 */
