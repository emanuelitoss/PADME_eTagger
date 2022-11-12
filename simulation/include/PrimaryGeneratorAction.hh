//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <fstream>
using namespace std;

#include "../include/OutputColors.hh"
#include "../include/DetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4Box.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4Sphere;
class G4GenericMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{

  public:

    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event*);
    
    // particle kinematic
    void ParticleKinematicsGenerator();

    // getters
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
    G4double GetInitialX() const { return fInitial_X; }
    G4double GetInitialY() const { return fInitial_Y; }

  private:

    void DefineCommands();

    void SetIncomingPositionX(G4double x);
    void SetIncomingPositionY(G4double y);
    void SetIncomingPositionsXY(G4double x, G4double y);
    
    G4ParticleGun* fParticleGun;

    G4Box* fEnvelope;
    G4Box* fTagger;

    G4GenericMessenger* fMessenger;
    G4double fInitial_X = 0.;
    G4double fInitial_Y = 0.;

    G4bool randomic_generation = 0.;

};

inline void PrimaryGeneratorAction::SetIncomingPositionX(G4double x) { fInitial_X = x; }
 
inline void PrimaryGeneratorAction::SetIncomingPositionY(G4double y) { fInitial_Y = y; }

inline void PrimaryGeneratorAction::SetIncomingPositionsXY(G4double x, G4double y){
  SetIncomingPositionX(x);
  SetIncomingPositionY(y);
}

#endif