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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "../include/PrimaryGeneratorAction.hh"

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "g4root.hh"
#include "G4UnitsTable.hh"
#include "G4GenericMessenger.hh"

#define EPS 0.1

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr),
  fEnvelope(nullptr),
  fTagger(nullptr),
  fMessenger(nullptr)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete fParticleGun;
  // delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.

  if (!fEnvelope)
  {
    G4LogicalVolume* envLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelope = dynamic_cast<G4Box*>(envLV->GetSolid());

    G4LogicalVolume* tagLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Tagger LV");
    if ( tagLV ) fTagger = dynamic_cast<G4Box*>(tagLV->GetSolid());
  }

  ParticleKinematicsGenerator();

  fParticleGun->GeneratePrimaryVertex(anEvent);

}

void PrimaryGeneratorAction::ParticleKinematicsGenerator(){

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(21.*MeV);
  
  const G4double initial_z = fEnvelope->GetZHalfLength();
  G4double max_x = (fTagger->GetXHalfLength()), max_y = (fTagger->GetYHalfLength());

  // direction of the beam
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));
 
  // set position of the particle

  if(randomic_generation)
  {
    fInitial_X = (G4UniformRand()-0.5)*(2-EPS);
    fInitial_Y = (G4UniformRand()-0.5)*(2-EPS);
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(max_x*fInitial_X/100., max_y*fInitial_Y/100., initial_z));

}

void PrimaryGeneratorAction::DefineCommands() {
  
  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/PadmETag/Generation/", "Primary generator");

  // position Y
  auto& setRandomGeneration
  = fMessenger->DeclareProperty("SetRandomGeneration",randomic_generation);
  setRandomGeneration.SetParameterName("RandomGenerationBool", true);
  setRandomGeneration.SetCandidates("0 1");
  setRandomGeneration.SetDefaultValue("0");

  // position X
  auto& setXcommand
  = fMessenger->DeclareProperty("SetX",fInitial_X);
  setXcommand.SetParameterName("percentageX", true);
  setXcommand.SetRange("percentageX>=-100. && percentageX<100.");
  setXcommand.SetDefaultValue("0.");

  // position Y
  auto& setYcommand
  = fMessenger->DeclareProperty("SetY",fInitial_Y);
  setYcommand.SetParameterName("percentageY", true);
  setYcommand.SetRange("percentageY>=-100. && percentageY<100.");
  setYcommand.SetDefaultValue("0.");

  // positions X and Y
  auto& setXYcommand
    = fMessenger->DeclareMethod("SetXY",
                                &PrimaryGeneratorAction::SetIncomingPositionsXY, 
                                "Set (X,Y) position of incoming beam");
  setXYcommand.SetParameterName("position", true);
  setXYcommand.SetRange("percentageX>=-100. && percentageX<100. && percentageY>=-100. && percentageY<100.");
  setXYcommand.SetDefaultValue("0., 0.");

}