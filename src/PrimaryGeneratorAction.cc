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

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr), 
  fEnvelope(nullptr),
  fTagger(nullptr)
{
 
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event

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

// this function sets kinematic and specifics of the incoming ray
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
  fParticleGun->SetParticlePosition(G4ThreeVector(-0.63*max_x, -0.32*max_y, initial_z));

}
