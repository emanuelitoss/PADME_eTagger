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
#include "G4Sphere.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "g4root.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String outputName)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0), 
  fEnvelopeSphere(0)
{

  fUseRandomicGeneration = false;

  this->setOutput(outputName);
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
  
  G4double envSizeR = 0;

  if (!fEnvelopeSphere)
  {
    G4LogicalVolume* envLV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if ( envLV ) fEnvelopeSphere = dynamic_cast<G4Sphere*>(envLV->GetSolid());
  }

  if ( fEnvelopeSphere ) {
    // in order to correct the envelope radius (+40% with respect minimal radius) we divide by 140%.
    envSizeR = fEnvelopeSphere->GetOuterRadius()/1.4;
  }
  else {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()", "MyCode0002",JustWarning,msg);
  }

  ParticleKinematicsGenerator();

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

// this function sets kinematic and specifics of the incoming ray
void PrimaryGeneratorAction::ParticleKinematicsGenerator(){

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="mu-");
  fParticleGun->SetParticleDefinition(particle);
  // typical energy of a cosmic muon at sea level (see PDG reference)
  fParticleGun->SetParticleEnergy(3.*GeV);

  // generation of randomic angles
  double phi = G4UniformRand() * 2 * M_PI * radian;
  double theta;
  do{
    theta = acos(pow(G4UniformRand(),1./3));
  } while(theta > 0.6*radian); // for theta > 0.6 rad we don't detect any muon, so I generate again.

  // in order to correct the envelope radius (+40% with respect minimal radius)
  //const double radius = fEnvelopeSphere->GetOuterRadius();///1.4;
  //  theta = acos(pow(G4UniformRand(),1./3));
  //}while(theta > 0.7*radian); // contraint based on angular analysis
  
  const double radius = fEnvelopeSphere->GetOuterRadius()/1.4;

  // direction of the beam
  G4ThreeVector* Direction_Beam = new G4ThreeVector(0, 0, -radius);
  Direction_Beam->rotateY(theta);
  Direction_Beam->rotateZ(phi);
  if(fUseRandomicGeneration) fParticleGun->SetParticleMomentumDirection(*Direction_Beam);
  else fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));
  
  // tangent plane position generation
  double position_x;
  double position_y;

  position_x = (G4UniformRand() - 0.5) * 2 * radius;
  position_y = (G4UniformRand() - 0.5) * 2 * radius;

  G4ThreeVector* Position_Beam = new G4ThreeVector(position_x, position_y, radius);
  Position_Beam->rotateY(theta);
  Position_Beam->rotateZ(phi);
 
  // set position of the particle
  if(fUseRandomicGeneration)  fParticleGun->SetParticlePosition(*Position_Beam);
  else  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,radius));

  //delete output;
  delete Position_Beam;
  delete Direction_Beam;

}
