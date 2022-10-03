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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "../include/DetectorConstruction.hh"
#include "../include/OutputColors.hh"

#include "G4RotationMatrix.hh"

#include "G4RunManager.hh" 
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

const double hPlanck = 4.135655e-15;
const double c_light = 3e+8;
const double meters_to_nanometers = 1e9;

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fBGOcrystal(nullptr),
  fPlasticScintillator_1(nullptr),
  fPlasticScintillator_2(nullptr),
  fCerenkovPMT(nullptr),
  fScintillatorPMT(nullptr),
  fStepLimit(nullptr),
  fScoringVolume(0)
{}

DetectorConstruction::~DetectorConstruction(){
  delete fStepLimit;
}

G4VPhysicalVolume* DetectorConstruction::Construct(){ 
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Dimensions of BGO and plastic scintillators
  const G4double shape_bgoXZ = 2.2*cm, shape_bgoY = 18*cm;
  const G4double shape_plasticX = 5*cm, shape_plasticY = 1*cm, shape_plasticZ = 10*cm;
  const G4double distance_scintillators = 0.5*cm;
  const G4double distance_BGOscintillators = 10*cm;

  // Dimensions of BGO photomultipliers
  const G4double shape_PMT_XZ = shape_bgoXZ;
  const G4double shape_PMT_Y = 0.5*cm; // the order of absorption length in BG crystals
  
  // minimal radius such that the experiment is circumscribed in a sphere
  // we willl use minimal radius + 40%
  G4double Yheight = 2*shape_PMT_Y + shape_bgoY;
  G4double H = 2 * shape_plasticZ + distance_BGOscintillators + distance_scintillators;
  minimal_radius = (4 * H*H + Yheight*Yheight) / (8*H);
  G4double radius_sphere = minimal_radius * 1.4; // < sqrt(3) = 1.732..
  
  std::cout << OBOLDWHITE
    << "Minimal radius of the envelope sphere:\t" << G4BestUnit(minimal_radius,"Length") << "\n"
    << "Effective radius of the envelope sphere:\t" << G4BestUnit(radius_sphere,"Length")
    << ORESET << std::endl;
  
  // Envelope material
  G4Material* envelope_material = CreateOpticalAir();
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4Material* world_material = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Sphere* solidWorld = new G4Sphere("World", 0, sqrt(3)*radius_sphere, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_material,      //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                   
  //     
  // Envelope
  //  
  G4Sphere* solidEnv = new G4Sphere("World", 0, radius_sphere, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        envelope_material,   //its material
                        "Envelope");         //its name

  G4VPhysicalVolume* physEnvelope = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  // limit stepLength to 4.8cm in order to optimize rejecting logic of non-detected events 
  // (see SteppingAction.cc SteppingAction::EstinguishParticleIfNotTrigger )
  G4double maxStep = 4.8*cm;
  fStepLimit = new G4UserLimits(maxStep);
  logicEnv->SetUserLimits(fStepLimit);

  //     
  // BGO Crystal + Photomultipliers
  //  
  G4Material* bgo_material = this->CreateBismuthGermaniumOxygen();
  G4Material* borosilicate = this->CreatePyrex();

  // position: choosen in order to minimize the envelope sphere
  G4double pos_BGO = minimal_radius - 2 * shape_plasticZ - distance_BGOscintillators - distance_scintillators - shape_bgoXZ / 2;
  G4ThreeVector bgo_position = G4ThreeVector(0, 0, pos_BGO);
  G4ThreeVector pmt1_position = G4ThreeVector(0, 0.5*(shape_bgoY+shape_PMT_Y), pos_BGO);
  G4ThreeVector pmt2_position = G4ThreeVector(0, -0.5*(shape_bgoY+shape_PMT_Y), pos_BGO);
        
  // BGO & PMTs shape
  G4Box* BGOShape = new G4Box("BGO Box", 0.5*shape_bgoXZ, 0.5*shape_bgoY, 0.5*shape_bgoXZ);
  G4Box* PMTsShape = new G4Box("PMT Box", 0.5*shape_PMT_XZ, 0.5*shape_PMT_Y, 0.5*shape_PMT_XZ);

  // logical volumes
  G4LogicalVolume* BGO_LogicalVolume = 
    new G4LogicalVolume(BGOShape,                 //its solid
                        bgo_material,             //its material
                        "BGO Crystal");           //its name
  
  G4LogicalVolume* PMT_LogicalVolume = 
    new G4LogicalVolume(PMTsShape,                               //its solid
                        borosilicate,                            //its material
                        "PMT made up of borosilicate glass");    //its name

  // in order to ensure at least one step in the PMT
  PMT_LogicalVolume->SetUserLimits(new G4UserLimits(shape_PMT_Y/2));

  // delta vector useful to translate PMTs
  G4ThreeVector original_vector = G4ThreeVector(0, 0.5*(shape_bgoY+shape_PMT_Y), 0);
  G4ThreeVector translated_vector = original_vector;

  G4ThreeVector delta_vect = G4ThreeVector();
  delta_vect = (translated_vector - original_vector);

  fBGOcrystal = new G4PVPlacement(0,
                    bgo_position,                 //at position
                    BGO_LogicalVolume,            //its logical volume
                    "BGO Crystal",                //its name
                    logicEnv,                     //its mother volume
                    false,                        //no boolean operation
                    0,                            //copy number
                    checkOverlaps);               //overlaps checking
  
  fCerenkovPMT = new G4PVPlacement(0,
                    pmt1_position + delta_vect,   //at position
                    PMT_LogicalVolume,            //its logical volume
                    "Cherenkov PMT",              //its name
                    logicEnv,                     //its mother volume
                    false,                        //no boolean operation
                    0,                            //copy number
                    checkOverlaps);               //overlaps checking

  fScintillatorPMT = new G4PVPlacement(0,
                    pmt2_position - delta_vect,   //at position
                    PMT_LogicalVolume,            //its logical volume
                    "Scintill. PMT",              //its name
                    logicEnv,                     //its mother volume
                    false,                        //no boolean operation
                    0,                            //copy number
                    checkOverlaps);               //overlaps checking

  // optical properties of BGO
  OpticalSurfaceBGO(fBGOcrystal, physEnvelope);
  OpticalSurfaceBGO_PMT(fBGOcrystal, fCerenkovPMT);
  OpticalSurfaceBGO_PMT(fBGOcrystal, fScintillatorPMT);

  //
  // Plastic scintillator
  //
  // the position in choosen in order to minimize the envelope sphere
  G4double pos_scint1 = minimal_radius - distance_scintillators - 1.5 * shape_plasticZ;
  G4double pos_scint2 = minimal_radius - shape_plasticZ / 2;
  G4Material* plastic_material = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  // Plastic scintillators shape
  G4Box* PlasticShape = new G4Box("BGO Box", 0.5*shape_plasticX, 0.5*shape_plasticY, 0.5*shape_plasticZ);

  // Plastic scintillators positions
  G4ThreeVector trigger_1_position = G4ThreeVector(0, 0, pos_scint1);
  G4ThreeVector trigger_2_position = G4ThreeVector(0, 0, pos_scint2);

  G4LogicalVolume* Plastic_LogicalVolume =                         
    new G4LogicalVolume(PlasticShape,              //its solid
                        plastic_material,          //its material
                        "Plastic Scintillator");   //its name
               
  fPlasticScintillator_1 = new G4PVPlacement(0,    //no rotation
                    trigger_1_position,            //at position
                    Plastic_LogicalVolume,         //its logical volume
                    "Plastic Scintillator",        //its name
                    logicEnv,                      //its mother volume
                    false,                         //no boolean operation
                    0,                             //copy number
                    checkOverlaps);                //overlaps checking

  fPlasticScintillator_2 = new G4PVPlacement(0,    //no rotation
                    trigger_2_position,            //at position
                    Plastic_LogicalVolume,         //its logical volume
                    "Plastic Scintillator",        //its name
                    logicEnv,                      //its mother volume
                    false,                         //no boolean operation
                    0,                             //copy number
                    checkOverlaps);                //overlaps checking

  fScoringVolume = logicEnv;

  //always return the physical World
  return physWorld;
}

// optical surface of BGO
void DetectorConstruction::OpticalSurfaceBGO(G4VPhysicalVolume* BGO_PV, G4VPhysicalVolume* Envelope_PV) const {

  G4OpticalSurface* opBGOSurface = new G4OpticalSurface("BGO Surface");
  // see here: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/BackupVersions/V10.7/html/TrackingAndPhysics/physicsProcess.html
  opBGOSurface->SetType(dielectric_LUTDAVIS); // dielectric_LUTDAVIS->get table BGO, dielectric_dielectric
  opBGOSurface->SetModel(DAVIS);  // fig 1.7 (link above) or glisur, DAVIS, unified
  opBGOSurface->SetFinish(Rough_LUT); // Rough_LUT, Polished_LUT, ground

  G4LogicalBorderSurface* LogicalBGOSurface = new G4LogicalBorderSurface(
    "BGO Surface", BGO_PV, Envelope_PV, opBGOSurface);
  
  G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
    LogicalBGOSurface->GetSurface(BGO_PV, Envelope_PV)->GetSurfaceProperty());

  if(opticalSurface)
    opticalSurface->DumpInfo();

  // Generate & Add Material Properties Table attached to the optical surface
  // energy = hPlanck * (light speed) / wavelength
  // reference for reflectivity https://aip.scitation.org/doi/10.1063/1.3272909
  // reference for transmittance https://www.researchgate.net/publication/264828576_The_Radiation_Hard_BGO_Crystals_for_Astrophysics_Applications
  G4int n10 = 10;
  G4int n6 = 6;
  G4double energies_photons[] = { 2.*eV, 2.5*eV, 3.*eV, 3.5*eV, 4.*eV,
                                4.5*eV, 4.75*eV, 4.9*eV, 5.*eV, 5.4*eV };
  G4double reflectivity[] = { 0.125, 0.13, 0.14, 0.155, 0.175,
                              0.25, 0.29, 0.243, 0.21, 0.22 };
  G4double energies_photons_bis[] = { hPlanck*c_light*meters_to_nanometers/310.*eV,
                                      hPlanck*c_light*meters_to_nanometers/350.*eV,
                                      hPlanck*c_light*meters_to_nanometers/400.*eV,
                                      hPlanck*c_light*meters_to_nanometers/500.*eV,
                                      hPlanck*c_light*meters_to_nanometers/600.*eV,
                                      hPlanck*c_light*meters_to_nanometers/700.*eV };                      
  //G4double transm[] = { 0.1, 0.70, 0.78, 0.79, 0.80, 0.80 };
  G4double transm[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };

  G4MaterialPropertiesTable* SurfaceTable = new G4MaterialPropertiesTable();

  SurfaceTable->AddProperty("REFLECTIVITY", energies_photons, reflectivity, n10);
  SurfaceTable->AddProperty("TRANSMITTANCE", energies_photons_bis, transm, n6);
  
  SurfaceTable->DumpTable();
  opBGOSurface->SetMaterialPropertiesTable(SurfaceTable);

}

// optical between of BGO and PMTs
void DetectorConstruction::OpticalSurfaceBGO_PMT(G4VPhysicalVolume* BGO_PV, G4VPhysicalVolume* TheOtherPV) const {

  G4OpticalSurface* opBGOSurface = new G4OpticalSurface("BGO Surface");
  // see here: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/BackupVersions/V10.7/html/TrackingAndPhysics/physicsProcess.html
  opBGOSurface->SetType(dielectric_LUTDAVIS);
  opBGOSurface->SetModel(DAVIS);
  opBGOSurface->SetFinish(PolishedESRGrease_LUT); // Using optical grease between PMTs and BGO

  G4LogicalBorderSurface* LogicalBGOSurface = new G4LogicalBorderSurface(
    "BGO Surface", BGO_PV, TheOtherPV, opBGOSurface);
  
  G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
    LogicalBGOSurface->GetSurface(BGO_PV, TheOtherPV)->GetSurfaceProperty());

  if(opticalSurface) opticalSurface->DumpInfo();

  // Generate & Add Material Properties Table attached to the optical surface
  // energy = hPlanck * (light speed) / wavelength
  // reference for reflectivity https://aip.scitation.org/doi/10.1063/1.3272909
  // reference for transmittance https://www.researchgate.net/publication/264828576_The_Radiation_Hard_BGO_Crystals_for_Astrophysics_Applications
  G4int n10 = 10;
  // G4int n6 = 6;
  G4double energies_photons[] = { 2.*eV, 2.5*eV, 3.*eV, 3.5*eV, 4.*eV,
                                4.5*eV, 4.75*eV, 4.9*eV, 5.*eV, 5.4*eV };
  G4double reflectivity[] = { 0.125, 0.13, 0.14, 0.155, 0.175,
                              0.25, 0.29, 0.243, 0.21, 0.22 };

  G4MaterialPropertiesTable* SurfaceTable = new G4MaterialPropertiesTable();

  SurfaceTable->AddProperty("REFLECTIVITY", energies_photons, reflectivity, n10);
  
  SurfaceTable->DumpTable();
  opBGOSurface->SetMaterialPropertiesTable(SurfaceTable);

}

// BGO - material
G4Material* DetectorConstruction::CreateBismuthGermaniumOxygen() const {
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // BGO material definition
  // references: Bi4-Ge3-O12
  // https://iopscience.iop.org/article/10.1088/1361-6560/aa6a49/pdf
  // https://www.gammadata.se/assets/Uploads/BGO-data-sheet.pdf

  // Composition
  G4Material* bgo_basic = nist->FindOrBuildMaterial("G4_BGO");
  G4Material* bgo_material = new G4Material("BismuthGermaniumOxygen Crystal", 7.13*g/cm3, bgo_basic);

  // Optical properties
  const G4int n10 = 10;

  // energy = hPlanck * (light speed) / wavelength
  G4double energies_photons[n10] = {hPlanck*c_light*meters_to_nanometers/300.*eV,
                                    hPlanck*c_light*meters_to_nanometers/350.*eV,
                                    hPlanck*c_light*meters_to_nanometers/400.*eV,
                                    hPlanck*c_light*meters_to_nanometers/450.*eV,
                                    hPlanck*c_light*meters_to_nanometers/500.*eV,
                                    hPlanck*c_light*meters_to_nanometers/550.*eV,
                                    hPlanck*c_light*meters_to_nanometers/600.*eV,
                                    hPlanck*c_light*meters_to_nanometers/650.*eV,
                                    hPlanck*c_light*meters_to_nanometers/700.*eV,
                                    hPlanck*c_light*meters_to_nanometers/750.*eV,};

  // from BGO FermiLab .pptx
  G4double rindex[n10] = {2.75, 2.42, 2.30, 2.25, 2.21, 2.17, 2.16, 2.15, 2.15, 2.15};
  G4double absorption[n10] = {0.1*mm, 75.*mm, 82.*mm, 90.*mm, 100.*mm, 105.*mm, 112.*mm, 120.*mm, 125.*mm, 130.*mm};
  G4double scintillation_spectrum[n10] = {0, 0.1/100, 3.4/100, 11.5/100, 14.5/100, 10.8/100, 5.5/100, 3/100, 1.7/100, 0.7/100}; // percent

  // new instance of Material Properties
  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
 
  // properties independent of energy
  MPT->AddConstProperty("FASTTIMECONSTANT", 300.*ns);
  MPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  MPT->AddConstProperty("SCINTILLATIONYIELD", 8200./MeV);

  // properties that depend on energy
  MPT->AddProperty("RINDEX", energies_photons, rindex, n10)->SetSpline(true);
  MPT->AddProperty("ABSLENGTH", energies_photons, absorption, n10);
  MPT->AddProperty("FASTCOMPONENT", energies_photons, scintillation_spectrum, n10)->SetSpline(true);

  // bgo material
  bgo_material->SetMaterialPropertiesTable(MPT);
  bgo_material->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  return bgo_material;
}

// Air - material
G4Material* DetectorConstruction::CreateOpticalAir() const {
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Air material definition
  G4Material* air_basic = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* air_optical = new G4Material("Air", 1.204*kg/m3, air_basic);
  
  // otpical properties
  const G4int n3 = 3;
  // energy = hPlanck * (light speed) / wavelength
  G4double photonenergy[n3] = {hPlanck*c_light*meters_to_nanometers/320.*eV,  // lower wavelength cutoff 320 nm
                            hPlanck*c_light*meters_to_nanometers/400.*eV,     // intermediate energy
                            hPlanck*c_light*meters_to_nanometers/480.*eV};    // maximum emission at 480 nm
  G4double rindex[n3]     = {1.000293, 1.000293, 1.000293};                   // google
  // new instance of Material Properties
  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  // properties that depend on energy
  MPT->AddProperty("RINDEX", photonenergy, rindex, n3);

  // material
  air_optical->SetMaterialPropertiesTable(MPT);

  return air_optical;
}

// Birosilicate glass - material
G4Material* DetectorConstruction::CreatePyrex() const {
  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // BGO material definition
  G4Material* pyrex_basic = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  G4Material* pyrex = new G4Material("Borosilicate_Glass", 2.23*g/cm3, pyrex_basic);
  
  const G4int n3 = 3;
  // energy = hPlanck * (light speed) / wavelength
  G4double photonenergy[n3] = {hPlanck*c_light*meters_to_nanometers/320.*eV,  // lower wavelength cutoff 320 nm
                            hPlanck*c_light*meters_to_nanometers/400.*eV,     // intermediate energy
                            hPlanck*c_light*meters_to_nanometers/480.*eV};    // maximum emission at 480 nm
  G4double rindex[n3]     = {1.471, 	1.471, 	1.471};                   // google
  // new instance of Material Properties
  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();

  // properties that depend on energy
  MPT->AddProperty("RINDEX", photonenergy, rindex, n3);

  // material
  pyrex->SetMaterialPropertiesTable(MPT);

  return pyrex;
}