# Simulation of electron tagger in the PADME experiment.
PadmeExperiment http://padme.lnf.infn.it/experiment/
@LNF Laboratori Nazionali di Frascati https://w3.lnf.infn.it/

Summerstudent INFN: Emanuele Rosi > rosiemanuele99@gmail.com
Supervisor LNF: Emanuele Leonardi > ************************

## Description

The electron/positron tagger has the purpose to detect these particles through scintillation mechanism and, using the time differences and distributions of the signals on 8 silicon photomultipliers, reconstruct the incidence position of the incoming particle.

# Geometry
    
    (1) The world and the envelope are boxes filled with low vacuum atmospheric air.
        Some optical parameters as the refractive index are given.

    (2) The main componenti is the Tagger: it is a plastic scintillator by Saint Gobain;
        its features are given by datasheets and properly implemented. Its surface is entirely polished
        and it is covered with a high reflective optical paint.

    (3) At the end of the tagger, 8 silicon photomultipliers (SiPM) by Hamamatsu are fixed through optical glue.
        A SiPM is schematized as a little 3x3 mm^2 layer of pyrex glass. Its purpose is to detect optical photons.
        Again, optical properties of the Pyrex window are given. Moreover, the quantum efficiency of the SiPM in implemented in detection time in SteppingAction.cc.

# Particle generation

    The simulations generates electrons and/or positrons. The initial energy of an electron is 21. MeV and its direction is 
    perpendicular to one face of the tagger; the position (x,y) is a free choice of the user.

# Photon detection strategy

    We state that a photon is detected if:
    - it is an Optical Photon and its wavelength is in the range [300,900] nm, that's the detection range of a SiPM.
    - if (statistically) it is engaged by the SiPM (quantum efficiency).
    - if it is not absorbed by the medium or the surfaces.
    - if its PreStepPoint is in the tagger volume and its PostStepPoint in the SiPM volume. In order to avoid that a photon
      cross the SiPM w/o any step point, the StepLength is limited in SiPM volume.
    After that, the photon is killed (fStopAndKill).

# Data storage

    We open only one file named data_eTag<sign><number>.root, where <sign> and <number> indicate the percentage of the X position of the incoming electron.
    