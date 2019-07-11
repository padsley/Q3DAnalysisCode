#include "Kinematics.h"

//mass 3He  = 2809.413242
//mass d    = 1876.123740
//mass 30Si = 27920.38734
//mass 31P  = 28851.87339

//40./100./2.65*0.5/cos(TargetAngle*TMath::Pi()/180.)
//40./100./2.65*0.5/cos((Q3DAngle-TargetAngle)*TMath::Pi()/180.)
//20./100./2.253/cos((Q3DAngle-TargetAngle)*TMath::Pi()/180.)

double ExToBrho(double TBeam, double Q3DAngle, double TargetAngle, double TargetThickness, double CarbonThickness, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Ex)
{
    //TargetThickness is the thickness of the actual target, CarbonThickness is that of any backing
    
    //double DoEnergyLoss(double TInitial, double Thickness, double a, double b, double c) - just here for reference
    
    //Correct the beam energy for half of the target thickness - assume that the reactions happen in the centre of the target
    //Target is 40 ug/cm^2, density of SiO2 is 2.65 g/cm**3, 40/1e6 to go to g/cm^2, *1e4 to convert to um, so thickness =  areal density/100/density with areal density in ug/cm**2, density in g/cm^3 and thickness in um
    TBeam = DoEnergyLoss(TBeam,
                         TargetThickness*0.5/cos(TargetAngle*TMath::Pi()/180.),
                         0.320508,
                         5.034675,
                         -3.82995
    );
    
    double TEjectile = ComputeEjectileEnergy(TBeam,
                                             Q3DAngle,
                                             MassBeam,
                                             MassTarget,
                                             MassEjectile,
                                             MassRecoil,
                                             Ex
    );
    cout << "TEjectile: " << TEjectile << endl;
    
    //Take ejectile information and do the energy losses
    //Energy loss for half target
    TEjectile = DoEnergyLoss(TEjectile,
                             TargetThickness*0.5/cos((Q3DAngle-TargetAngle)*TMath::Pi()/180.),
                             1.786181,
                             27.04590,
                             -25.44124
    );
    cout << "TEjectile: " << TEjectile << endl;
    
    //Energy loss for carbon backing
    TEjectile = DoEnergyLoss(TEjectile,
                             CarbonThickness/cos((Q3DAngle-TargetAngle)*TMath::Pi()/180.),
                             1.96859,
                             26.78673,
                             -26.50826
    );
    cout << "TEjectile: " << TEjectile << endl;
    
    return ConvertKEToBrho(TEjectile, MassEjectile); //Convert the kinematic energy of the ejectile to a magnetic ridigity value
}

double BrhoToEx(double TBeam, double Q3DAngle, double TargetAngle, double TargetThickness, double CarbonThickness, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Brho)
{
    TBeam = DoEnergyLoss(TBeam,
                         TargetThickness*0.5/cos(TargetAngle*TMath::Pi()/180.),
                         0.320508,
                         5.034675,
                         -3.82995
    );
    
    double TEjectile = ConvertBrhoToKE(Brho, MassEjectile);
    
    //Correct for the energy loss in the carbon layer - NOTE THAT WE USE A NEGATIVE THICKNESS FOR THIS WHICH IS THE INVERSE PROCESS!
    TEjectile = DoEnergyLoss(TEjectile,
                             -CarbonThickness/cos((Q3DAngle-TargetAngle)*TMath::Pi()/180.),
                             1.96859,
                             26.78673,
                             -26.50826
    );
    
    //Correct for the half-target energy loss - AGAIN NEGATIVE THICKNESS
    TEjectile = DoEnergyLoss(TEjectile,
                             -TargetThickness*0.5/cos((Q3DAngle-TargetAngle)*TMath::Pi()/180.),
                             1.786181,
                             27.04590,
                             -25.44124
    );

    return ComputeExcitationEnergy(TBeam, MassBeam, MassTarget, MassEjectile, MassRecoil, TEjectile, Q3DAngle);
}

double ComputeEjectileEnergy(double TBeam, double Q3DAngle, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Ex)
{
    MassRecoil += Ex;
    
    double s = MassBeam*MassBeam + MassTarget*MassTarget + 2*MassTarget*(TBeam + MassBeam);//The Mandelstam 's' variable = the centre-of-mass energy squared
    TLorentzVector fTotalEnergyImpulsionCM = TLorentzVector(0,0,0,sqrt(s)); //4-vector giving the energy in the colliding system - this is conserved!
    
    //Compute centre-of-mass energies
    double ECM_1 = (s + MassBeam*MassBeam - MassTarget*MassTarget)/(2*sqrt(s));
    double ECM_2 = (s + MassTarget*MassTarget - MassBeam*MassBeam)/(2*sqrt(s));
    double ECM_3 = (s + MassEjectile*MassEjectile - MassRecoil*MassRecoil)/(2*sqrt(s));
    double ECM_4 = (s + MassRecoil*MassRecoil - MassEjectile*MassEjectile)/(2*sqrt(s));
    
    //Compute centre-of-mass momenta
    double pCM_1 = sqrt(ECM_1*ECM_1 - MassBeam*MassBeam);
    double pCM_2 = sqrt(ECM_2*ECM_2 - MassTarget*MassTarget);
    double pCM_3 = sqrt(ECM_3*ECM_3 - MassEjectile*MassEjectile);
    double pCM_4 = sqrt(ECM_4*ECM_4 - MassRecoil*MassRecoil);
    
    //Initial lab 3-momenta
    TVector3 fImpulsionLab_1 = TVector3(0,0,sqrt(TBeam*TBeam + 2*TBeam*MassBeam));
    TVector3 fImpulsionLab_2 = TVector3(0,0,0);
    
    //Initial lab 4-momenta
    TLorentzVector fEnergyImpulsionLab_1 = TLorentzVector(fImpulsionLab_1,MassBeam+TBeam);
    TLorentzVector fEnergyImpulsionLab_2 = TLorentzVector(fImpulsionLab_2,MassTarget);
    
    //Total lab 4-momentum
    TLorentzVector fTotalEnergyImpulsionLab = fEnergyImpulsionLab_1 + fEnergyImpulsionLab_2;
    
    //Boost between the lab and CoM frames
    double BetaCM = fTotalEnergyImpulsionLab.Beta();
    
    //Work out CoM 4-momentum of the beam
    TLorentzVector fEnergyImpulsionCM_1 = fEnergyImpulsionLab_1;
    fEnergyImpulsionCM_1.Boost(0,0,-BetaCM);
    
    //...and the target
    TLorentzVector fEnergyImpulsionCM_2 = fEnergyImpulsionLab_2;
    fEnergyImpulsionCM_2.Boost(0,0,-BetaCM);
    
    //Awkward poking around to get the centre-of-mass angle because we set the ***lab*** angle with the Q3D
    double thetaCM = 0;
    double rapidity = TMath::ATanH(BetaCM);
    if(pCM_3 - pCM_4 > 1e-6)cout << "I don't understand this: " << pCM_3 << "\t" << pCM_4 << endl;
    double p3 = (sqrt(MassEjectile*MassEjectile + pCM_3*pCM_3) * cos(Q3DAngle*TMath::Pi()/180.) * sinh(rapidity) + cosh(rapidity) * sqrt(pCM_3*pCM_3 - MassEjectile*MassEjectile*sin(Q3DAngle*TMath::Pi()/180.)*sin(Q3DAngle*TMath::Pi()/180.)*sinh(rapidity)*sinh(rapidity))) / (1.0 + sin(Q3DAngle*TMath::Pi()/180.)*sin(Q3DAngle*TMath::Pi()/180.) * sinh(rapidity)*sinh(rapidity));
    thetaCM = asin(p3 * sin(Q3DAngle*TMath::Pi()/180.) / pCM_3);
    cout << "thetaCM: " << thetaCM*180./TMath::Pi() << endl;
    
    //Centre-of-mass 4-momenta
    TLorentzVector fEnergyImpulsionCM_3  = TLorentzVector(pCM_3*sin(thetaCM),0,pCM_3*cos(thetaCM),ECM_3);
    TLorentzVector fEnergyImpulsionCM_4  = fTotalEnergyImpulsionCM - fEnergyImpulsionCM_3;
    
    //Work out the lab 4-momenta
    TLorentzVector fEnergyImpulsionLab_3 = fEnergyImpulsionCM_3;
    fEnergyImpulsionLab_3.Boost(0,0,BetaCM);
    TLorentzVector fEnergyImpulsionLab_4 = fEnergyImpulsionCM_4;
    fEnergyImpulsionLab_4.Boost(0,0,BetaCM);
    
    // Angle in the lab frame
    double ThetaLab3 = fEnergyImpulsionLab_3.Angle(fEnergyImpulsionLab_1.Vect());
    if (ThetaLab3 < 0) ThetaLab3 += M_PI;
    //cout << "ThetaLab3: " << ThetaLab3*180./TMath::Pi() << endl;
    
    //Get the recoil angle - don't actually need this but meh
    double ThetaLab4 = fEnergyImpulsionLab_4.Angle(fEnergyImpulsionLab_1.Vect());
    if (fabs(ThetaLab3) < 1e-6) ThetaLab3 = 0;
    ThetaLab4 = fabs(ThetaLab4);
    if (fabs(ThetaLab4) < 1e-6) ThetaLab4 = 0;
    
    // Kinetic Energies in the lab frame
    double KineticEnergyLab3 = fEnergyImpulsionLab_3.E() - MassEjectile;
    double KineticEnergyLab4 = fEnergyImpulsionLab_4.E() - MassRecoil;
    //cout << "KineticEnergyLab3: " << KineticEnergyLab3 << endl;
    
    // test for total energy conversion
    if (fabs(fTotalEnergyImpulsionLab.E() - (fEnergyImpulsionLab_3.E()+fEnergyImpulsionLab_4.E())) > 1e-6)
        cout << "Problem for energy conservation" << endl;
    
    //return the ejectile KE
    return KineticEnergyLab3;
}

double ComputeExcitationEnergy(double TBeam, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double TEjectile, double Q3DAngle)
{
    
    double pEjectile = sqrt(TEjectile * (TEjectile + 2*MassEjectile));
    
    double pBeam = sqrt(TBeam * (TBeam + 2*MassBeam));
    
    double pRecoil = sqrt(pBeam*pBeam - 2*pBeam*pEjectile*cos(Q3DAngle*TMath::Pi()/180.) + pEjectile*pEjectile);
    //cout << "p4: " << p4 << endl;
    
    double T4 = sqrt(pRecoil*pRecoil + MassRecoil*MassRecoil) - MassRecoil;
    //cout << "T4: " << T4 << endl;
    
    Ex = sqrt( pow(TBeam + MassBeam + MassTarget - TEjectile - MassEjectile,2.) - pRecoil*pRecoil) - MassRecoil;
    
    return Ex; 
}

void TestKinematics(double TBeam, double Q3DAngle, double TargetAngle, double TargetThickness, double CarbonThickness, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Ex)
{
    //The following line will return a Brho value
    //ExToBrho(TBeam, Q3DAngle, TargetAngle, TargetThickness, CarbonThickness, MassBeam, MassTarget, MassEjectile, MassRecoil, Ex);
    double Ex2 =  BrhoToEx(TBeam,
                           Q3DAngle,
                           TargetAngle,
                           TargetThickness,
                           CarbonThickness,
                           MassBeam,
                           MassTarget,
                           MassEjectile,
                           MassRecoil,
                           ExToBrho(TBeam,
                                    Q3DAngle,
                                    TargetAngle,
                                    TargetThickness,
                                    CarbonThickness,
                                    MassBeam,
                                    MassTarget,
                                    MassEjectile,
                                    MassRecoil,
                                    Ex
                           ));
    
    std::cout << "Delta Ex = " << 1000*(Ex - Ex2) << " keV" << std::endl;
}