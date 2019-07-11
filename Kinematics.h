#include <TLorentzVector.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>

//p0, p1 and p2 are the parameters that come from fitting E vs range data for the target/ejectile
double DoEnergyLoss(double TInitial, double Thickness, double a, double b, double c)
{
    double Range = a*pow(TInitial,2.) + b*TInitial + c;
    Range -= Thickness;
    return (-b + sqrt(b*b - 4*a*(c-Range)))/2./a;
}

double ExToBrho(double TBeam, double Q3DAngle, double TargetAngle, double TargetThickness, double CarbonThickness, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Ex);

double ConvertKEToBrho(double TEjectile, double MassEjectile)
{
    return 1.e6/TMath::C() * sqrt(TEjectile * (TEjectile + 2*MassEjectile));
}

double ConvertBrhoToKE(double Brho, double MassEjectile)
{
    return sqrt(pow(TMath::C()/1e6 * Brho,2.) + pow(MassEjectile,2.)) - MassEjectile;
}

double BrhoToEx(double TBeam, double Q3DAngle, double TargetAngle, double TargetThickness, double CarbonThickness, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Brho);

double ComputeEjectileEnergy(double TBeam, double Q3DAngle, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double Ex);

double ComputeExcitationEnergy(double TBeam, double MassBeam, double MassTarget, double MassEjectile, double MassRecoil, double TEjectile, double Q3DAngle);

