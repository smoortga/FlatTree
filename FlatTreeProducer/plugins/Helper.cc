#include "FlatTree/FlatTreeProducer/interface/Helper.hh"
#include <iostream>
#include "TMath.h"
#include "TLorentzVector.h"

float GetDeltaR(float eta1,float phi1,float eta2,float phi2)
{
   float DeltaPhi = TMath::Abs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
