#ifndef HELPER_H
#define HELPER_H

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <algorithm>

namespace SelfVetoPolicy
{
   enum SelfVetoPolicy
     {
	selfVetoNone=0, selfVetoAll=1, selfVetoFirst=2
     };
}

float GetDeltaR(float,float,float,float);

#endif
