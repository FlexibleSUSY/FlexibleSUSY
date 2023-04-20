

ParticleDefinitions[GaugeES] = {
      {H0,  {    PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "H^0",
                 OutputName -> "H0" }},                         
      
      
      {Hp,  {    PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 2,
                 LaTeX -> "H^+",
                 OutputName -> "Hp" }}, 

      {phi13, {  PDG -> {0},
                 Width -> 0,
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "\\phi_{1/3}",
                 OutputName -> "phi1"}},

      {phi23, {  PDG -> {0},
                 Width -> 0,
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> {"\\phi_{-2/3}","\\phi_{4/3}"},
                 OutputName -> "phi2"}},
      {phi43, {  PDG -> {0},
                 Width -> 0,
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "\\phi_{4/3}",
                 OutputName -> "phi3"}},

      {VB,   { Description -> "B-Boson"}},                                                   
      {VG,   { Description -> "Gluon"}},          
      {VWB,  { Description -> "W-Bosons"}},          
      {gB,   { Description -> "B-Boson Ghost"}},                                                   
      {gG,   { Description -> "Gluon Ghost" }},          
      {gWB,  { Description -> "W-Boson Ghost"}}
      
      };
      
      
      
      
  ParticleDefinitions[EWSB] = {
            
      
    {hh   ,  {  Description -> "Higgs",
                 PDG -> {25},
                 PDG.IX -> {101000001} }}, 
                 
     {Ah   ,  {  Description -> "Pseudo-Scalar Higgs",
                 PDG -> {0},
                 PDG.IX ->{0},
                 Mass -> {0},
                 Width -> {0} }},                       
      
      
     {Hp,     { Description -> "Charged Higgs", 
                 PDG -> {0},
                 PDG.IX ->{0},
                 Width -> {0}, 
                 Mass -> {0},
                 LaTeX -> {"H^+","H^-"},
                 OutputName -> {"Hp","Hm"},
                 ElectricCharge->1
                 }},  

      {phi13,  { Description -> "Phi 1/3",
                 PDG -> {100},
                 LaTeX -> {"\\phi_{1/3}","\\phi_{1/3}"},
                 OutputName -> "phi1",
                 ElectricCharge->1/3
                 }},

      {phi43,  { Description -> "Charged phi43",
                 PDG -> {101},
                 LaTeX -> {"\\phi_{4/3}","\\phi_{-4/3}"},
                 OutputName -> {"phi2","conjphi2"},
                 ElectricCharge->{4/3,-4/3}
                 }},

      {phi23,  { Description -> "Charged phi23",
                 PDG -> {102},
                 LaTeX -> {"\\phi_{2/3}","\\phi_{-2/3}"},
                 OutputName -> "phi23",
                 ElectricCharge->{2/3,-2/3}
                 }},

      {VP,   { Description -> "Photon"}}, 
      {VZ,   { Description -> "Z-Boson",
      			 Goldstone -> Ah }}, 
      {VG,   { Description -> "Gluon" }},          
      {VWp,  { Description -> "W+ - Boson",
      			Goldstone -> Hp }},         
      {gP,   { Description -> "Photon Ghost"}},                                                   
      {gWp,  { Description -> "Positive W+ - Boson Ghost"}}, 
      {gWpC, { Description -> "Negative W+ - Boson Ghost" }}, 
      {gZ,   { Description -> "Z-Boson Ghost" }},
      {gG,   { Description -> "Gluon Ghost" }},          
                               
                 
      {Fd,   { Description -> "Down-Quarks"}},   
      {Fu,   { Description -> "Up-Quarks"}},   
      {Fe,   { Description -> "Leptons" }},
      {Fv,   { Description -> "Neutrinos" }}                                                              
     
        };    
        
        
        
 WeylFermionAndIndermediate = {
     
    {H,      {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 LaTeX -> "H",
                 OutputName -> "" }},

   {dR,     {LaTeX -> "d_R" }},
   {eR,     {LaTeX -> "e_R" }},
   {lep,     {LaTeX -> "l" }},
   {uR,     {LaTeX -> "u_R" }},
   {q,     {LaTeX -> "q" }},
   {eL,     {LaTeX -> "e_L" }},
   {dL,     {LaTeX -> "d_L" }},
   {uL,     {LaTeX -> "u_L" }},
   {vL,     {LaTeX -> "\\nu_L" }},

   {DR,     {LaTeX -> "D_R" }},
   {ER,     {LaTeX -> "E_R" }},
   {UR,     {LaTeX -> "U_R" }},
   {EL,     {LaTeX -> "E_L" }},
   {DL,     {LaTeX -> "D_L" }},
   {UL,     {LaTeX -> "U_L" }}
        };       


