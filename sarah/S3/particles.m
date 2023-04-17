

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
     
      {psi1, {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "\\psi_{1/3}",
                 OutputName -> "psi1"}},
      
      {psi2, {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> {"\\psi_{-2/3}","\\psi_{4/3}"},
                 OutputName -> "psi2"}},
      {psi3, {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "\\psi_{4/3}",
                 OutputName -> "psi3"}},
    
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
                 
                 
      {psi1,  { Description -> "Charged psi1", 
                 PDG -> {100},
                 PDG.IX ->{0},
                 LaTeX -> {"\\psi_{1/3}","\\psi_{1/3}"},
                 OutputName -> "psi1",
                 ElectricCharge->1/3
                 }},   
      
      {psi2,  { Description -> "Charged psi2", 
                 PDG -> {101},
                 PDG.IX ->{0},
                 LaTeX -> {"\\psi_{4/3}","\\psi_{-4/3}"},
                 OutputName -> {"psi2","conjpsi2"},
                 ElectricCharge->{4/3,-4/3}
                 }},      
    
      {psi3,  { Description -> "Charged psi3", 
                 PDG -> {102},
                 PDG.IX ->{0},
                 LaTeX -> {"\\psi_{-2/3}","\\psi_{2/3}"},
                 OutputName -> "psi3",
                 ElectricCharge->{-2/3,2/3}
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


