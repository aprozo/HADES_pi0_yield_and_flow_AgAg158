Int_t GetIonPdg(Int_t z, Int_t a, Int_t i, Int_t l=0, Int_t sign=+1) 
{
// taken from  [geant3.git] / TGeant3 / TGeant3.cxx
// modified to convention taken from PDG: 
// http://pdg.lbl.gov/2015/reviews/rpp2015-rev-monte-carlo-numbering.pdf

/*  
Nuclear codes are given as 10-digit numbers ±10LZZZAAAI.
For a (hyper)nucleus consisting of np protons, nn neutrons and
nΛ Λ’s, A = np + nn + nΛ gives the total baryon number, Z = np
the total charge and L = nΛ the total number of strange quarks.
I gives the isomer level, with I = 0 corresponding to the ground
state and I > 0 to excitations, see [4], where states denoted
m, n, p, q translate to I = 1 − 4. As examples, the deuteron
is 1000010020 and 235U is 1000922350. To avoid ambiguities,
nuclear codes should not be applied to a single hadron, like p, n
or Λ0, where quark-contents-based codes already exist.
*/
//       ±10LZZZAAAI
  return (1000000000 
       +    10000000*l
       +       10000*z
       +          10*a
       +           1*i)*sign;
}  

void loadGeantIons(){

//--------------------------------------------------------------------------------------
    cout << endl;
    cout<<"-------------------------------------------"<<endl;
    cout <<"Adding new Particles(Fragments) - remove Patricle 61-69 -- Re-Intialize HPhysicsConstants with GEANT3 numbering ... \n" <<endl;
    for(UInt_t jj = 61 ; jj < 70; jj ++) HPhysicsConstants::removeParticle(jj);
    
    // Mass from GEANT3 gpions.F
    // and information from shield to evt convert.sh

//    HPhysicsConstants::addParticle(Int_t id, TString name, Double_t mass,Int_t chrg,Int_t leptchrg,Int_t barychrg,Int_t strange)
    HPhysicsConstants::addParticle(61,"Li 3 6",      5603.05,      3,  0,  6, 0);  //
    HPhysicsConstants::addParticle(62,"Li 3 7",      6535.36,      3,  0,  7, 0);  //
    HPhysicsConstants::addParticle(63,"Be 4 7",      6536.22,      4,  0,  7, 0);  //
    HPhysicsConstants::addParticle(64,"Be 4 9",      8394.79,      4,  0,  9, 0);  //
    HPhysicsConstants::addParticle(65,"Be 5 10",     9326.99,      5,  0, 10, 0);  //
    HPhysicsConstants::addParticle(66,"Be 5 11",    10255.10,      5,  0, 11, 0);  //
    HPhysicsConstants::addParticle(67,"C  6 12",    11177.93,      6,  0, 12, 0);  //
    HPhysicsConstants::addParticle(68,"N  7 14",    13043.78,      7,  0, 14, 0);  //
    HPhysicsConstants::addParticle(69,"O  8 16",    14899.17,      8,  0, 16, 0);  //
    HPhysicsConstants::addParticle(70,"F  9 19",    17696.90,      9,  0, 19, 0);  //70: F  9 19
    HPhysicsConstants::addParticle(71,"Ne 10 20",   18622.84,     10,  0, 20, 0);  //71: Ne 10 20
    HPhysicsConstants::addParticle(72,"Na 11 23",   21414.83,     11,  0, 23, 0);  //Na 11 23
    HPhysicsConstants::addParticle(73,"Mg 12 24",   22341.93,     12,  0, 24, 0);  //Mg 12 24
    HPhysicsConstants::addParticle(74,"Al 13 27",   25133.14,     13,  0, 27, 0);  //Al 13 27
    HPhysicsConstants::addParticle(75,"Si 14 28",   26060.34,     14,  0, 28, 0);  //Si 14 28
    HPhysicsConstants::addParticle(76,"P  15 31",   28851.88,     15,  0, 31, 0);  //P  15 31
    HPhysicsConstants::addParticle(77,"S  16 32",   29781.80,     16,  0, 32, 0);  //S  16 32
    HPhysicsConstants::addParticle(78,"Cl 17 35",   32573.28,     17,  0, 35, 0);  //Cl 17 35
    HPhysicsConstants::addParticle(79,"Ar 18 36",   33503.56,     18,  0, 36, 0);  //Ar 18 36
    HPhysicsConstants::addParticle(80,"K  19 39",   36294.47,     19,  0, 39, 0);  //K  19 39
    HPhysicsConstants::addParticle(81,"Ca 20 40",   37224.92,     20,  0, 40, 0);  //Ca 20 40
    HPhysicsConstants::addParticle(82,"Sc 21 45",   41876.17,     21,  0, 45, 0);  //Sc 21 45
    HPhysicsConstants::addParticle(83,"Ti 22 48",   44663.24,     22,  0, 48, 0);  //Ti 22 48
    HPhysicsConstants::addParticle(84,"V  23 51",   47454.01,     23,  0, 51, 0);  //V  23 51
    HPhysicsConstants::addParticle(85,"Cr 24 52",   48382.28,     24,  0, 52, 0);  //Cr 24 52
    HPhysicsConstants::addParticle(86,"Mn 25 55",   51174.47,     25,  0, 55, 0);  //Mn 25 55
    HPhysicsConstants::addParticle(87,"Fe 26 56",   52103.07,     26,  0, 56, 0);  //Fe 26 56
    HPhysicsConstants::addParticle(88,"Co 27 59",   54895.93,     27,  0, 59, 0);  //Co 27 59
    HPhysicsConstants::addParticle(89,"Ni 28 58",   53966.44,     28,  0, 58, 0);  //Ni 28 58
    HPhysicsConstants::addParticle(90,"Cu 29 63",   58618.56,     29,  0, 63, 0);  //Cu 29 63
    HPhysicsConstants::addParticle(91,"Zn 30 64",   59549.63,     30,  0, 64, 0);  //Zn 30 64
    HPhysicsConstants::addParticle(92,"Ge 32 74 ",  68857.15,     32,  0, 74, 0);  //
    HPhysicsConstants::addParticle(93,"Se 34 80 ",  74441.78,     34,  0, 80, 0);  //
    HPhysicsConstants::addParticle(94,"Kr 36 84 ",  78163.09,     36,  0, 84, 0);  //
    HPhysicsConstants::addParticle(95,"Sr 38 88 ",  81883.58,     38,  0, 88, 0);  //
    HPhysicsConstants::addParticle(96,"Zr 40 90 ",  83745.71,     40,  0, 90, 0);  //
    HPhysicsConstants::addParticle(97,"Mo 42 98 ",  91198.32,     42,  0, 98, 0);  //
    HPhysicsConstants::addParticle(98,"Pd 46 106",  98649.97,     46,  0, 106, 0);  //
    HPhysicsConstants::addParticle(99,"Cd 48 114", 106109.97,     48,  0, 114, 0);  //
    HPhysicsConstants::addParticle(100,"Sn 50 120",111688.21,     50,  0, 120, 0);  //
    HPhysicsConstants::addParticle(101,"Xe 54 132",122867.96,     54,  0, 132, 0);  //
    HPhysicsConstants::addParticle(102,"Ba 56 138",128457.93,     56,  0, 138, 0);  //
    HPhysicsConstants::addParticle(103,"Ce 58 140",130321.11,     58,  0, 140, 0);  //
    HPhysicsConstants::addParticle(104,"Sm 62 152",141512.36,     62,  0, 152, 0);  //
    HPhysicsConstants::addParticle(105,"Dy 66 164",152699.09,     66,  0, 164, 0);  //
    HPhysicsConstants::addParticle(106,"Yb 70 174",162022.45,     70,  0, 174, 0);  //
    HPhysicsConstants::addParticle(107,"W  74 184",171349.24,     74,  0, 184, 0);  //
    HPhysicsConstants::addParticle(108,"Pt 78 194",180675.13,     78,  0, 194, 0);  //
    HPhysicsConstants::addParticle(109,"Au 79 197",183473.24,     79,  0, 197, 0);  //
    HPhysicsConstants::addParticle(110,"Hg 80 202",188134.51,     80,  0, 110, 0);  //
    HPhysicsConstants::addParticle(111,"Pb 82 208",193729.07,     82,  0, 111, 0);  //
    HPhysicsConstants::addParticle(112,"U  92 238",221742.95,     92,  0, 112, 0);  //

    HPhysicsConstants::print();
    cout<<"-------------------------------------------"<<endl;
    cout<< endl;


    return;
} 
