

Int_t getTargetBin() {
    HVertex PrimVertexReco = ((HEventHeader*)(gHades->getCurrentEvent()->getHeader()))->getVertexReco();
    Float_t vz=PrimVertexReco.getZ();

    Int_t target = -1;
    if (vz >= -63.0 && vz < -60.0) target = 0;
    if (vz >= -60.0 && vz < -56.5) target = 1;
    if (vz >= -56.5 && vz < -52.5) target = 2;
    if (vz >= -52.5 && vz < -49.5) target = 3;
    if (vz >= -49.5 && vz < -45.5) target = 4;
    if (vz >= -45.0 && vz < -42.0) target = 5;
    if (vz >= -42.0 && vz < -39.0) target = 6;
    if (vz >= -39.0 && vz < -35.5) target = 7;
    if (vz >= -35.5 && vz < -31.5) target = 8;
    if (vz >= -31.5 && vz < -28.0) target = 9;
    if (vz >= -28.0 && vz < -24.0) target = 10;
    if (vz >= -24.0 && vz < -21.5) target = 11;
    if (vz >= -21.5 && vz < -17.0) target = 12;
    if (vz >= -17.0 && vz < -13.0) target = 13;
    if (vz >= -13.0) target = 14;
    

    return target;
}


