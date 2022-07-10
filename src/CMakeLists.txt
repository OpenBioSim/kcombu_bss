

cmake_minimum_required(VERSION 3.00)

project(KCOMBU_BSS C)

set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_RELEASE})

add_definitions("-D_USE_MATH_DEFINES")

set(COMMON_SOURCES
    match.c
    io_match.c
    buildup.c
    clique.c
    c_clique.c
    MrgSrtMATCH.c
    ioLib.c
    ioLINE.c
    2DMAP.c
    ioSDF.c
    ioKCF.c
    ioMOL2.c
    ioSMILES.c
    kcombu_io_CIF.c
    MOLECULE_from_CCD_CIF.c
    gen2D.c
    disgeo2D.c
    qRMS.c
    2dRMS.c
    PCAfit.c
    molprop.c
    molring.c
    moltopodis.c
    molpermu.c
    options.c
    ringblock.c
    RingDesc.c
    vector3.c
    outPS.c
    energies.c
    stereo_check.c
    ReadOpt.c
    ioPDB.c
)

set(CKCOMBU_SOURCES
    ckcombu.c
    ${COMMON_SOURCES}
)

set(DKCOMBU_SOURCES
    dkcombu.c
    AtmPairDesc.c
    DualAtmPair.c
    MrgSrtLIBMOL.c
    stereo_check.c
    ${COMMON_SOURCES}
)

set(FKCOMBU_SOURCES
    fkcombu.c
    transform.c
    energies.c
    steep_descen.c
    tpdisrest.c
    stamp_transf.c
    stereo_adjust.c
    ${COMMON_SOURCES}
)

set(LKCOMBU_SOURCES
    lkcombu.c
    MrgSrtLIBMOL.c
    ${COMMON_SOURCES}
)

set(PKCOMBU_SOURCES
    pkcombu.c
    ${COMMON_SOURCES}
)

set(SMICOMP_SOURCES
    SmiComp.c
)

set(TREEUN_SOURCES
    TreeUN.c
)

set(PROPSEARCH_SOURCES
    propsearch.c
)

add_executable(ckcombu_bss ${CKCOMBU_SOURCES})
add_executable(dkcombu_bss ${DKCOMBU_SOURCES})
add_executable(fkcombu_bss ${FKCOMBU_SOURCES})
add_executable(lkcombu_bss ${LKCOMBU_SOURCES})
add_executable(pkcombu_bss ${PKCOMBU_SOURCES})
add_executable(smicomp_bss ${SMICOMP_SOURCES})
add_executable(treeun_bss ${TREEUN_SOURCES})
add_executable(propsearch_bss ${PROPSEARCH_SOURCES})

# installation
install(TARGETS ckcombu_bss dkcombu_bss fkcombu_bss
                lkcombu_bss pkcombu_bss)
