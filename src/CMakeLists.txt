cmake_minimum_required(VERSION 3.15)
project(Arbores)

set(CMAKE_C_STANDARD 99)

include_directories(.)
include_directories(../include)

add_executable(Arbores
        main.c
        ARGdecomposer.c
        backtracking.c
        commandlineoutput.c
        constants.c
        data.c
        debugging.c
        exhaustiveSearch.c
        fileaccess.c
        freeTimes.c
        graph2tikz.c
        initialisation.c
        initialTimeProposal.c
        jittering.c
        likelihood.c
        MCMCutils.c
        mt19937ar.c
        pathutils.c
        randomness.c
        recombinationTimes.c
        results.c
        shrub.c
        smcPrior.c
        sorting.c
        timeadjustment.c
        treeutils.c
        utils.c)
