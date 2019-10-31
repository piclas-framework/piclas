#!/bin/bash

MatchesToIgnore=(
    "PartState(dir,"
    "PartState(iDir"
    "PartState(iDir,"
    "PartState(iPartState,"
    "PartState(React1Inx"
    "PartState(ind,"
    "PartState(dir("
    "PartState(dir+3,"
    "PartState(iPartState"
)

for i in "${MatchesToIgnore[@]}"; do
  # Concacenate items (use first item if variable excludes is empty or non-existent)
    excludes=${excludes:+$excludes"\|"}$i
done

#echo $excludes
#echo ""
clear

grep -nri --include=*.f90 --color "FieldAtParticle([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "Pt([[:alpha:]]" ../src | grep -inv "PartLambdaAccept(\|BGK_octree_adapt(\|BGK_quadtree_adapt(\|PartInsSubSidesAdapt(" | grep -in --color "Pt("

grep -nri --include=*.f90 --color "lastpartpos([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "partstate([[:alpha:]]" ../src | grep -inv "${excludes}" | grep -in --color PartState

grep -nri --include=*.f90 --color "DSMC_RHS([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "ExtPartState([[:alpha:]]" ../src | grep -inv "ExtPartState(ind," | grep -in --color ExtPartState

grep -nri --include=*.f90 --color "PartStateIntEn([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "PartStateN([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "Pt_temp([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "PartStage([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "PEM%NormVec([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "PartPosRef([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "PartData([[:alpha:]]" ../src | grep -inv "ConvertPartData(InputStateFile\|SurfPartData(offsetnSurfPart\|SurfPartData(iPart\|PartData(PartDataSize\|SurfPartData(iOffset\|PartData(INT(PartDataSize,IK)," | grep -in --color PartData

# octree stuff
grep -nri --include=*.f90 --color "MappedPartStates([[:alpha:]]" ../src

grep -nri --include=*.f90 --color "MappedPart_ChildNode([[:alpha:]]" ../src
















