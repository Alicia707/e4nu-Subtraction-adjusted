#!/bin/bash

GenieChoice=1
ROTATIONS=100
percentEvents=100

./genie_analysis C12 1161 $GenieChoice $ROTATIONS $percentEvents
./genie_analysis 4He 2261 $GenieChoice $ROTATIONS $percentEvents
./genie_analysis C12 2261 $GenieChoice $ROTATIONS $percentEvents
./genie_analysis 56Fe 2261 $GenieChoice $ROTATIONS $percentEvents
./genie_analysis 4He 4461 $GenieChoice $ROTATIONS $percentEvents
./genie_analysis C12 4461 $GenieChoice $ROTATIONS $percentEvents
./genie_analysis 56Fe 4461 $GenieChoice $ROTATIONS $percentEvents
