#!/usr/bin/env bash
java -Xmx5120m -Xms2048m -classpath /opt/gurobi652/linux64/lib/gurobi.jar:.: $1 $2 $3 $4 $5 $6
