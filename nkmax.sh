#!/bin/bash
#SBATCH -p cosc6326
#SBATCH -t 01:00:00
#SBATCH -N 8
#SBATCH --ntasks-per-node 4

export TS="$(date +%Y-%m-%d.%H.%M.%S)"
export DOUT="./results/nkmax/${TS}"
mkdir -p ${DOUT}

for ((k = 1; k <= 32; k *= 2)); do
    for ((n = $[2**10]; n <= $[2**28]; n *= 2)); do
	mpirun -np $k ./nkmax -- -n $n 1> ${DOUT}/nkmax.$k.$n.log 2> /dev/null
    done
done
