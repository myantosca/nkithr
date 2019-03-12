#!/bin/bash
#SBATCH -p cosc6326
#SBATCH -t 01:00:00
#SBATCH -N 8
#SBATCH --ntasks-per-node 4

export TS="$(date +%Y-%m-%d.%H.%M.%S)"
export DOUT="./results/nkith/${TS}"
mkdir -p ${DOUT}

for ((k = 1; k <= 32; k *= 2)); do
    for ((n = $[2**10]; n <= $[2**28]; n *=2)); do
	i=$n
	mpirun -np $k ./nkith -- -n $n -i $i 1> ${DOUT}/nkith.max.$k.$n.log 2> /dev/null
	i=$[$n/2]
	mpirun -np $k ./nkith -- -n $n -i $i 1> ${DOUT}/nkith.med.$k.$n.log 2> /dev/null
	i=1
	mpirun -np $k ./nkith -- -n $n -i $i 1> ${DOUT}/nkith.min.$k.$n.log 2> /dev/null
    done
done
