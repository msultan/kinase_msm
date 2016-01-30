#!/bin/env python
import numpy as np
import mdtraj as mdt
from .tica_utils import _map_tic_component
import argparse

from msmbuilder.utils import verboseload
import pandas as pd

VMDSCRIPT = """

    source /Users/muneeb/Documents/research/sscache.tcl
    set mol [mol new {traj_fn} step {step} waitfor all ]
    mol addfile {top_fn} waitfor all
    set nf [molinfo $mol get numframes]
    animate delete beg [expr $nf-1] end $nf
    set mid [molinfo top]


    animate goto 0
    start_sscache
    for {{set i 0}} {{$i < $nf}} {{incr i}} {{
    puts $i
    animate goto $i
    sleep 0.15
    }}


    set sel [atomselect $mol all]
    set nf [molinfo $mol get numframes]
    set fp [open {dat_fn} "r"]
    set line ""
    gets $fp line
    for {{set i 0}} {{$i < $nf}} {{incr i}} {{
        puts $i
        $sel frame $i
        $sel set user $line
    }}

    close $fp
    $sel delete

    mol delrep 0 top
    mol color User User
    mol representation NewCartoon 0.3 10.0 4.1 0
    mol selection {{protein}}
    mol addrep top
    mol smoothrep top 0 5

    mol representation Licorice 0.3 12.0 12.0
    mol color User
    mol selection {{user > {cutoff} and noh}}
    mol addrep top
    mol selupdate 1 top 1
    mol colupdate 1 top 1
    mol scaleminmax top 1 auto

    """


def write_vmd_file(tcl_fn, traj_fn, top_fn, stride, dat_fn, cutoff):
    with open(tcl_fn, 'w') as f:
        f.write(str(VMDSCRIPT.format(traj_fn=traj_fn,
                                     top_fn=top_fn,
                                     step=stride,
                                     dat_fn=dat_fn,
                                     cutoff=cutoff)))

    return


def tica_to_vmd(df, tica_mdl, tic_index, traj_fn,
                    top_fn, trj=None, stride=1,
                    dat_fn=None, out_file=None,
                    cutoff=0.75):

    if trj is None:
        trj = mdt.load(traj_fn, top=top_fn)

    tic_component = tica_mdl.components_[tic_index,:]
    assert(len(df) == len(tic_component))

    atom_imp, res_imp = _map_tic_component(tic_component, df, trj)

    if dat_fn is None:
        dat_fn = "important%d.txt"%tic_index

    np.savetxt(dat_fn, atom_imp, fmt="%.5f")

    if out_file is None:
        out_file = "ent%d.tcl"%tic_index

    cutoff = cutoff * np.max(res_imp)

    write_vmd_file(tcl_fn=out_file,
                   traj_fn= traj_fn,
                   top_fn = top_fn,
                   stride=stride,
                   dat_fn=dat_fn,
                   cutoff=cutoff)

    return


def parse_commandline():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--trj', dest='t',
                        default='./traj.xtc',
          help='Trajectory file')
    parser.add_argument('-p', '--top', dest='p',
                        default="./prot.pdb",
          help='Reference Top file')
    parser.add_argument('-c', '--tica_file', dest='c',
                        default="./tica.h5",
          help='Tica object')
    parser.add_argument('-i', '--index', dest='i',
                        type=int,default=0,
          help='Tic to plot')
    parser.add_argument('-u', '--cutoff', dest='u',
                        type=float, default=0.75,
          help='Cutoff for tic')
    parser.add_argument('-s', '--stride', dest='s',
                        type=int,default="1",
          help='Stride to use while loading and writing')
    parser.add_argument('-o', '--out_file', dest='o',
                        type=str,default="ent.tcl",
          help='Name of output tcl file')
    parser.add_argument('-d', '--describer', dest='d',
                        type=str,help='List of dictionaries')

    args = parser.parse_args()
    return args

def main():
    args = parse_commandline()
    traj_file = args.t
    top_file  = args.p
    tica_file = args.c
    tic_index = args.i
    out_file = args.o
    stride = args.s
    describer = args.d
    cutoff = args.u
    #load stuff
    trj = mdt.load(traj_file,top=top_file,stride=stride)
    tica_mdl = verboseload(tica_file)
    df = pd.DataFrame(verboseload(describer))


    dat_fn = "importances_{}.txt".format(out_file)
    tcl_fn = "{}.tcl".format(out_file)

    tica_to_vmd(df, tica_mdl, tic_index, traj_file, top_file,
                trj, stride, dat_fn, tcl_fn, cutoff)



if __name__ == "__main__":
    main()
