if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    #---
    lnums = [ 36, 108  ]
    string=open('simulations.py').readlines() #--- python script
#    lnums = [ 39, 116, 121 ]
#    string=open('simulations-ncbj-slurm.py').readlines() #--- python script
    #---
    Temps=			{
                0:400,
                1:500,
                2:600,
                3:700,
                4:800,
                5:900,
                6:1000,
            }
    #---
    count = 0
    for keys_t in Temps:
                temp = Temps[keys_t]
            #---	densities
                inums = lnums[ 0 ] - 1
                string[ inums ] = "\t5:\'pureHydrogenDiffusion_md_MultipleTemps/temp%s\',\n"%(keys_t) #--- change job name
#                string[ inums ] = "\t7:\'biCrystalMultipleTemp2nd/temp%s\',\n"%(keys_t) #--- change job name
            #---
                inums = lnums[ 1 ] - 1
#                string[ inums ] = "\t\'p3\':\' data_minimized.txt init_xyz.conf %%s %s\'%%(os.getcwd()+\'/lmpScripts\'),\n"%temp
    #				string[ inums ] = "\t\'p21\':\' %%s 3.52 %s 18.0 26.0 data_init.txt 2 2 1.0 0.0\'%%(os.getcwd()+\'/lmpScripts\'),\n"%temp
    #				string[ inums ] = "\t\'p7\':\' sortieproc.0 %s Topo_ignore\',\n"%(temp)

                string[ inums ] = "\t72:\' -var seed %%s -var buff 0.0 -var buffy 0.0 -var Tinit %s -var T %s -var nevery 100 -var ParseData 1 -var DataFile data_minimized.txt -var DumpFile dumpThermalized.xyz -var WriteData equilibrated.dat\'%%np.random.randint(1001,9999),\n"%(temp,temp)

#                 inums = lnums[ 2 ] - 1
#                 string[ inums ] = "\t12:\' -var buff 0.0 -var buffy 5.0 -var T %s -var swap_every 100 -var swap_atoms 267 -var rn %%s -var dump_every 100 -var ParseData 1 -var DataFile equilibrated.dat -var DumpFile traj.dump\'%%np.random.randint(1001,100000),\n"%(temp)
            #---
                sfile=open('junk%s.py'%count,'w');sfile.writelines(string);sfile.close()
                os.system( 'python2 junk%s.py'%count )
                os.system( 'rm junk%s.py'%count )
                count += 1
