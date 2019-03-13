from subprocess import call

fstem = "../Results/Data/sim/instance"
for i in range(100):
    if i == 51:
        continue
    cmd = "java -jar cal2.jar"
    #infile = "../Results/Data/sim/unclust/instance%d_readcounts.txt" % i
    infile = "../Results/Data/sim/instance%d_cl_readcounts.txt" % i

    call(cmd.split() + [infile])
