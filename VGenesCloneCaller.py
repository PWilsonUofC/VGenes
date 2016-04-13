__author__ = 'wilsonp'
import VGenesSeq
from operator import itemgetter
from math import ceil
import os
def CloneCaller(DataList, Duplicates):


    # 'SeqName', 'Sequence', 'GermlineSequence', 'VLocus', 'DLocus', 'JLocus', 'CDR3Length', 'CDR1From', 'CDR1To', 'CDR2From', 'CDR2To', 'CDR3beg', 'CDR3end', 'Mutations', 'IDEvent'

    DataList.sort(key=itemgetter(3,5,6)) #sort data list by Vlocus, Jlocus, and CDR3 length
    i = 0
    StartNew = False
    NewList = []
    PrePools = []
    NewPools = []
    for record in DataList:
        Sequence = record[1]
        CDR3 = Sequence[int(record[11]):int(record[12])]
        CDR1from = int(record[7])
        if len(CDR3) > 0 and CDR1from > 0:
            if i+1 < len(DataList):
                if StartNew == False:
                    NewList.clear()
                    NewList.append(record)
                    # test1 = record[3]
                    # test2 = DataList[i+1][3]
                if record[3] == DataList[i+1][3] and record[5] == DataList[i+1][5] and record[6] == DataList[i+1][6]: #has same V, J and CDR3 length and CDR3 is >2
                    StartNew = True
                    NewList.append(DataList[i+1])
                else:
                    if len(NewList) > 1:

                        PrePools.append(tuple(NewList))
                        NewList.clear()
                    StartNew = False

            i += 1
    if len(NewList) > 1:
        PrePools.append(tuple(NewList))

    SeqDict = []
    i =0
    MutFreqDict = {}

    ClonalPools = []

    for pool in PrePools:
        SeqDict.clear()
        NearPool = []
        for record in pool:
            Seqname = record[0]
            Sequence = record[1]
            NewList.clear()


            CDR3 = Sequence[int(record[11]):int(record[12])]
            # CDR3 = CDR3.upper()
            # CountNuc = CDR3.count('A') + CDR3.count('G') + CDR3.count('C') + CDR3.count('T')
            NewList.append(Seqname)
            NewList.append(CDR3)

            SeqDict.append(tuple(NewList))



            GSequence = record[2]
            CDRSeq = Sequence[int(record[8]):int(record[9])] + Sequence[int(record[10]):int(record[11])]
            GCDRSeq = GSequence[int(record[8]):int(record[9])] + Sequence[int(record[10]):int(record[11])]
            CDRMuts = 0
            for i in range(0, len(CDRSeq)-1):
                try:
                    if CDRSeq[i] != GCDRSeq[i]:
                        CDRMuts +=1
                except:
                    print('miss')
            if CDRMuts > 0:
                CDRMutFreq =  CDRMuts/len(CDRSeq)#+1
            else:
                CDRMutFreq = 0.0
            # the +1 is for seq errors
            Mutations = record[13]
            MutList = Mutations.split(',')
            Vend = int(record[14])

            MutFreqDict[Seqname] = (CDRMutFreq, CDR3, MutList, Vend)

        ClustalOut = VGenesSeq.ClustalO(SeqDict, 1000, False)

        outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes','ClustalOmega', 'my-out-seqs.fa')
        Aligned = VGenesSeq.readClustalOutput(outfilename)
        Similars = {}

        SimPool = []
        Matches = []
        Different = []

        i = 0
        LastSeqName = ''

        for seq in range(0,len(Aligned)-1):
            SeqName = Aligned[seq][0]
            if SeqName == '':
                print('bad')
            CDR3MutFreq = MutFreqDict[SeqName][0]
            CDR3 = MutFreqDict[SeqName][1]
            CDR3dif = CDR3MutFreq * len(CDR3)
            SimPool.clear()
            SimPool.append(SeqName)
            # DupPool.clear()
            # DupPool.append(SeqName)
            for comp in range(seq+1, len(Aligned)):
                CSeqName = Aligned[comp][0] #comparison SeqName, etc
                CCDR3MutFreq = MutFreqDict[CSeqName][0]
                CCDR3 = MutFreqDict[CSeqName][1]
                if Duplicates == False:
                    CCDR3dif = CCDR3MutFreq * len(CCDR3)
                    FiftyPercent = int(0.5*len(CDR3))

                    differences = 0
                    AllowedDif = CDR3dif + CCDR3dif
                    if AllowedDif < 1 and CDR3dif >0 and CCDR3dif >0: AllowedDif = 1
                    AllowedDif+=2
                    AllowedDif = ceil(AllowedDif)

                    for j in range(0, len(CDR3)-1):
                        try:
                            if CDR3[j] != CCDR3[j]:
                                differences += 1
                        except:
                            print('stop')
                else:
                    if CDR3 == CCDR3:

                        differences = 1
                    else:
                        differences = 0


                if Duplicates == False:
                    if differences <= AllowedDif:
                        # SimPool.append(SeqName)
                        SimPool.append(CSeqName)
                    elif differences <= FiftyPercent:
                        MutList = MutFreqDict[SeqName][2]
                        # NumMuts = len(MutList)
                        # Vend = MutFreqDict[SeqName][3]
                        MutListComp = MutFreqDict[CSeqName][2]
                        # NumMutsComp = len(MutListComp)
                        # VendComp = MutFreqDict[CSeqName][3]
                        Matches.clear()
                        Matches = [element for element in MutList if element in MutListComp]
                        if len(Matches)>3:  #need get average shared + 2SD from clones for this instead of arbitrary 3
                            SimPool.append(CSeqName)
                else:
                    if differences == 0:
                        MutList = MutFreqDict[SeqName][2]
                        MutListComp = MutFreqDict[CSeqName][2]
                        if MutList == MutListComp:
                            SimPool.append(CSeqName)





                # i+=1
            if len(SimPool) > 1:

                # because aligned now can check if this one has everything in previous
                # already and only add to Similars dict if something differs
                # only need to find when stuff is different from list that already exists
                if len(Similars) >0:
                    SetLast = list(Similars[LastSeqName]) #set(Similars[LastSeqName])
                    Matches.clear()
                    # Different.clear()
                    Matches = [element for element in SimPool if element in SetLast]
                    if len(Matches) > 0: #if some match
                        SetLast += [element for element in SimPool if element not in SetLast]
                        Entry = tuple(SetLast)
                        Similars[LastSeqName] = Entry  #update

                    else:

                        Entry = tuple(SimPool)  #Start  new one
                        Similars[SeqName] = Entry
                        LastSeqName = SeqName

                elif len(Similars) == 0:
                    Entry = tuple(SimPool)
                    Similars[SeqName] = Entry
                    LastSeqName = SeqName


        Different.clear()

        Pools = list(Similars.keys())
        # i=0
        KeepList = []
        if len(Similars) >1:

            for i in range(0,len(Pools)-1):
                SimPoolName = Pools[i]
                SimPool = list(Similars[SimPoolName])
                if len(SimPool) != 0:
                    for j in range(i+1,len(Pools)):
                        NextName = Pools[j]
                        nextset = list(Similars[NextName])
                        Matches = [element for element in SimPool if element in nextset]
                        if len(Matches) > 0:
                            Different = [element for element in nextset if element not in SimPool]

                            SimPool = SimPool + Different
                            Similars[SimPoolName] = tuple(SimPool)
                            Similars[NextName] = ''

            # could marge  test name and make rest merged empty so no further matches
            # in end iterate through and make list of tuples of all final merged lists
            # Different.append(list(set))

        NewPools.clear()
        for i in range(0,len(Pools)):  #code to remove empty, merged pools and covert from dictionary to lists
            SimPoolName = Pools[i]
            SimPool = list(Similars[SimPoolName])
            if len(SimPool) != 0:
                NewPools.append(SimPool)


        for pool in NewPools:  #removes duplicates from list, faster then finding records and remarking as clonal later
            NearPool = pool
            NearPool.sort()
            Lennear = len(NearPool)-1
            for k in range (0, Lennear):
                try:
                    if k < (len(NearPool)-1):
                        if k > 0 and NearPool[k] == NearPool[k-1]:
                            NearPool.remove(NearPool[k])
                except:
                    print('Stop')

            if len(NearPool) >1: ClonalPools.append(tuple(NearPool))  #after dups removed if still pool add to final list of clonal pools



    return ClonalPools

# calculates CDR mutation frequency and uses that as the cutoff for
# mutation frequency allowed in the CDR3
#todo need to calculat standard deviation as well and allow within
# it then orders sequences by clustalO and scores CDR3 differences
# in order to see if different from CDR frequnecies

