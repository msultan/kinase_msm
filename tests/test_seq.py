#!/bin/env/python
'''
Tests to make sure the sequences were copied right
'''

from kinase_msm.kinases import *

btk_test = "GSWEIDPKDLTFLKELGTGQFGVVKYGKWRGQYDVAIKMIKEGSMSEDEFIEEAKVMMNLSHEKLVQLY\
GVCTKQRPIFIITEYMANGCLLNYLREMRHRFQTQQLLEMCKDVCEAMEYLESKQFLHRDLAARNCLVNDQGVVKVSDFG\
LSRYVLDDEYTSSVGSKFPVRWSPPEVLMYSKFSSKSDIWAFGVLMWEIYSLGKMPYERFTNSETAEHIAQGLRLYRPHLA\
SEKVYTIMYSCWHEKADERPTFKILLSNILDVMDEES"

her2_0p3_test = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikvlrentspkankei\
ldeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiakgmsfledvrlvhrdla\
arnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvwsygvtvwelmtfgakpydgip\
areipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqrfvviqned".upper()

itk_test = "GKWVIDPSELTFVQEIGSGQFGLVHLGYWLNKDKVAIKTIREGAMSEEDFIEEAEVMMKLSHPKLVQ\
LYGVCLEQAPICLVFEFMEHGCLSDYLRTQRGLFAAETLLGMCLDVCEGMAYLEEACVIHRDLAARNCLVGENQVIKVS\
DFGMTRFVLDDQYTSSTGTKFPVKWASPEVFSFSRYSSKSDVWSFGVLMWEVFSEGKIPYENRSNSEVVEDISTGFRLY\
KPRLASTHVYQIMNHCWKERPEDRPAFSRLLRQLAEIAESGL"

yes_test = "KDAWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTKVAIKTLKPGTMMPEAFLQEAQIMKKLRHDKLV\
PLYAVVSEEPIYIVTEFMSKGSLLDFLKEGDGKYLKLPQLVDMAAQIADGMAYIERMNYIHRDLRAANILVGENLVCKIA\
DFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILQTELVTKGRVPYPGMVNREVLEQVERGYRMPC\
PQGCPESLHELMNLCWKKDPDERPTFEYIQSFLEDY"


def test_btk_seq():
    assert kinase_seq["btk_asp"] == kinase_seq["btk_ash"] == btk_test
    return


def test_itk_seq():
    assert kinase_seq["btk_asp"] == kinase_seq["btk_ash"] == btk_test
    return


def test_her2_seq():
    for kinase in kinase_set["her2"]:
        if kinase == "her2_0p3x":
            assert kinase_seq[kinase] == her2_0p3_test
        else:
            assert kinase_seq[kinase] != her2_0p3_test
    return


def test_sfk_seq():
    for kinase in kinase_set["sfk"]:
        if kinase == "sfk_yes":
            assert kinase_seq[kinase] == yes_test
        else:
            assert kinase_seq[kinase] != yes_test
            assert kinase_seq[kinase] != her2_0p3_test
    return
