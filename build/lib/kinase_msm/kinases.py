#!/bin/env/python

series_set = ["sfk", "mek", "her2", "btk", "itk"]

# dictionary of all kinase sequences that I have run and their corresponding \
# project ids
kinase_set = {}
kinase_seq = {}
proj_list = {}

ligand_set = ["apo", "ATP", "LIG"]

kinase_set["her2"] = ["her2_wt", "her2_40x", "her2_22x", "her2_5x",
                      "her2_4p7x", "her2_4p7x", "her2_3p3x", "her2_2p2x",
                      "her2_1p2x", "her2_1x", "her2_lap", "her2_0p3x"]
kinase_set["sfk"] = ["sfk_fyn", "sfk_fgr", "sfk_hck", "sfk_lck",
                     "sfk_lyn", "sfk_yes", "sfk_blk"]
kinase_set["itk"] = ["itk"]
kinase_set["btk"] = ["btk_asp", "btk_ash"]
kinase_set["mek"] = ["mek2"]


kinase_seq["btk_asp"] = "gsweidpkdltflkelgtgqfgvvkygkwrgqydvaikmikegsmsede\
fieeakvmmnlsheklvqlygvctkqrpifiiteymangcllnylremrhrfqtqqllemckdvceameyleskq\
flhrdlaarnclvndqgvvkvsdfglsryvlddeytssvgskfpvrwsppevlmyskfssksdiwafgvlmweiys\
lgkmpyerftnsetaehiaqglrlyrphlasekvytimyscwhekaderptfkillsnildvmdees".upper()

kinase_seq["btk_ash"] = kinase_seq["btk_asp"]

kinase_seq["itk"] = "gkwvidpseltfvqeigsgqfglvhlgywlnkdkvaiktiregamseedfieeae\
vmmklshpklvqlygvcleqapiclvfefmehgclsdylrtqrglfaaetllgmcldvcegmayleeacvihrdlaa\
rnclvgenqvikvsdfgmtrfvlddqytsstgtkfpvkwaspevfsfsryssksdvwsfgvlmwevfsegkipyenrs\
nsevvedistgfrlykprlasthvyqimnhcwkerpedrpafsrllrqlaeiaesgl".upper()

kinase_seq["her2_wt"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikvlr\
entspkankeildeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiakgm\
syledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvwsygv\
tvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqrfv\
viqned".upper()

kinase_seq["her2_40x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeildeayvmagvgspgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcm\
qiakgmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqs\
dvwsygvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrm\
ardpqrfvviqned".upper()

kinase_seq["her2_22x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeilheayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqia\
kgmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvw\
sygvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmard\
pqrfvviqned".upper()

kinase_seq["her2_5x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeilheayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqia\
kgmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvw\
sygvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmard\
pqrfvviqned".upper()

kinase_seq["her2_4p7x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeilyeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiak\
gmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvwsy\
gvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqr\
fvviqned".upper()

kinase_seq["her2_3p3x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeildeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiak\
gmsyledvrlihrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvwsy\
gvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqr\
fvviqned".upper()

kinase_seq["her2_2p2x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeildeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiak\
gmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrcrfthqsdvwsy\
gvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqr\
fvviqne".upper()

kinase_seq["her2_1p2x"] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankemldeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiak\
gmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvwsy\
gvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqr\
fvviqned".upper()

kinase_seq['her2_1x'] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentapkankeildeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqia\
kgmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvw\
sygvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardp\
qrfvviqned".upper()

kinase_seq['her2_0p3x'] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
lrentspkankeildeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqiak\
gmsfledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvwsy\
gvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpqr\
fvviqned".upper()

kinase_seq['her2_lap'] = "ampnqaqmrilketelrkvkvlgsgafgtvykgiwipdgenvkipvaikv\
srentspkankeildeayvmagvgspyvsrllgicltstvqlvtqlmpygclldhvrenrgrlgsqdllnwcmqia\
kgmsyledvrlvhrdlaarnvlvkspnhvkitdfglarlldideteyhadggkvpikwmalesilrrrfthqsdvws\
ygvtvwelmtfgakpydgipareipdllekgerlpqppictidvymimvkcwmidsecrprfrelvsefsrmardpq\
rfvviqned".upper()

kinase_seq["sfk_fyn"] = "kdvweipreslqlikrlgngqfgevwmgtwngntkvaiktlkpgtmspesflee\
aqimkklkhdklvqlyavvseepiyivteymnkgslldflkdgegralklpnlvdmaaqvaagmayiermnyihrdlrsa\
nilvgnglickiadfglarliedneytarqgakfpikwtapeaalygrftiksdvwsfgilltelvtkgrvpypgmnnrev\
leqvergyrmpcpqdcpislhelmihcwkkdpeerptfeylqsfledy".upper()

kinase_seq["sfk_fgr"] = "kdaweisrssitlerrlgtgcfgdvwlgtwngstkvavktlkpgtmspkafleea\
qvmkllrhdklvqlyavvseepiyivtefmchgslldflknpegqdlrlpqlvdmaaqvaegmaymermnyihrdlraani\
lvgerlackiadfglarlikddeynpcqgskfpikwtapeaalfgrftiksdvwsfgilltelitkgripypgmnkrevle\
qveqgyhmpcppgcpaslyeameqtwrldpeerptfeylqsfledy".upper()

kinase_seq["sfk_hck"] = "kdaweipreslklekklgagqfgevwmatynkhtkvavktmkpgsmsveaflaea\
nvmktlqhdklvklhavvtkepiyiitefmakgslldflksdegskqplpklidfsaqiaegmafieqrnyihrdlraani\
lvsaslvckiadfglarviedneytaregakfpikwtapeainfgsftiksdvwsfgillmeivtygripypgmsnpevir\
alergyrmprpencpeelynimmrcwknrpeerptfeyiqsvlddf".upper()

kinase_seq["sfk_lck"] = "edewevpretlklverlgagqfgevwmgyynghtkvavkslkqgsmspdaflaea\
nlmkqlqhqrlvrlyavvtqepiyiiteymengslvdflktpsgikltinklldmaaqiaegmafieernyihrdlraani\
lvsdtlsckiadfglarliedneytaregakfpikwtapeainygtftiksdvwsfgillteivthgripypgmtnpeviq\
nlergyrmvrpdncpeelyqlmrlcwkerpedrptfdylrsvledf".upper()

kinase_seq["sfk_lyn"] = "kdaweipresiklvkrlgagqfgevwmgyynnstkvavktlkpgtmsvqafleean\
lmktlqhdklvrlyavvtreepiyiiteymakgslldflksdeggkvllpklidfsaqiaegmayierknyihrdlraanv\
lvseslmckiadfglarviedneytaregakfpikwtapeainfgcftiksdvwsfgillyeivtygkipypgrtnadvmta\
lsqgyrmprvencpdelydimkmcwkekaeerptfdylqsvlddf".upper()

kinase_seq["sfk_yes"] = "kdaweipreslrlevklgqgcfgevwmgtwngttkvaiktlkpgtmmpeaflqeaq\
imkklrhdklvplyavvseepiyivtefmskgslldflkegdgkylklpqlvdmaaqiadgmayiermnyihrdlraanilv\
genlvckiadfglarliedneytarqgakfpikwtapeaalygrftiksdvwsfgilqtelvtkgrvpypgmvnrevleqve\
rgyrmpcpqgcpeslhelmnlcwkkdpderptfeyiqsfledy".upper()

kinase_seq["sfk_blk"] = "qdeweiprqslrlvrklgsgqfgevwmgyyknnmkvaiktlkegtmspeaflgean\
vmkalqherlvrlyavvtkepiyivteymargclldflktdegsrlslprlidmsaqiaegmayiermnsihrdlraanil\
vsealcckiadfglariidseytaqegakfpikwtapeaihfgvftikadvwsfgvllmevvtygrvpypgmsnpevirnle\
rgyrmprpdtcppelyrgviaecwrsrpeerptfeflqsvledf".upper()

proj_list["her2_wt"] = [9114, 9120]
proj_list["her2_40x"] = [9104, 9121]
proj_list["her2_0p3x"] = [9112, 9119]
proj_list["her2_22x"] = [9105, 9122]
proj_list["her2_5x"] = [9106, 9123]
proj_list["her2_4p7x"] = [9107, 9124]
proj_list["her2_3p3x"] = [9108, 9125]
proj_list["her2_2p2x"] = [9109, 9126]
proj_list["her2_1p2x"] = [9110, 9127]
proj_list["her2_1x"] = [9111, 9128]
proj_list["her2_lap"] = [9113, 9129]
proj_list["btk_asp"] = [9145, 9151]
proj_list["btk_ash"] = [9146, 9152]
