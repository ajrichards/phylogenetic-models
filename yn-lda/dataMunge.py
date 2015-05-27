#!/usr/bin/env python
"""
get data into a format that LDA understands
https://pythonhosted.org/ete2/tutorial/tutorial_trees.html#reading-newick-trees
"""

import os,re,sys

dirName = "vieira-high"
parse = "vieira-high5-alpha_sine_endo-dc_codon.phylip-16x1338-v1.5a-ms-c-uni-"

for fileName in os.listdir(dirName):
    if re.match(parse,fileName):
        os.rename(os.path.join(dirName,fileName), os.path.join(dirName,re.sub(parse,"",fileName)))

dirName = "vieira-low"
parse = "vieira-low5-alpha_sine_endo-dc_codon.phylip-16x1649-v1.5a-ms-c-uni-"    
for fileName in os.listdir(dirName):
    if re.match(parse,fileName):
        os.rename(os.path.join(dirName,fileName), os.path.join(dirName,re.sub(parse,"",fileName)))



## debug
constrained = "(Rhodospirillum_rubrum_3.8_ACC:1.21718:ACC:0.795739:ATC:0.423212:GTC,Gluconobacter_oxydans_621H_0.94_ACC:2.74115:ACC:0.799534:AAC:0.627617:ACC:0.146405:ATC:1.33201:GTC,((Zymomonas_mobilis_2_ATT:1.09065:ATT:0.849922:GTT:2.71093:GTC:0.449868:CTC:0.131212:GTC,Sphingopyxis_alaskensis_2.39_ATC:0.74794:ATC:0.594103:GTC)_GTC:1.0974:GTC:0.38102:TTC:0.913789:GTC:0.0701186:GCC,((Rhodobacter_sphaeroides_3_ACC:1.08263:ACC,Paracoccus_denitrificans_2.1_ACC:1.30088:ACC)_ACC:1.38045:ACC:0.18737:GCC,(Caulobacter_crescentus_1.5_TAT:0.238576:TAT:1.60043:TTT:0.695896:TTC:0.327243:TCC:0.0261362:GCC,((Nitrobacter_winogradskyi_8_GTC:1.1677:GTC,(Rhodopseudomonas_palustris_9_ATC:0.510936:ATC:0.392843:GTC,Bradyrhizobium_japonicum_USDA110_20_GTC:0.691677:GTC)_GTC:0.0708203:GTC)_GTC:1.28532:GTC,((Mesorhizobium_loti_2.4_GTC:1.38965:GTC,(Brucella_suis_1330_2_ATC:0.427971:ATC,(Bartonella_quintana_Toulouse_3_ATT:0.486504:ATT,Bartonella_henselae_Houston-1_3_ATT:0.363287:ATT)_ATT:0.68998:ATT:0.184105:ACT:0.0937438:ACG:0.672405:ACC:0.682118:ACT:0.630906:ACC:0.647151:ACT:0.8754:ATT:0.498107:ATC)_ATC:0.766973:ATC:0.629549:GTC)_GTC:0.292427:GTC,(Sinorhizobium_meliloti_1021_1.5_GTC:0.998178:GTC,Agrobacterium_tumefaciens_C58_3_GTC:1.33291:GTC)_GTC:0.73082:GTC)_GTC:1.26457:GTC)_GTC:0.506179:GTC:0.365595:GCC)_GCC:0.0466269:GCC)_GCC:0.0418115:GCC)_GCC:0.443192:GCC:0.196449:GTC)_GTC;"

unconstrained = "(Rhodospirillum_rubrum_3.8_GTC:0.70614:GTC:0.895494:ATC:0.214485:ACC:0.620016:ATC,Gluconobacter_oxydans_621H_0.94_ATC:5.64671:ATC,((Zymomonas_mobilis_2_TTC:1.77575:TTC:1.81786:TTT:0.468152:ATT:1.17081:ATC,Sphingopyxis_alaskensis_2.39_ATT:0.302737:ATT:0.409838:ATC:0.520213:TTC:0.109255:ATC)_ATC:2.46233:ATC,((Rhodobacter_sphaeroides_3_ACC:1.08263:ACC,Paracoccus_denitrificans_2.1_ACC:1.30088:ACC)_ACC:1.44124:ACC:0.126574:ATC,(Caulobacter_crescentus_1.5_ACC:2.53826:ACC:0.350027:ATC,((Nitrobacter_winogradskyi_8_GAA:0.442367:GAA:0.144013:GAG:0.0699881:GAC:0.511332:GTC,(Rhodopseudomonas_palustris_9_GTT:0.6072:GTT:0.296579:GTC,Bradyrhizobium_japonicum_USDA110_20_GTC:0.691677:GTC)_GTC:0.0708203:GTC)_GTC:0.515999:GTC:0.671378:GTG:0.0979445:GTC,((Mesorhizobium_loti_2.4_GTC:1.38965:GTC,(Brucella_suis_1330_2_ATC:0.427971:ATC,(Bartonella_quintana_Toulouse_3_ATT:0.486504:ATT,Bartonella_henselae_Houston-1_3_ATT:0.363287:ATT)_ATT:3.81122:ATT:1.1627:ATC)_ATC:0.845855:ATC:0.550666:GTC)_GTC:0.292427:GTC,(Sinorhizobium_meliloti_1021_1.5_GTC:0.998178:GTC,Agrobacterium_tumefaciens_C58_3_GTC:1.33291:GTC)_GTC:0.73082:GTC)_GTC:1.26457:GTC)_GTC:0.141463:GTC:0.730312:ATC)_ATC:0.0466269:ATC)_ATC:0.0418115:ATC)_ATC:0.63964:ATC)_ATC;"

from ete2 import Tree
t = Tree(constrained)

#print dir(t)
parents = []
level = 0
for node in t.traverse("postorder"):
    #if level == 0 and node.name != 'NoName':
    #    parents.append(node.name)
    print node.name

    #for a in node.children:
    #    print '...',a.name
t.show()
