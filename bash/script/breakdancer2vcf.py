#!/usr/bin/env python
# -*- coding: utf-8 -*-
# converts breakdancer output into pseudo-VCF

import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Convert BreakDancer output to VCF format')
    parser.add_argument('-i', '--input', required=True, help='Input BreakDancer file')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    parser.add_argument('-s', '--sample', required=True, help='sample')
    return parser.parse_args()

def write_vcf_header(output, sample):
    header = """##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##INFO=<ID=Chr1,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Pos1,Number=1,Type=Integer,Description="BreakDancer output">
##INFO=<ID=Orient1,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Chr2,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Pos2,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Orient2,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Type,Number=1,Type=String,Description="BreakDancer output">
##INFO=<ID=Size,Number=1,Type=Integer,Description="BreakDancer output">
##INFO=<ID=Score,Number=1,Type=Float,Description="BreakDancer output">
##INFO=<ID=num_Reads,Number=1,Type=Integer,Description="BreakDancer output">
##contig=<ID=1,length=274330532>
##contig=<ID=2,length=151935994>
##contig=<ID=3,length=132848913>
##contig=<ID=4,length=130910915>
##contig=<ID=5,length=104526007>
##contig=<ID=6,length=170843587>
##contig=<ID=7,length=121844099>
##contig=<ID=8,length=138966237>
##contig=<ID=9,length=139512083>
##contig=<ID=10,length=69359453>
##contig=<ID=11,length=79169978>
##contig=<ID=12,length=61602749>
##contig=<ID=13,length=208334590>
##contig=<ID=14,length=141755446>
##contig=<ID=15,length=140412725>
##contig=<ID=16,length=79944280>
##contig=<ID=17,length=63494081>
##contig=<ID=18,length=55982971>
##contig=<ID=X,length=125939595>
##contig=<ID=Y,length=43547828>
##contig=<ID=AEMK02000452.1,length=3843259>
##contig=<ID=AEMK02000698.1,length=2641821>
##contig=<ID=AEMK02000361.1,length=2125600>
##contig=<ID=AEMK02000598.1,length=2117578>
##contig=<ID=AEMK02000682.1,length=2088057>
##contig=<ID=AEMK02000449.1,length=1437513>
##contig=<ID=AEMK02000677.1,length=1284139>
##contig=<ID=AEMK02000408.1,length=1187208>
##contig=<ID=AEMK02000569.1,length=1159966>
##contig=<ID=AEMK02000393.1,length=914028>
##contig=<ID=AEMK02000510.1,length=899810>
##contig=<ID=AEMK02000328.1,length=734522>
##contig=<ID=AEMK02000465.1,length=725547>
##contig=<ID=AEMK02000133.1,length=663541>
##contig=<ID=AEMK02000661.1,length=658525>
##contig=<ID=AEMK02000261.1,length=631996>
##contig=<ID=FPKY02000002.1,length=492890>
##contig=<ID=AEMK02000324.1,length=463454>
##contig=<ID=AEMK02000258.1,length=463424>
##contig=<ID=AEMK02000271.1,length=457951>
##contig=<ID=AEMK02000537.1,length=452860>
##contig=<ID=AEMK02000522.1,length=450457>
##contig=<ID=AEMK02000639.1,length=441855>
##contig=<ID=FPKY02000001.1,length=438769>
##contig=<ID=AEMK02000694.1,length=437940>
##contig=<ID=AEMK02000289.1,length=437422>
##contig=<ID=AEMK02000451.1,length=414619>
##contig=<ID=AEMK02000478.1,length=413771>
##contig=<ID=AEMK02000496.1,length=404434>
##contig=<ID=AEMK02000390.1,length=383558>
##contig=<ID=AEMK02000175.1,length=368575>
##contig=<ID=AEMK02000602.1,length=362670>
##contig=<ID=AEMK02000370.1,length=353404>
##contig=<ID=AEMK02000213.1,length=352148>
##contig=<ID=AEMK02000630.1,length=338089>
##contig=<ID=AEMK02000264.1,length=333895>
##contig=<ID=AEMK02000600.1,length=332772>
##contig=<ID=AEMK02000410.1,length=328665>
##contig=<ID=AEMK02000525.1,length=305256>
##contig=<ID=AEMK02000159.1,length=304882>
##contig=<ID=AEMK02000423.1,length=303938>
##contig=<ID=AEMK02000566.1,length=292505>
##contig=<ID=AEMK02000278.1,length=250081>
##contig=<ID=AEMK02000577.1,length=248197>
##contig=<ID=AEMK02000292.1,length=241003>
##contig=<ID=AEMK02000697.1,length=238141>
##contig=<ID=AEMK02000532.1,length=234155>
##contig=<ID=AEMK02000189.1,length=231906>
##contig=<ID=AEMK02000400.1,length=230773>
##contig=<ID=AEMK02000253.1,length=229168>
##contig=<ID=AEMK02000649.1,length=218479>
##contig=<ID=AEMK02000256.1,length=196394>
##contig=<ID=AEMK02000234.1,length=195359>
##contig=<ID=AEMK02000368.1,length=192446>
##contig=<ID=AEMK02000640.1,length=191720>
##contig=<ID=AEMK02000703.1,length=186897>
##contig=<ID=AEMK02000341.1,length=185816>
##contig=<ID=AEMK02000495.1,length=185666>
##contig=<ID=AEMK02000497.1,length=182797>
##contig=<ID=AEMK02000310.1,length=177715>
##contig=<ID=AEMK02000297.1,length=177295>
##contig=<ID=AEMK02000171.1,length=176234>
##contig=<ID=AEMK02000517.1,length=175079>
##contig=<ID=AEMK02000286.1,length=172180>
##contig=<ID=AEMK02000425.1,length=170852>
##contig=<ID=AEMK02000591.1,length=169331>
##contig=<ID=AEMK02000485.1,length=168947>
##contig=<ID=FPKY02000005.1,length=167521>
##contig=<ID=AEMK02000658.1,length=166067>
##contig=<ID=FPKY02000004.1,length=162579>
##contig=<ID=AEMK02000435.1,length=162163>
##contig=<ID=AEMK02000215.1,length=161540>
##contig=<ID=AEMK02000589.1,length=159309>
##contig=<ID=AEMK02000204.1,length=158593>
##contig=<ID=AEMK02000612.1,length=156208>
##contig=<ID=AEMK02000493.1,length=154740>
##contig=<ID=AEMK02000579.1,length=153825>
##contig=<ID=AEMK02000476.1,length=153176>
##contig=<ID=AEMK02000301.1,length=153020>
##contig=<ID=AEMK02000394.1,length=152378>
##contig=<ID=AEMK02000434.1,length=151582>
##contig=<ID=AEMK02000477.1,length=151529>
##contig=<ID=AEMK02000156.1,length=151356>
##contig=<ID=AEMK02000526.1,length=151330>
##contig=<ID=AEMK02000468.1,length=149177>
##contig=<ID=AEMK02000336.1,length=148079>
##contig=<ID=AEMK02000660.1,length=147927>
##contig=<ID=AEMK02000240.1,length=145869>
##contig=<ID=AEMK02000320.1,length=145609>
##contig=<ID=AEMK02000162.1,length=141821>
##contig=<ID=AEMK02000444.1,length=136983>
##contig=<ID=AEMK02000333.1,length=136362>
##contig=<ID=AEMK02000680.1,length=134991>
##contig=<ID=AEMK02000248.1,length=133611>
##contig=<ID=AEMK02000254.1,length=131902>
##contig=<ID=AEMK02000173.1,length=131066>
##contig=<ID=AEMK02000453.1,length=130110>
##contig=<ID=AEMK02000617.1,length=128219>
##contig=<ID=AEMK02000431.1,length=127758>
##contig=<ID=AEMK02000576.1,length=127286>
##contig=<ID=AEMK02000620.1,length=126623>
##contig=<ID=AEMK02000378.1,length=126168>
##contig=<ID=AEMK02000276.1,length=125997>
##contig=<ID=AEMK02000137.1,length=125348>
##contig=<ID=AEMK02000514.1,length=122898>
##contig=<ID=AEMK02000679.1,length=122690>
##contig=<ID=AEMK02000277.1,length=122329>
##contig=<ID=AEMK02000238.1,length=121275>
##contig=<ID=AEMK02000318.1,length=121216>
##contig=<ID=AEMK02000588.1,length=116223>
##contig=<ID=AEMK02000399.1,length=114674>
##contig=<ID=AEMK02000200.1,length=113982>
##contig=<ID=AEMK02000182.1,length=111963>
##contig=<ID=AEMK02000383.1,length=111191>
##contig=<ID=AEMK02000558.1,length=111136>
##contig=<ID=AEMK02000692.1,length=110875>
##contig=<ID=AEMK02000295.1,length=107759>
##contig=<ID=AEMK02000678.1,length=107223>
##contig=<ID=AEMK02000584.1,length=106523>
##contig=<ID=AEMK02000131.1,length=105554>
##contig=<ID=AEMK02000186.1,length=105175>
##contig=<ID=AEMK02000364.1,length=104663>
##contig=<ID=AEMK02000696.1,length=101917>
##contig=<ID=AEMK02000438.1,length=100920>
##contig=<ID=AEMK02000534.1,length=100697>
##contig=<ID=AEMK02000299.1,length=100207>
##contig=<ID=AEMK02000621.1,length=99288>
##contig=<ID=AEMK02000440.1,length=99018>
##contig=<ID=AEMK02000528.1,length=98791>
##contig=<ID=AEMK02000494.1,length=98741>
##contig=<ID=AEMK02000623.1,length=98239>
##contig=<ID=AEMK02000585.1,length=97889>
##contig=<ID=AEMK02000369.1,length=97098>
##contig=<ID=AEMK02000245.1,length=96200>
##contig=<ID=AEMK02000554.1,length=95954>
##contig=<ID=AEMK02000541.1,length=95750>
##contig=<ID=AEMK02000665.1,length=95440>
##contig=<ID=AEMK02000593.1,length=95045>
##contig=<ID=AEMK02000574.1,length=94815>
##contig=<ID=AEMK02000632.1,length=93881>
##contig=<ID=AEMK02000503.1,length=92946>
##contig=<ID=AEMK02000637.1,length=92764>
##contig=<ID=AEMK02000223.1,length=92180>
##contig=<ID=AEMK02000460.1,length=92164>
##contig=<ID=AEMK02000672.1,length=92067>
##contig=<ID=AEMK02000152.1,length=91927>
##contig=<ID=AEMK02000482.1,length=91570>
##contig=<ID=AEMK02000472.1,length=91429>
##contig=<ID=AEMK02000208.1,length=91308>
##contig=<ID=AEMK02000346.1,length=89764>
##contig=<ID=AEMK02000407.1,length=89753>
##contig=<ID=AEMK02000164.1,length=89749>
##contig=<ID=AEMK02000153.1,length=87954>
##contig=<ID=AEMK02000456.1,length=87638>
##contig=<ID=AEMK02000170.1,length=87541>
##contig=<ID=AEMK02000191.1,length=86406>
##contig=<ID=AEMK02000662.1,length=85322>
##contig=<ID=AEMK02000572.1,length=85184>
##contig=<ID=AEMK02000210.1,length=84782>
##contig=<ID=AEMK02000437.1,length=83739>
##contig=<ID=AEMK02000291.1,length=83331>
##contig=<ID=AEMK02000563.1,length=83258>
##contig=<ID=AEMK02000644.1,length=83016>
##contig=<ID=AEMK02000507.1,length=82813>
##contig=<ID=AEMK02000415.1,length=81645>
##contig=<ID=FPKY02000007.1,length=81402>
##contig=<ID=AEMK02000250.1,length=80194>
##contig=<ID=AEMK02000513.1,length=80124>
##contig=<ID=AEMK02000515.1,length=79602>
##contig=<ID=AEMK02000422.1,length=78723>
##contig=<ID=AEMK02000474.1,length=78397>
##contig=<ID=AEMK02000168.1,length=77398>
##contig=<ID=AEMK02000685.1,length=77218>
##contig=<ID=AEMK02000683.1,length=75957>
##contig=<ID=AEMK02000512.1,length=75122>
##contig=<ID=AEMK02000389.1,length=75082>
##contig=<ID=AEMK02000568.1,length=74675>
##contig=<ID=AEMK02000559.1,length=74153>
##contig=<ID=AEMK02000366.1,length=74125>
##contig=<ID=AEMK02000656.1,length=73391>
##contig=<ID=AEMK02000235.1,length=73291>
##contig=<ID=FPKY02000006.1,length=73062>
##contig=<ID=AEMK02000227.1,length=72950>
##contig=<ID=AEMK02000287.1,length=72946>
##contig=<ID=AEMK02000427.1,length=72884>
##contig=<ID=AEMK02000251.1,length=72613>
##contig=<ID=AEMK02000616.1,length=72435>
##contig=<ID=AEMK02000699.1,length=72425>
##contig=<ID=AEMK02000704.1,length=72037>
##contig=<ID=AEMK02000188.1,length=72014>
##contig=<ID=AEMK02000506.1,length=71707>
##contig=<ID=AEMK02000655.1,length=71706>
##contig=<ID=AEMK02000529.1,length=71279>
##contig=<ID=AEMK02000414.1,length=71268>
##contig=<ID=AEMK02000232.1,length=71095>
##contig=<ID=AEMK02000331.1,length=70479>
##contig=<ID=AEMK02000646.1,length=69502>
##contig=<ID=AEMK02000172.1,length=68315>
##contig=<ID=AEMK02000146.1,length=67873>
##contig=<ID=AEMK02000207.1,length=67854>
##contig=<ID=AEMK02000332.1,length=67596>
##contig=<ID=AEMK02000421.1,length=67144>
##contig=<ID=AEMK02000521.1,length=66883>
##contig=<ID=AEMK02000221.1,length=66807>
##contig=<ID=AEMK02000225.1,length=66683>
##contig=<ID=AEMK02000391.1,length=66671>
##contig=<ID=AEMK02000163.1,length=66393>
##contig=<ID=AEMK02000300.1,length=65977>
##contig=<ID=AEMK02000628.1,length=65930>
##contig=<ID=AEMK02000668.1,length=65664>
##contig=<ID=AEMK02000380.1,length=65658>
##contig=<ID=AEMK02000177.1,length=65518>
##contig=<ID=AEMK02000259.1,length=65455>
##contig=<ID=AEMK02000334.1,length=65042>
##contig=<ID=AEMK02000203.1,length=64878>
##contig=<ID=AEMK02000501.1,length=64639>
##contig=<ID=AEMK02000641.1,length=64268>
##contig=<ID=AEMK02000344.1,length=64041>
##contig=<ID=AEMK02000398.1,length=63919>
##contig=<ID=AEMK02000263.1,length=63916>
##contig=<ID=AEMK02000260.1,length=63393>
##contig=<ID=AEMK02000273.1,length=63293>
##contig=<ID=AEMK02000375.1,length=63078>
##contig=<ID=FPKY02000008.1,length=62281>
##contig=<ID=AEMK02000311.1,length=61746>
##contig=<ID=AEMK02000327.1,length=61624>
##contig=<ID=AEMK02000412.1,length=61586>
##contig=<ID=AEMK02000305.1,length=61450>
##contig=<ID=AEMK02000229.1,length=61233>
##contig=<ID=AEMK02000626.1,length=61233>
##contig=<ID=AEMK02000636.1,length=61061>
##contig=<ID=AEMK02000676.1,length=60710>
##contig=<ID=AEMK02000247.1,length=60458>
##contig=<ID=AEMK02000135.1,length=60359>
##contig=<ID=AEMK02000382.1,length=60335>
##contig=<ID=AEMK02000561.1,length=59899>
##contig=<ID=AEMK02000659.1,length=59742>
##contig=<ID=AEMK02000445.1,length=59679>
##contig=<ID=AEMK02000571.1,length=59490>
##contig=<ID=AEMK02000627.1,length=59096>
##contig=<ID=AEMK02000595.1,length=59052>
##contig=<ID=AEMK02000280.1,length=58992>
##contig=<ID=AEMK02000592.1,length=58548>
##contig=<ID=AEMK02000169.1,length=58416>
##contig=<ID=AEMK02000339.1,length=58201>
##contig=<ID=AEMK02000275.1,length=58047>
##contig=<ID=AEMK02000631.1,length=57870>
##contig=<ID=AEMK02000363.1,length=57548>
##contig=<ID=AEMK02000509.1,length=57178>
##contig=<ID=AEMK02000647.1,length=56960>
##contig=<ID=AEMK02000179.1,length=56626>
##contig=<ID=FPKY02000003.1,length=56573>
##contig=<ID=AEMK02000439.1,length=56510>
##contig=<ID=AEMK02000403.1,length=56236>
##contig=<ID=AEMK02000183.1,length=55954>
##contig=<ID=AEMK02000542.1,length=55808>
##contig=<ID=AEMK02000304.1,length=55662>
##contig=<ID=AEMK02000309.1,length=55519>
##contig=<ID=AEMK02000130.1,length=55490>
##contig=<ID=AEMK02000315.1,length=55444>
##contig=<ID=AEMK02000282.1,length=55403>
##contig=<ID=AEMK02000543.1,length=55144>
##contig=<ID=AEMK02000673.1,length=55089>
##contig=<ID=AEMK02000611.1,length=55069>
##contig=<ID=AEMK02000174.1,length=54656>
##contig=<ID=AEMK02000176.1,length=54489>
##contig=<ID=AEMK02000523.1,length=54345>
##contig=<ID=AEMK02000142.1,length=54115>
##contig=<ID=AEMK02000228.1,length=53671>
##contig=<ID=AEMK02000606.1,length=53195>
##contig=<ID=AEMK02000674.1,length=52974>
##contig=<ID=AEMK02000443.1,length=52466>
##contig=<ID=AEMK02000340.1,length=52374>
##contig=<ID=AEMK02000461.1,length=52294>
##contig=<ID=AEMK02000518.1,length=52274>
##contig=<ID=AEMK02000550.1,length=52058>
##contig=<ID=AEMK02000489.1,length=51932>
##contig=<ID=AEMK02000274.1,length=51665>
##contig=<ID=AEMK02000480.1,length=51625>
##contig=<ID=AEMK02000313.1,length=50759>
##contig=<ID=AEMK02000125.1,length=50575>
##contig=<ID=AEMK02000385.1,length=50567>
##contig=<ID=AEMK02000583.1,length=50335>
##contig=<ID=AEMK02000702.1,length=50243>
##contig=<ID=AEMK02000241.1,length=50081>
##contig=<ID=AEMK02000181.1,length=49406>
##contig=<ID=AEMK02000197.1,length=49329>
##contig=<ID=AEMK02000416.1,length=49126>
##contig=<ID=AEMK02000471.1,length=49050>
##contig=<ID=AEMK02000283.1,length=48980>
##contig=<ID=AEMK02000457.1,length=48413>
##contig=<ID=AEMK02000426.1,length=47889>
##contig=<ID=AEMK02000144.1,length=47643>
##contig=<ID=AEMK02000404.1,length=47524>
##contig=<ID=AEMK02000552.1,length=47289>
##contig=<ID=AEMK02000129.1,length=47123>
##contig=<ID=AEMK02000147.1,length=47030>
##contig=<ID=AEMK02000575.1,length=46903>
##contig=<ID=AEMK02000138.1,length=46849>
##contig=<ID=AEMK02000395.1,length=46845>
##contig=<ID=AEMK02000345.1,length=46301>
##contig=<ID=AEMK02000549.1,length=46115>
##contig=<ID=AEMK02000483.1,length=46114>
##contig=<ID=AEMK02000322.1,length=46061>
##contig=<ID=AEMK02000657.1,length=46053>
##contig=<ID=AEMK02000255.1,length=45734>
##contig=<ID=AEMK02000511.1,length=45692>
##contig=<ID=AEMK02000161.1,length=45522>
##contig=<ID=AEMK02000308.1,length=45474>
##contig=<ID=AEMK02000540.1,length=45163>
##contig=<ID=AEMK02000266.1,length=44947>
##contig=<ID=AEMK02000429.1,length=44863>
##contig=<ID=AEMK02000596.1,length=44774>
##contig=<ID=AEMK02000386.1,length=44543>
##contig=<ID=FPKY02000009.1,length=44425>
##contig=<ID=AEMK02000622.1,length=44010>
##contig=<ID=AEMK02000614.1,length=43888>
##contig=<ID=AEMK02000565.1,length=43773>
##contig=<ID=AEMK02000279.1,length=43530>
##contig=<ID=AEMK02000555.1,length=43188>
##contig=<ID=AEMK02000367.1,length=43143>
##contig=<ID=AEMK02000590.1,length=42628>
##contig=<ID=AEMK02000504.1,length=42619>
##contig=<ID=AEMK02000373.1,length=42492>
##contig=<ID=AEMK02000546.1,length=42190>
##contig=<ID=AEMK02000687.1,length=42170>
##contig=<ID=AEMK02000479.1,length=42029>
##contig=<ID=AEMK02000180.1,length=42006>
##contig=<ID=AEMK02000411.1,length=41893>
##contig=<ID=AEMK02000392.1,length=41701>
##contig=<ID=AEMK02000330.1,length=41407>
##contig=<ID=AEMK02000487.1,length=41262>
##contig=<ID=AEMK02000226.1,length=40713>
##contig=<ID=AEMK02000466.1,length=40414>
##contig=<ID=AEMK02000556.1,length=40390>
##contig=<ID=AEMK02000402.1,length=40091>
##contig=<ID=AEMK02000681.1,length=40048>
##contig=<ID=AEMK02000520.1,length=39997>
##contig=<ID=AEMK02000202.1,length=39923>
##contig=<ID=AEMK02000268.1,length=39828>
##contig=<ID=AEMK02000187.1,length=39765>
##contig=<ID=AEMK02000325.1,length=39598>
##contig=<ID=AEMK02000433.1,length=39522>
##contig=<ID=AEMK02000140.1,length=39503>
##contig=<ID=AEMK02000638.1,length=39448>
##contig=<ID=AEMK02000671.1,length=39267>
##contig=<ID=AEMK02000619.1,length=39236>
##contig=<ID=AEMK02000469.1,length=39098>
##contig=<ID=AEMK02000244.1,length=39063>
##contig=<ID=AEMK02000209.1,length=38942>
##contig=<ID=AEMK02000374.1,length=38887>
##contig=<ID=AEMK02000222.1,length=38859>
##contig=<ID=AEMK02000139.1,length=38775>
##contig=<ID=AEMK02000548.1,length=38669>
##contig=<ID=AEMK02000246.1,length=38453>
##contig=<ID=AEMK02000505.1,length=38289>
##contig=<ID=AEMK02000624.1,length=38262>
##contig=<ID=AEMK02000645.1,length=38247>
##contig=<ID=AEMK02000442.1,length=38108>
##contig=<ID=AEMK02000675.1,length=38099>
##contig=<ID=AEMK02000252.1,length=38060>
##contig=<ID=AEMK02000536.1,length=37985>
##contig=<ID=AEMK02000500.1,length=37867>
##contig=<ID=AEMK02000243.1,length=37675>
##contig=<ID=AEMK02000267.1,length=37674>
##contig=<ID=AEMK02000450.1,length=37546>
##contig=<ID=AEMK02000127.1,length=37414>
##contig=<ID=AEMK02000581.1,length=37293>
##contig=<ID=AEMK02000150.1,length=37252>
##contig=<ID=AEMK02000365.1,length=37050>
##contig=<ID=AEMK02000459.1,length=36996>
##contig=<ID=AEMK02000420.1,length=36859>
##contig=<ID=AEMK02000467.1,length=36784>
##contig=<ID=AEMK02000603.1,length=36608>
##contig=<ID=AEMK02000570.1,length=36425>
##contig=<ID=AEMK02000664.1,length=36264>
##contig=<ID=AEMK02000269.1,length=36105>
##contig=<ID=AEMK02000607.1,length=36098>
##contig=<ID=AEMK02000387.1,length=35978>
##contig=<ID=AEMK02000201.1,length=35839>
##contig=<ID=AEMK02000377.1,length=35800>
##contig=<ID=AEMK02000618.1,length=35796>
##contig=<ID=AEMK02000157.1,length=35768>
##contig=<ID=AEMK02000355.1,length=35698>
##contig=<ID=AEMK02000597.1,length=35610>
##contig=<ID=AEMK02000216.1,length=35435>
##contig=<ID=AEMK02000337.1,length=35413>
##contig=<ID=AEMK02000475.1,length=35399>
##contig=<ID=AEMK02000165.1,length=35256>
##contig=<ID=AEMK02000134.1,length=35132>
##contig=<ID=AEMK02000218.1,length=35090>
##contig=<ID=AEMK02000379.1,length=35033>
##contig=<ID=AEMK02000553.1,length=34945>
##contig=<ID=AEMK02000604.1,length=34712>
##contig=<ID=AEMK02000635.1,length=34704>
##contig=<ID=AEMK02000316.1,length=34675>
##contig=<ID=AEMK02000533.1,length=34572>
##contig=<ID=AEMK02000233.1,length=34570>
##contig=<ID=AEMK02000356.1,length=34299>
##contig=<ID=AEMK02000360.1,length=34294>
##contig=<ID=AEMK02000348.1,length=34187>
##contig=<ID=AEMK02000231.1,length=33498>
##contig=<ID=AEMK02000605.1,length=33397>
##contig=<ID=AEMK02000643.1,length=33108>
##contig=<ID=AEMK02000199.1,length=32990>
##contig=<ID=AEMK02000490.1,length=32989>
##contig=<ID=AEMK02000192.1,length=32780>
##contig=<ID=AEMK02000265.1,length=32772>
##contig=<ID=AEMK02000650.1,length=32614>
##contig=<ID=AEMK02000257.1,length=32545>
##contig=<ID=AEMK02000613.1,length=32495>
##contig=<ID=AEMK02000634.1,length=32492>
##contig=<ID=AEMK02000473.1,length=32348>
##contig=<ID=AEMK02000582.1,length=32249>
##contig=<ID=AEMK02000524.1,length=32207>
##contig=<ID=AEMK02000262.1,length=32067>
##contig=<ID=AEMK02000448.1,length=32062>
##contig=<ID=AEMK02000123.1,length=32000>
##contig=<ID=AEMK02000578.1,length=31967>
##contig=<ID=AEMK02000249.1,length=31957>
##contig=<ID=AEMK02000303.1,length=31897>
##contig=<ID=AEMK02000237.1,length=31858>
##contig=<ID=AEMK02000319.1,length=31818>
##contig=<ID=AEMK02000609.1,length=31792>
##contig=<ID=AEMK02000211.1,length=31659>
##contig=<ID=AEMK02000498.1,length=31599>
##contig=<ID=AEMK02000516.1,length=31562>
##contig=<ID=AEMK02000141.1,length=31401>
##contig=<ID=AEMK02000143.1,length=31175>
##contig=<ID=AEMK02000326.1,length=30837>
##contig=<ID=AEMK02000499.1,length=30824>
##contig=<ID=AEMK02000293.1,length=30558>
##contig=<ID=AEMK02000418.1,length=30547>
##contig=<ID=AEMK02000610.1,length=30466>
##contig=<ID=AEMK02000376.1,length=30397>
##contig=<ID=AEMK02000347.1,length=30289>
##contig=<ID=AEMK02000314.1,length=30249>
##contig=<ID=AEMK02000560.1,length=30120>
##contig=<ID=AEMK02000601.1,length=30048>
##contig=<ID=AEMK02000455.1,length=30003>
##contig=<ID=AEMK02000372.1,length=29770>
##contig=<ID=AEMK02000484.1,length=29734>
##contig=<ID=AEMK02000329.1,length=29544>
##contig=<ID=AEMK02000359.1,length=29315>
##contig=<ID=AEMK02000470.1,length=29154>
##contig=<ID=AEMK02000160.1,length=29131>
##contig=<ID=AEMK02000428.1,length=28902>
##contig=<ID=AEMK02000195.1,length=28777>
##contig=<ID=AEMK02000288.1,length=28312>
##contig=<ID=AEMK02000463.1,length=28245>
##contig=<ID=AEMK02000298.1,length=28220>
##contig=<ID=AEMK02000436.1,length=28062>
##contig=<ID=AEMK02000405.1,length=28054>
##contig=<ID=AEMK02000538.1,length=27978>
##contig=<ID=AEMK02000321.1,length=27922>
##contig=<ID=AEMK02000633.1,length=27743>
##contig=<ID=AEMK02000158.1,length=27693>
##contig=<ID=AEMK02000196.1,length=27681>
##contig=<ID=AEMK02000388.1,length=27663>
##contig=<ID=AEMK02000488.1,length=27635>
##contig=<ID=AEMK02000615.1,length=27455>
##contig=<ID=AEMK02000224.1,length=27355>
##contig=<ID=AEMK02000666.1,length=27124>
##contig=<ID=AEMK02000290.1,length=26999>
##contig=<ID=AEMK02000648.1,length=26953>
##contig=<ID=AEMK02000629.1,length=26916>
##contig=<ID=AEMK02000691.1,length=26904>
##contig=<ID=AEMK02000317.1,length=26654>
##contig=<ID=AEMK02000132.1,length=26642>
##contig=<ID=AEMK02000126.1,length=26602>
##contig=<ID=AEMK02000599.1,length=26507>
##contig=<ID=AEMK02000242.1,length=26351>
##contig=<ID=AEMK02000446.1,length=26304>
##contig=<ID=AEMK02000527.1,length=26146>
##contig=<ID=AEMK02000684.1,length=26145>
##contig=<ID=AEMK02000154.1,length=26018>
##contig=<ID=AEMK02000212.1,length=25912>
##contig=<ID=AEMK02000447.1,length=25908>
##contig=<ID=AEMK02000136.1,length=25867>
##contig=<ID=AEMK02000670.1,length=25784>
##contig=<ID=AEMK02000178.1,length=25711>
##contig=<ID=AEMK02000653.1,length=25690>
##contig=<ID=AEMK02000625.1,length=25476>
##contig=<ID=AEMK02000194.1,length=25467>
##contig=<ID=AEMK02000167.1,length=25463>
##contig=<ID=AEMK02000338.1,length=25373>
##contig=<ID=AEMK02000193.1,length=25358>
##contig=<ID=AEMK02000454.1,length=25241>
##contig=<ID=AEMK02000567.1,length=25217>
##contig=<ID=AEMK02000128.1,length=25182>
##contig=<ID=AEMK02000312.1,length=25057>
##contig=<ID=AEMK02000335.1,length=24880>
##contig=<ID=AEMK02000667.1,length=24873>
##contig=<ID=AEMK02000354.1,length=24819>
##contig=<ID=AEMK02000530.1,length=24729>
##contig=<ID=AEMK02000198.1,length=24675>
##contig=<ID=AEMK02000688.1,length=24670>
##contig=<ID=AEMK02000587.1,length=24651>
##contig=<ID=AEMK02000371.1,length=24648>
##contig=<ID=AEMK02000413.1,length=24645>
##contig=<ID=AEMK02000296.1,length=24616>
##contig=<ID=AEMK02000547.1,length=24575>
##contig=<ID=AEMK02000545.1,length=24571>
##contig=<ID=AEMK02000464.1,length=24501>
##contig=<ID=AEMK02000396.1,length=24422>
##contig=<ID=AEMK02000669.1,length=24413>
##contig=<ID=AEMK02000419.1,length=24396>
##contig=<ID=AEMK02000652.1,length=24338>
##contig=<ID=AEMK02000205.1,length=24313>
##contig=<ID=AEMK02000693.1,length=24178>
##contig=<ID=AEMK02000557.1,length=24118>
##contig=<ID=AEMK02000149.1,length=24060>
##contig=<ID=AEMK02000353.1,length=24015>
##contig=<ID=AEMK02000651.1,length=23994>
##contig=<ID=AEMK02000352.1,length=23971>
##contig=<ID=AEMK02000166.1,length=23908>
##contig=<ID=AEMK02000145.1,length=23883>
##contig=<ID=AEMK02000700.1,length=23737>
##contig=<ID=AEMK02000306.1,length=23710>
##contig=<ID=AEMK02000462.1,length=23690>
##contig=<ID=AEMK02000502.1,length=23675>
##contig=<ID=AEMK02000155.1,length=23493>
##contig=<ID=AEMK02000214.1,length=23279>
##contig=<ID=AEMK02000424.1,length=23240>
##contig=<ID=AEMK02000351.1,length=23110>
##contig=<ID=AEMK02000384.1,length=22885>
##contig=<ID=AEMK02000458.1,length=22707>
##contig=<ID=AEMK02000492.1,length=22673>
##contig=<ID=AEMK02000686.1,length=22650>
##contig=<ID=AEMK02000357.1,length=22457>
##contig=<ID=AEMK02000608.1,length=22393>
##contig=<ID=AEMK02000544.1,length=22327>
##contig=<ID=AEMK02000358.1,length=22269>
##contig=<ID=AEMK02000349.1,length=22254>
##contig=<ID=AEMK02000217.1,length=22216>
##contig=<ID=AEMK02000401.1,length=22208>
##contig=<ID=AEMK02000486.1,length=21821>
##contig=<ID=AEMK02000690.1,length=21690>
##contig=<ID=AEMK02000406.1,length=21463>
##contig=<ID=AEMK02000206.1,length=21439>
##contig=<ID=AEMK02000397.1,length=21409>
##contig=<ID=AEMK02000539.1,length=21264>
##contig=<ID=AEMK02000236.1,length=21107>
##contig=<ID=AEMK02000190.1,length=20938>
##contig=<ID=AEMK02000230.1,length=20938>
##contig=<ID=AEMK02000508.1,length=20819>
##contig=<ID=AEMK02000307.1,length=20771>
##contig=<ID=AEMK02000551.1,length=20589>
##contig=<ID=AEMK02000564.1,length=20419>
##contig=<ID=AEMK02000594.1,length=20300>
##contig=<ID=AEMK02000430.1,length=20230>
##contig=<ID=AEMK02000219.1,length=19883>
##contig=<ID=AEMK02000481.1,length=19728>
##contig=<ID=AEMK02000148.1,length=19496>
##contig=<ID=AEMK02000285.1,length=19482>
##contig=<ID=AEMK02000491.1,length=19466>
##contig=<ID=AEMK02000294.1,length=19443>
##contig=<ID=AEMK02000343.1,length=19362>
##contig=<ID=AEMK02000562.1,length=19341>
##contig=<ID=AEMK02000185.1,length=19277>
##contig=<ID=AEMK02000284.1,length=19274>
##contig=<ID=AEMK02000362.1,length=19234>
##contig=<ID=AEMK02000184.1,length=19005>
##contig=<ID=AEMK02000281.1,length=18788>
##contig=<ID=AEMK02000705.1,length=18345>
##contig=<ID=AEMK02000220.1,length=18189>
##contig=<ID=AEMK02000151.1,length=18134>
##contig=<ID=AEMK02000381.1,length=18081>
##contig=<ID=AEMK02000441.1,length=17962>
##contig=<ID=AEMK02000350.1,length=17924>
##contig=<ID=AEMK02000701.1,length=17846>
##contig=<ID=AEMK02000580.1,length=17747>
##contig=<ID=AEMK02000531.1,length=17708>
##contig=<ID=AEMK02000535.1,length=17699>
##contig=<ID=AEMK02000272.1,length=17331>
##contig=<ID=AEMK02000239.1,length=17318>
##contig=<ID=AEMK02000642.1,length=17027>
##contig=<ID=AEMK02000432.1,length=17025>
##contig=<ID=MT,length=16613>
##contig=<ID=AEMK02000409.1,length=16500>
##contig=<ID=AEMK02000586.1,length=16304>
##contig=<ID=AEMK02000695.1,length=16266>
##contig=<ID=AEMK02000342.1,length=15919>
##contig=<ID=AEMK02000124.1,length=15784>
##contig=<ID=AEMK02000689.1,length=15689>
##contig=<ID=AEMK02000573.1,length=15624>
##contig=<ID=AEMK02000663.1,length=15547>
##contig=<ID=AEMK02000302.1,length=15543>
##contig=<ID=AEMK02000270.1,length=15496>
##contig=<ID=AEMK02000323.1,length=15455>
##contig=<ID=AEMK02000417.1,length=15107>
##contig=<ID=AEMK02000654.1,length=15107>
##contig=<ID=AEMK02000519.1,length=15096>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}
""".format(SAMPLE=sample)
    output.write(header)

def toVcfBreakend(localChr, localPos, localPositive, remoteChr, remotePos, remotePositive):
    if remotePositive:
        remote = "]" + remoteChr + ":" + str(remotePos) + "]"
    else:
        remote = "[" + remoteChr + ":" + str(remotePos) + "["
    return "N" + remote if localPositive else remote + "N"

def isRefToBreakend(orientation):
    parts = orientation.replace("-", "+").split("+")
    return int(parts[0]) >= int(parts[1])

def main(args):
    with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
        write_vcf_header(fout, args.sample)

        line_count = 0
        for line in fin:
            line_count += 1
            if line.startswith('#'): continue
            
            fields = line.strip().split('\t')
            if len(fields) < 10:
                print(f"Skipping malformed line (L{line_count}): {line.strip()}", file=sys.stderr)
                continue
            
            try:
                # 解析核心字段
                leftChr = fields[0]
                leftPos = int(fields[1])
                rightPos = int(fields[4])
                sv_type = fields[6]
                size = int(fields[7])
                score = fields[8]
                numreads = fields[9]
                event_id = f"line{line_count}"

                # 强制修正END值
                if sv_type == "INV":
                    end_value = max(rightPos, leftPos + 1)
                else:
                    end_value = leftPos + abs(size)
                end_value = max(end_value, leftPos + 1)  # 确保END > POS

                # 构建基础INFO字段
                base_info = [
                    f"SVTYPE={sv_type if sv_type not in ('ITX', 'CTX') else 'BND'}",
                    f"END={end_value}",
                    f"SVLEN={-size if sv_type in ('DEL','INS') else size}",
                    f"Chr1={fields[0]}",
                    f"Pos1={fields[1]}",
                    f"Orient1={fields[2]}",
                    f"Chr2={fields[3]}",
                    f"Pos2={fields[4]}", 
                    f"Orient2={fields[5]}",
                    f"Size={size}",
                    f"Score={score}",
                    f"num_Reads={numreads}"
                ]

                # 生成标准记录模板
                base_record = [
                    leftChr,               # CHROM
                    str(leftPos),          # POS
                    event_id,              # ID
                    "N",                   # REF
                    "",                    # ALT (根据类型填充)
                    score,                 # QUAL
                    "PASS",                # FILTER
                    "",                    # INFO (根据类型填充)
                    "GT",                  # FORMAT
                    "0/1"                  # SAMPLE
                ]

                # 处理不同类型
                if sv_type in ("ITX", "CTX"):
                    # 将 sv_type 转换为 BND
                    sv_type = "BND"  # 修改 sv_type 为 BND
                    alt1 = toVcfBreakend(leftChr, leftPos, isRefToBreakend(fields[2]),
                                          fields[3], rightPos, isRefToBreakend(fields[5]))
                    alt2 = toVcfBreakend(fields[3], rightPos, isRefToBreakend(fields[5]),
                                          leftChr, leftPos, isRefToBreakend(fields[2]))
                    
                    for idx, alt in enumerate([alt1, alt2], 1):
                        record = base_record.copy()
                        record[4] = alt
                        record[2] = f"{event_id}_bnd{idx}"
                        record[7] = ";".join(base_info + ["IMPRECISE", "SVTYPE=BND"])  # 确保类型为 BND
                        fout.write("\t".join(record) + "\n")
                else:
                    # 标准SV类型
                    record = base_record.copy()
                    record[4] = f"<{sv_type}>"
                    record[7] = ";".join(base_info + ["IMPRECISE"])
                    fout.write("\t".join(record) + "\n")

            except Exception as e:
                print(f"Error processing line {line_count}: {str(e)}", file=sys.stderr)
                continue

if __name__ == '__main__':
    args = parse_args()
    main(args)