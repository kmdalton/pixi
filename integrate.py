import experiment

nxds_params = {
    'UNIT_CELL_CONSTANTS=' : ' 50 60 70 90 90 90',
    'SPACE_GROUP_NUMBER='  : ' 19',
}

reference = 'off1' #master dataset to use for indexing

denominator = [
    'off1',
    'off2',
    'off3',
    'off4',
]

numerator = [
    'on1',
    'on2',
    'on3',
    'on4',
    'on5',
    'on5',
    'on7',
    'on8',
]

data = {
    'off1' : '../off/1_00001.tiff', 
    'off2' : '../off/2_00001.tiff', 
    'off3' : '../off/3_00001.tiff', 
    'off4' : '../off/4_00001.tiff', 
    'on1'  : '../pumped1/1_00001.tiff', 
    'on2'  : '../pumped1/2_00001.tiff', 
    'on3'  : '../pumped1/3_00001.tiff', 
    'on4'  : '../pumped1/4_00001.tiff', 
    'on5'  : '../pumped1/5_00001.tiff', 
    'on6'  : '../pumped1/6_00001.tiff', 
    'on7'  : '../pumped1/7_00001.tiff', 
    'on8'  : '../pumped1/8_00001.tiff', 
}

scale = {
    'off1' : '../off/metadata.tsv', 
    'off2' : '../off/metadata.tsv', 
    'off3' : '../off/metadata.tsv', 
    'off4' : '../off/metadata.tsv', 
    'on1'  : '../pumped1/metadata.tsv', 
    'on2'  : '../pumped1/metadata.tsv', 
    'on3'  : '../pumped1/metadata.tsv', 
    'on4'  : '../pumped1/metadata.tsv', 
    'on5'  : '../pumped1/metadata.tsv', 
    'on6'  : '../pumped1/metadata.tsv', 
    'on7'  : '../pumped1/metadata.tsv', 
    'on8'  : '../pumped1/metadata.tsv', 
}

e = experiment.experiment()

c = experiment.crystal()

for k,v in data.items():
    scaledict = {i.split()[0]: float(i.split()[1]) for i in open(scale[k]).readlines()[1:]}
    c[k] = experiment.image_series(v).scale(scaledict)
    if k == reference:
        n = c[k].generate_nxdsin()
n.update(nxds_params)

e.append(c)
e.integrate(reference, n)

g = e.doeke_scaling(numerator, denominator)
#print g
