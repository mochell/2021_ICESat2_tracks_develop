#%matplotlib inline

execfile(os.environ['PYTHONSTARTUP'])
execfile(STARTUP_2019_surfacemom)

import m_colormanager as m_col

mconfig['paths']
path=mconfig['paths']['config']
A = m_col.ase_to_json(path+'color_def.ase')

B=dict()
for i in A[0]['swatches']:
    B[i['name']] = i['data']['values']
    print(i['name'] + '  ' + str(i['data']['values']))

rels=dict()


rels['plus']=B['red']
rels['minus']=B['blue']

rels['blue']=B['blue']
rels['lightblue']=B['cascade3']
rels['darkblue']=B['cascade1']


rels['white']=B['white']
rels['gridcolor']=B['gridcolor']
rels['grey']=B['gray']

rels['orange']=B['orange']
rels['red']=B['red']
rels['green']=B['green']

rels['cascade1']=B['cascade1']
rels['cascade2']=B['cascade2']
rels['cascade3']=B['cascade3']
rels['cascade4']=B['gridcolor']

rels['rascade1']=B['rascade2']
rels['rascade2']=B['rascade1']
rels['rascade3']=B['rascade3']

rels['aug1']=B['orange']
rels['aug2']=B['green']

B['rels']=rels

m_col.json_save('color_def', path, B)
