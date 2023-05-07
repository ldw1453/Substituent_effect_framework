import numpy as np
from rdkit import Chem
from sklearn.preprocessing import OneHotEncoder
def get_oh(tem_smiles):    
    oh_enc = OneHotEncoder()
    tem_oh_array=np.array([i if isinstance(i,str) else '' for i in tem_smiles]).reshape(-1,1)
    tem_oh_desc=oh_enc.fit_transform(tem_oh_array).toarray()
    return tem_oh_desc 

def des_std(des_array):
    react_feat_all=des_array[:,des_array.max(axis=0) != des_array.min(axis=0)]
    react_feat_all=(react_feat_all-react_feat_all.min(axis=0))/(react_feat_all.min(axis=0)-react_feat_all.max(axis=0))
    return react_feat_all

