from rdkit import Chem
import numpy as np
import re

def sta_smi(smi):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))
def sta_list_smi(list_smi):
    return list(map(sta_smi,list_smi))

#分子读取与操作
class Molecule:
    def __init__(self,smi,main_ring_smi):
        self.mol = Chem.MolFromSmiles(sta_smi(smi))
        self.pat_main_ring = Chem.MolFromSmarts(sta_smi(main_ring_smi))
        self.a2m,self.m2a = self.get_aro_ring_num(self.mol,self.pat_main_ring)
        self.M_nei = self.get_M_nei(main_ring_smi,self.a2m)
        self.node = self.get_node(self.M_nei)#芳环节点
        self.dic_sub_gro = self.get_dic_sub_gro(self.mol,self.a2m,self.m2a)
        self.atom_feature = self.get_atom_feature(self.mol,self.a2m,self.m2a)


    def get_aro_ring_num(self,mol,pat_main_ring):
        sub_gro = ''
        atom_matches = mol.GetSubstructMatches(pat_main_ring)
        if len(atom_matches)>1:
            pat_check = Chem.MolFromSmarts('*C(=O)O')
            check_matches = self.mol.GetSubstructMatches(pat_check)
            check_matches = list(map(int,np.concatenate(check_matches)))
            for i in range(len(atom_matches)):
                if set(check_matches)&set(atom_matches[i])!=set():
                    a2m = list(atom_matches[i])
                    break
        else: 
            #列表（芳环编号-原编号）
            a2m = list(map(int,np.concatenate(atom_matches)))
        #print(a2m)
        #字典（原编号-芳环编号）
        m2a = {}
        
        for atom in self.mol.GetAtoms():
            if atom.GetIdx() in a2m:
                a_num = a2m.index(atom.GetIdx())
                m2a[atom.GetIdx()] = a_num
                atom.SetProp("atomNote",str(a_num))

        return a2m,m2a
    
    #邻接矩阵
    def get_M_nei(self,main_ring_smi,a2m):
        M_nei = []
        tmp_mol = Chem.MolFromSmiles(sta_smi(main_ring_smi))
        for atom in tmp_mol.GetAtoms():
            tmp_nei = [0]*len(a2m)
            for j in atom.GetNeighbors():
                tmp_nei[j.GetIdx()] = 1
            M_nei.append(tmp_nei)
        return M_nei

    #节点判定
    def get_node(self,M_nei):
        node = [0]*len(M_nei)
        #节点判定
        for i in range(len(M_nei)):
            c = np.sum(np.array(M_nei[i]))
            if c!=2:
                node[i] = 1
        return node
    
    def get_dic_sub_gro(self,mol,a2m,m2a):
        sub_site = []
        sub_gro_site = []
        sub_num = []
        dic_sub_gro = {}
        def get_m2a(m):
            return m2a[m]
        
        for atom in mol.GetAtoms():
            if atom.GetIdx() in a2m:
                for atom_nei in atom.GetNeighbors():
                    if atom_nei.GetIdx() not in a2m:
                        sub_site.append(m2a[atom.GetIdx()])
                        sub_gro_site.append(atom_nei.GetIdx())
                        sub_num.append([atom.GetIdx(),atom_nei.GetIdx()])
        bonds_id = [mol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in sub_num]
        if bonds_id==[]:
            list_frags = []
        else:
            list_frags = Chem.MolToSmiles(Chem.FragmentOnBonds(mol,bonds_id
                        )).split('.')
        #print(sub_gro_site,sub_site,sub_num)
        #print(list_frags)
        for tmp_frag in list_frags:
            if tmp_frag[0]=='*':
                tmp_frag = '[0*]'+tmp_frag[1:]
            tmp_list_num = list((map(int,''.join(re.findall('\d+\*',
            tmp_frag)).split('*')[:-1])))
            new_frag = re.sub('\[\d+\*\]','*',tmp_frag)
            if set(tmp_list_num)!=set(sub_gro_site):
                if new_frag not in dic_sub_gro.keys():
                    dic_sub_gro[new_frag] = [list(map(get_m2a,tmp_list_num))]
                else:
                    dic_sub_gro[new_frag].append(list(map(get_m2a,tmp_list_num)))               
        return dic_sub_gro
    
    def get_atom_feature(self,mol,a2m,m2a):
        dic_atom_feature = {}
        for atom in mol.GetAtoms():
            if atom.GetIdx() in a2m:
                if atom.GetAtomicNum()==6:
                    tmp_AtomicNum = [1,0,0]
                elif atom.GetAtomicNum()==8:
                    tmp_AtomicNum = [0,1,0]
                elif atom.GetAtomicNum()==15:
                    tmp_AtomicNum = [0,0,1]
                if atom.GetHybridization()==Chem.HybridizationType.SP3:
                    tmp_Hybrid = [1,0,0]
                elif atom.GetHybridization()==Chem.HybridizationType.SP2:
                    tmp_Hybrid = [0,1,0]
                elif atom.GetHybridization()==Chem.HybridizationType.SP:
                    tmp_Hybrid = [0,0,1]
                dic_atom_feature[m2a[atom.GetIdx()]] = tmp_AtomicNum+\
                tmp_Hybrid+[atom.GetNumImplicitHs(),
                atom.GetDegree(),atom.GetIsAromatic()]
        return dic_atom_feature
    
#组合分子路径操作
def get_onestart_path(n_start,M_nei,dic_sub_gro):
    onestart_path = []
    for sub in list(dic_sub_gro.keys()):
        for i in dic_sub_gro[sub]:
            for j in i:
                #终点
                n_end = j
                onestart_path.append(get_onestart_end_path(n_start,
                                                           M_nei,n_end))
    return onestart_path
#分子路径操作
def get_onestart_end_path(n_start,M_nei,n_end):
    #上一步结束时已完成与未完成的路径
    l_path = [[n_start]]
    #暂存已完成的路径，暂存未完成的路径的下一步
    tmp_l_path = []
    #走完步数等于节点个数后，只有可能路径遇到重复原子，或走完路径
    for i in range(len(M_nei)):
        #对一条未完结路径
        for j in l_path:
            #取上一个末尾
            end_path = j[-1]
            if end_path==n_end:
                #发现一条路径已完结，保存至tmp_l_path
                tmp_l_path.append(j)
            else:
                #取未经过的相邻原子，走一步；如果没有，则该路径消失
                for k in range(len(M_nei[end_path])):
                    if M_nei[end_path][k]==1 and k not in j:
                        tmp_path = j[:]
                        #路径增长
                        tmp_path.append(k)
                        #增长后路径保存至d2
                        tmp_l_path.append(tmp_path)
        #路径更新
        l_path = tmp_l_path[:]
        #暂存清空
        tmp_l_path = []
    return l_path

#格式化路线
def get_formed_path(path,atom_feature):
    formed_path = []
    for tmp_path in path:
        tmp_formed_path = []
        for i in tmp_path:
            tmp_formed_path.append(atom_feature[i])
        formed_path.append(tmp_formed_path)
    return formed_path