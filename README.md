
# 本科毕业设计：芳环体系中取代基效应的描述符设计与应用

## 环境安装

```bash
python 环境
  - python==3.9.16
  - ipykernel==6.21.2
  - pandas==1.5.2
  - numpy==1.22.3
  - scikit-learn==1.0.2
  - py-xgboost==1.7.3
  - scipy==1.8.1
  - seaborn==0.12.2
  - matplotlib==3.6.2
  - morfeus==0.7.2
  - xtb-python==22.1
```

## 数据集

本工作用到的手性磷酸催化下的亚胺不对称加成反应[数据集](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/dataset/danmark.csv)在[dataset](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/dataset)文件夹中，涉及5种亚胺底物、5种硫醇底物及43种手性磷酸催化剂，共1075条反应数据，具体反应细节参见[原文献](https://www.science.org/doi/10.1126/science.aau5631)。

## 分子信息处理脚本

分子信息处理脚本为[Mol_script.py](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/Mol_script.py)，脚本提供“Molecule”类、“get_onestart_end_path”函数、“get_onestart_path”函数、“get_formed_path”函数接口用以提取分子中的取代基并生成相对取代位置路径，详见脚本代码。

## 描述符生成

取代基框架描述符的生成脚本见[Descriptor_generation.ipynb](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/Descriptor_generation.ipynb)，该脚本给出了中描述符生成的具体代码，生成的描述符均保存在[descriptor](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/descriptor)中。其中[des_danmark_path_onehot.csv](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/descriptor/des_danmark_path_onehot.csv)为One_Hot编码取代基框架描述符而[des_danmark_ste_ele.csv](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/descriptor/des_danmark_ste_ele.csv)为取代基立体电子效应编码取代基框架描述符。

## 描述符应用

取代基框架描述符的评估脚本见[Descriptor_application.ipynb](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/Descriptor_application.ipynb)，该脚本测试了取代基框架描述符的有效性与可解释性，前向特征筛选的[结果](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/result/ste_ele_des_sel.csv)保存在[result](http://10.72.200.35:13000/HongGroup/Substituent_effect_framework/src/branch/main/result)文件夹中。
