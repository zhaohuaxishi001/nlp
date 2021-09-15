import numpy as np
import pandas as pd
import random
import lightgbm
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler, MinMaxScaler

from pandas.core.frame import DataFrame
def cou(table, key, dict, sep=""):
	if not isinstance(key, list):
		key = [key]
	dict = table.groupby(key).aggregate(dict)
	# df = pd.DataFrame(pd.Series(dict), columns=[' delta_g '])
	# df.reset_index().rename(columns={'index':' delta_g '})

	#dict.columns = ["%s%s之%s%s" % (sep, "".join(key), k, A if isinstance(A, str) else A) for k, value in dict.items() for A in (value if isinstance(value, list) else [value])]
	return dict
# def count(table, key, dict, sep=""):
# 	if not isinstance(key, list):

# 		key = [key]
# 	newtable = table.groupby(key).aggregate(dict)
# 	newtable.columns = ["%s%s之%s%s" % (sep, "".join(key), k, a if isinstance(a, str) else a.__name__) for k, value in dict.items() for a in (value if isinstance(value, list) else [value])]
#     return newtable
train=pd.read_csv("final_dataset_train.tsv", sep="\t")
train["id"] = range(-1, -len(train) -1, -1)
# print(train)
test = pd.read_csv("final_dataset_testA.tsv", sep="\t")

test["delta_g"] = -1
# print(test)
train_feature= train.loc[:, ["id", "antibody_seq_a", "antibody_seq_b", "antigen_seq"]]
# print(train_feature)
for x in ["antibody_seq_a", "antibody_seq_b", "antigen_seq"]:
    train_feature["%s長度" % x] = train_feature[x].str.len()
  
    for y in [chr(65 + zimu) for zimu in range(26)]:
        train_feature["%s_%s" % (x, y)] = train_feature[x].str.count(y)
#         for z in [chr(65 + zimu) for zimu in range(26)]:
#             train_feature["%s_%s" % (x, y + z)] = train_feature[x].str.count(y + z)
#             for a in [chr(65 + zimu) for zimu in range(26)]:
#                 train_feature["%s_%s" % (x, y + z + a)] = train_feature[x].str.count(y + z + a)
#                 for b in [chr(65 + zimu) for zimu in range(26)]:
#                      train_feature["%s_%s" % (x, y + z + a + b)] = train_feature[x].str.count(y + z + a+b)
                    
              
train_feature = train_feature.drop(["antibody_seq_a", "antibody_seq_b", "antigen_seq"], axis=1)




test_feature = test.loc[:, ["id", "antibody_seq_a", "antibody_seq_b", "antigen_seq"]]
for x in ["antibody_seq_a", "antibody_seq_b", "antigen_seq"]:
	test_feature["%s長度" % x] = test_feature[x].str.len()
	for y in [chr(65 + zimu) for zimu in range(26)]:
		test_feature["%s_%s" % (x, y)] = test_feature[x].str.count(y)
# 		for z in [chr(65 + zimu) for zimu in range(26)]:
# 			test_feature["%s_%s" % (x, y + z)] = test_feature[x].str.count(y + z)
# 			for a in [chr(65 + zimu) for zimu in range(26)]:
# 				test_feature["%s_%s" % (x, y + z + a)] = test_feature[x].str.count(y + z + a)
# 				for b in [chr(65+i) for i in range(26)]:
# 					test_feature["%s_%s" %(x,y+z+a+ b)]=test_feature[x].str.count(y + z + a+b)
                 
                       
                   
test_feature = test_feature.drop(["antibody_seq_a", "antibody_seq_b", "antigen_seq"], axis=1)

def get(table, basetable,featuretable):
	newtable = table
	newtable = newtable.merge(basetable, on="id", how="left")
	a=cou(featuretable,"antibody_seq_a",{"delta_g": ["mean", "median", "min", "max"]})
	b=cou(featuretable,"antibody_seq_b",{"delta_g": ["mean", "median", "min", "max"]})
	c=cou(featuretable,"antigen_seq",{"delta_g": ["mean", "median", "min", "max"]})

	newtable=newtable.merge(a.reset_index(drop=False), on="antibody_seq_a",how="left")
	newtable=newtable.merge(b.reset_index(drop=False), on="antibody_seq_b",how="left")
	newtable=newtable.merge(c.reset_index(drop=False), on="antigen_seq",how="left")

	

 	# newtable = newtable.merge(featuretable.cou("antibody_seq_a", "delta_g").reset_index(), on="antibody_seq_a",how="left")
# 	newtable = newtable.merge(featuretable.count("antibody_seq_b", {"delta_g": ["mean", "median", "min", "max"]}).reset_index(), on="antibody_seq_b", how="left")
# 	newtable = newtable.merge(featuretable.count("antigen_seq", {"delta_g": ["mean", "median", "min", "max"]}).reset_index(), on="antigen_seq", how="left")
	newtable = newtable.drop(["pdb", "antibody_seq_a", "antibody_seq_b", "antigen_seq"], axis=1)
	
	newtable["label"] = newtable.delta_g.rank()
	newtable = newtable.loc[:, ["id", "delta_g", "label"] + [a for a in newtable.columns if a not in ["id", "delta_g", "label"]]]
	
	return newtable
# def count(table:DataFrame, key, dict, sep="")-> DataFrame:

#     if not isinstance(key, list):


#         key = [key]

#     newtable=table.groupby(key).aggregate(dict)
#     return table

# def count(table, key, dict, sep=""):
# 	if not isinstance(key, list):
# 		key = [key]
# 	newtable = table.groupby(key).aggregate(dict)
# 	# dict.columns = ["%s%s之%s%s" % (sep, "".join(key), k, A if isinstance(A, str) else A.__name__) for k, value in dict.items() for A in (value if isinstance(value, list) else [value])]
# 	return newtable
        

    
# 	newtable.columns = ["%s%s之%s%s" % (sep, "".join(key), k, a if isinstance(a, str) else a.__name__) for k, value in dict.items() for a in (value if isinstance(value, list) else [value])]
# return newtable
# spl = 6
# index = random.sample(range(len(train)), len(train))
# train_detail = None
# for cishu in range(spl):
# 	label = train.iloc[[a for a in range(len(index)) if a % spl == cishu]].reset_index(drop=True)
# 	fea = train.iloc[[a for a in range(len(index)) if a % spl != cishu]].reset_index(drop=True)
	
# 	cishu資料表 = 取得資料表(label, 訓練基礎特征表, fea)
# 	train_detail = pandas.concat([train_detail, cishu資料表], ignore_index=True)

# 輕模型 = lightgbm.train(train_set=lightgbm.Dataset(train_detail.iloc[:, 3:], label=train_detail.標籤)
# 	, num_boost_round=2048, params={"objective": "regression", "learning_rate": 0.05, "max_depth": 6, "num_leaves": 32, "bagging_fraction": 0.7, "feature_fraction": 0.7, "num_threads": 64, "verbose": -1}
# )
spl = 20
index = random.sample(range(len(train)), len(train))

train_detail = None
for cishu in range(spl):
	lab = train.iloc[[i for i in range(len(index)) if i % spl == cishu]].reset_index(drop=True)
	feat = train.iloc[[i for i in range(len(index)) if i % spl != cishu]].reset_index(drop=True)
	
	cishutable = get(lab, train_feature, feat)
	train_detail = pd.concat([train_detail, cishutable], ignore_index=True)
	print("次数为：",cishu)
param = {'num_leaves': 600,
         'min_data_in_leaf': 30,
         'objective': 'rmse',
         'max_depth': -1,
         'learning_rate': 0.05,
         "min_child_samples": 30,
         "boosting": "gbdt",
         "feature_fraction": 0.9,
         "bagging_freq": 1,
         "bagging_fraction": 0.9,
         "bagging_seed": 12,
         "metric": 'mse',
         "lambda_l2": 0.1,
         'is_unbalance': True,
         "verbosity": -1}

light= lightgbm.train(train_set=lightgbm.Dataset(train_detail.iloc[:, 3:], label=train_detail.label)
	, num_boost_round=2048, params=param)
	# 五折交叉验证
# folds = KFold(n_splits=10, shuffle=True, random_state=42)
# oof = np.zeros(len(train))
# predictions = np.zeros(len(test))
# X_train=train_detail.iloc[:, 3:]
# y_train=train_detail.label
# for fold_, (trn_idx, val_idx) in enumerate(folds.split(X_train, y_train)):
#     print("fold n°{}".format(fold_ + 1))
#     trn_data = lightgbm.Dataset(X_train[trn_idx], y_train[trn_idx])
#     val_data = lightgbm.Dataset(X_train[val_idx], y_train[val_idx])

    # num_round = 100000
    # clf = lightgbm.train(param,
    #                 trn_data,
    #                 num_round,
    #                 valid_sets=[trn_data, val_data],
    #                 verbose_eval=2000,
    #                 early_stopping_rounds=1000)
    # oof[val_idx] = clf.predict(X_train[val_idx], num_iteration=clf.best_iteration)
    # b = [round(i, 3) for i in oof]
    # predictions += clf.predict(X_test, num_iteration=clf.best_iteration) / folds.n_splits


out = get(test, test_feature, train)

predict = out.loc[:, ["id"]]
predict["delta_g"] = light.predict(out.iloc[:, 3:])

predict.to_csv("result4.csv", index=False)
